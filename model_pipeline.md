# How `ewars_plus_template` works

## High-level pipeline

```mermaid
flowchart TD
    subgraph TRAIN["train.R"]
        T_IN["training_data.csv + model_config.yaml"] --> T_RESOLVE{"n_lags<br/>provided?"}
        T_RESOLVE -- "yes" --> T_MANUAL["use n_lags directly"]
        T_RESOLVE -- "no" --> T_CV["select_lags_per_district<br/>(see detail below)"]
        T_MANUAL --> T_WRITE["write &lt;model&gt;_lags.rds"]
        T_CV --> T_WRITE
    end

    subgraph PREDICT["predict.R"]
        H["historic_data.csv"] --> RBIND
        F["future_data.csv<br/>covariates only, Cases = NA"] --> RBIND
        RBIND["rbind hist + future"] --> OFFSET["offset_years_and_weeks /<br/>offset_years_and_months"]

        READ["read &lt;model&gt;_lags.rds<br/>(fallback: manual / in-predict CV)"]
        READ --> BUILD
        OFFSET --> BUILD

        BUILD["add_lagged_columns<br/>cov_LAG_k per location"] --> GROUP["inla.group on each<br/>cov_LAG_k → cov_LAG_k_grp"]
        GROUP --> FORMULA["build formula"]
        FORMULA --> INLA["inla(NB, offset=log E)<br/>strategy=adaptive"]
        INLA --> SAMPLE["inla.posterior.sample × 1000<br/>+ rnbinom(mu, size)"]
        SAMPLE --> OUT["predictions.csv<br/>time_period, location,<br/>sample_0 … sample_999"]
    end

    T_WRITE -. "&lt;model&gt;_lags.rds" .-> READ
```

## Per-district lag selection (CV path only)

```mermaid
flowchart TD
    IN["historic_data per location<br/>(sorted by year, period)"] --> FOLDS["make_expanding_window_folds<br/>K folds: train = 1..k·N/(K+1),<br/>test = next N/(K+1) rows"]

    FOLDS --> LOOP["for each<br/>location × covariate × candidate_lag × fold:"]
    LOOP --> SHIFT["materialize single shifted column<br/>x_LAG_k[t] = x[t - k]"]
    SHIFT --> FIT["INLA NB:<br/>Cases ~ 1 + x_LAG_k<br/>(Cases = NA on test rows)"]
    FIT --> SCORE["sum log dnbinom(obs, mu, size)<br/>on test rows<br/>(higher = better)"]

    SCORE --> AGG_LOC["mean across folds<br/>→ score per (location, covariate, lag)"]
    AGG_LOC --> AGG_GLOBAL["mean across locations,<br/>argmax (smallest-lag tie-break)<br/>→ one lag per covariate"]
    AGG_GLOBAL --> NLAG["nlag vector aligned<br/>with covariates"]
```

## Final INLA formula

```
Cases ~ 1
  + f(ID_spat,        model = "iid",  replicate = ID_year)
  + f(ID_time_cyclic, model = "rw1",  cyclic = TRUE, scale.model = TRUE)
  + f(<cov_i>_LAG<k_i>_grp, model = "rw1", scale.model = TRUE)   # one term per covariate
  [+ f(ID_time_cyclic2, model = "rw1", cyclic = TRUE, scale.model = TRUE,
       replicate = ID_spat)]                                     # if region_seasonal = TRUE
```

Family: `nbinomial`. Offset: `log(E)` (population). The covariate enters
only through the RW1 smooth on its discretised (`inla.group`) shifted
column — no separate linear term, matching ewars_Plus's
`selected_Model_form_rw`.

## Structural correspondence with upstream ewars_Plus

| Concept | ewars_Plus | ewars_plus_template |
|---|---|---|
| Per-district lag candidate columns | `paste0(alarm_vars, "_LAG", Min_lag:Max_lag)` (`Lag_Model_selection_…R:130`) | `add_lagged_columns(df, covariates, lags)` |
| Lag selection criterion | INLA-based variable selection (`Sel_Vars`, ~L285) on the full multi-variable formula | per-(district, covariate) expanding-window CV log-score |
| Aggregation across districts | implicit (one lag fixed per variable for the joint fit) | mean log-score across districts, argmax with smallest-lag tie-break |
| Final formula | `selected_Model_form_rw` — RW1 smooth on `inla.group`'d selected shifted column | same shape: `f(cov_LAG_k_grp, model='rw1', scale.model=TRUE)` |
| Predictive sampling | INLA posterior sample → `rnbinom` | INLA posterior sample → `rnbinom` |
| Output rows | sparse subset of requested forecast weeks | one row per `Cases = NA` row of `rbind(historic, future)` |
| Endemic channel / alarms | computed and emitted | **dropped** |
| Container statefulness | trained RDS state in container | **none** — train.R is a no-op; predict.R retrains each call |
