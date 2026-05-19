# ewars_plus_template

A clean chap-compatible re-implementation of the central modelling ideas in
the upstream ewars_Plus model, layered on top of the `ewars_template`
backbone (Bayesian hierarchical NB regression with INLA).

## What this adds on top of ewars_template

- **Per-(location, covariate) adaptive lag selection** over a
  configurable candidate set (default `[7, 10, 12]` for weekly,
  configurable for monthly), using expanding-window cross-validation
  per location. Each location ends up with its own selected lag for
  each climate covariate — closer in spirit to ewars_Plus's per-district
  selection. Ties resolved to the smallest lag (more parsimonious).
- **Single shifted-lag column per covariate, with an INLA-grouped RW1
  smooth** for the exposure-response shape — matching ewars_Plus's
  production formula (`selected_Model_form_rw`). No dlnm crossbasis. The
  column carries different lags across locations, but the smooth shape
  is shared, so cross-location pooling on the exposure-response curve
  is preserved.
- **Optional location-specific deviation** (`location_specific_effects:
  true`): adds a per-location RW1 smooth on the same grouped covariate
  (`replicate = ID_spat`), so each location gets a deviation from the
  shared curve. Partial-pooled via a shared precision hyperprior.
  Default `false` keeps the original shared-smooth-only behaviour.

## Covariate + location-specific lags vs. the upstream ewars_Plus

This section unpacks what's the same, what's different, and why.

### Selection step

| Aspect | Upstream ewars_Plus | `ewars_plus_template` |
|---|---|---|
| Candidate columns | `paste0(var, "_LAG", Min_lag:Max_lag)` — single shifted columns per candidate `k` (`Lag_Model_selection_…R:130`) | Same idea: candidate `k` values come from `candidate_lags`; the scorer builds a shifted column per candidate `k` in-memory. |
| Scoring | Joint INLA fit on the full multi-variable formula, scored by DIC (`Sel_Vars`, ~L285), choosing the best joint combination of `(var1_LAG_k1, var2_LAG_k2, …)`. | Per-(location, covariate) **independent** expanding-window CV. Each cell fits a single-covariate NB+INLA on the training fold and scores the next fold by held-out NB log-density. Folds within a (location, covariate) are averaged. |
| Where the lag varies | Per-district. The brittle bit in upstream: different CV folds within one district could pick different `k`, producing columns like `_LAG12 / _LAG10 / _LAG7` that `foreach(.combine = rbind)` then failed to stack. | Per `(location, covariate)`. Folds within a (location, covariate) are aggregated by mean before `argmax` (smallest-lag tie-break), so each (location, covariate) cell ends up with exactly one selected lag — no across-fold divergence. |
| Output | One selected shifted column name per district per variable, e.g. `Selected_lag_Vars = c("rainfall_LAG10", "mean_temperature_LAG7")`. | `data.frame(location, covariate, lag)` — one row per (location, covariate). Cached as `<model>_lags.rds` by `train.R` so the CV runs once per backtest split, not once per predict call. |

### Materialisation in the design matrix

| Aspect | Upstream ewars_Plus | `ewars_plus_template` |
|---|---|---|
| Column naming | The selected lag is encoded in the column name itself (`rainfall_LAG10`). Different districts in the same fit use the same column name only when they happen to share the selected lag. | Column name is `<cov>_lag` and **doesn't encode** the lag. Each row's value is the covariate at the lag that this row's location chose. So location A's `rainfall_lag` rows hold `rainfall[t-2]`, location B's hold `rainfall[t-3]`, and so on — all in one column. |
| Why it matters | Per-district lag heterogeneity meant the upstream had to fit per district or contort the design matrix; the historical brittleness with `_LAG12 / _LAG10 / _LAG7` and `rbind` came from trying to stack frames with divergent column names. | Per-location lag heterogeneity is moved into the **values**, not the column names. The design matrix has consistent shape across locations regardless of how diverse the lag map is. |

### Final INLA formula

Both models use a `Cases ~ negative binomial, offset = log(E)` likelihood with an iid spatial random effect (replicated across years) and a cyclic RW1 seasonal effect. The covariate piece is what's interesting:

**Upstream `selected_Model_form_rw`** (per district, then INLA is called per district):

```
Cases ~ 1
  + f(ID_spat,    model = "iid", replicate = ID_year)
  + f(week,       model = "rw1", cyclic = TRUE, scale.model = TRUE)
  + f(Var1_Inla_group, model = "rw1", scale.model = TRUE)
  + f(Var2_Inla_group, model = "rw1", scale.model = TRUE)
  + …
```

where `VarN_Inla_group = inla.group(Selected_lag_Vars[N])` and each district's INLA call uses its own selected `Selected_lag_Vars`.

**This template's formula** (single joint INLA call across all locations):

```
Cases ~ 1
  + f(ID_spat,        model = "iid", replicate = ID_year)
  + f(ID_time_cyclic, model = "rw1", cyclic = TRUE, scale.model = TRUE)
  + f(<cov>_lag_grp,  model = "rw1", scale.model = TRUE)               # one per covariate
  [+ f(<cov>_lag_grp_loc, model = "rw1", scale.model = TRUE,
       replicate = ID_spat)]                                            # if location_specific_effects = TRUE
```

where `<cov>_lag_grp = inla.group(<cov>_lag)`. With `location_specific_effects = FALSE` (default) the exposure-response is a single shared smooth across all locations. With it `TRUE`, each location gets a partial-pooled deviation from that shared smooth — a hierarchical decomposition where the global RW1 captures the average exposure-response shape and the per-location RW1 captures location-specific deviation. The two precisions are independent hyperparameters, both estimated from the data.

### Equivalence vs. divergence vs. upstream

- **Lag *placement* is the same idea**: each location picks its own lag from a configurable candidate set, by CV.
- **Lag *carrier* differs**: upstream encodes the lag in the column name and fits per district; this template encodes it in the row's value and fits jointly. Result: same per-location heterogeneity, but one INLA call instead of N, with cross-location pooling for free on the smooth's precision hyperparameter.
- **Exposure-response shape differs by default**: upstream's per-district fits give each district its own RW1 smooth. This template defaults to a shared smooth (no per-location deviation). Set `location_specific_effects: true` to get the per-location RW1 back — and the shared smooth then plays the role of an average / prior mean.
- **What is dropped on purpose**: the joint multi-variable DIC-based selection from upstream, the HTTP service, the on-disk session state, the endemic-channel / outbreak-threshold side outputs.

## What this deliberately does not include from ewars_Plus

- No HTTP/plumber service. The MLproject contract uses simple Rscript
  entry points.
- No on-disk session state per district.
- No endemic-channel / outbreak-threshold side outputs. The model returns
  posterior-predictive samples only; chap-core consumes those and
  renders alarms downstream if desired.
- No `foreach(.combine = rbind)` per-fold stacking. `bind_rows` and
  `match()` are used where the upstream model used the brittle
  patterns we patched in CHAP-core's `external_models/ewars_plus_api_patch`.

## Configuration

```yaml
additional_continuous_covariates:
  - rainfall
  - mean_temperature
user_option_values:
  candidate_lags: [7, 10, 12]
  lag_selection_cv_folds: 3
  precision: 1
  region_seasonal: false
```

Set `n_lags: [7, 10]` (or `n_lags: 7`) to bypass selection and use a
manual lag per covariate.

## Layout

```
MLproject              entry points and adapters
example_config.yaml    default user options
train.R                resolves lags (manual override or CV) and writes
                       `<model>_lags.rds` as a companion file
predict.R              reads the cached lags, fits INLA, samples
lib.R                  period-offset helpers + lag-selection helpers
tests/                 testthat unit tests for the helpers
```

Lag resolution priority in `predict.R`: (1) cached `<model>_lags.rds`
written by `train.R`; (2) `n_lags` manual override from config;
(3) in-predict CV fallback. The cache means the CV loop runs once per
backtest split rather than once per predict call.

## Running the tests

```sh
Rscript -e "testthat::test_dir('tests')"
```

The lag-selection tests inject a stub scoring function and do not
require INLA. The end-to-end smoke test does require INLA + dlnm and
runs against `example_data_monthly`.
