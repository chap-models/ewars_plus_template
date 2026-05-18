# ewars_plus_template

A clean chap-compatible re-implementation of the central modelling ideas in
the upstream ewars_Plus model, layered on top of the `ewars_template`
backbone (Bayesian hierarchical NB regression with INLA).

## What this adds on top of ewars_template

- **Per-district adaptive lag selection** for climate covariates over a
  configurable candidate set (default `[7, 10, 12]` for weekly,
  configurable for monthly), using expanding-window cross-validation. The
  per-district CV scores are produced as a side artefact for inspection;
  the final model uses one lag per covariate, aggregated as the mean
  log-score across districts (argmax with smallest-lag tie-break).
- **Single shifted-lag column per covariate, with an INLA-grouped RW1
  smooth** for the exposure-response shape — matching ewars_Plus's
  production formula (`selected_Model_form_rw`). No dlnm crossbasis.

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
