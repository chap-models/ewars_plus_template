# Tests for resolve_lags's cache + override + CV-fallback priority.
# Run with: Rscript -e "testthat::test_file('tests/test_resolve_lags.R')"

library(testthat)
source(test_path("..", "lib.R"))

writeLines(c(
  "additional_continuous_covariates:",
  "  - rainfall",
  "  - mean_temperature",
  "user_option_values: {}"
), con = (config_no_override <- tempfile(fileext = ".yaml")))

test_that("resolve_lags reads cached lags when the companion file exists", {
  cache <- tempfile(fileext = ".rds")
  saveRDS(list(rainfall = 5L, mean_temperature = 9L), cache)
  result <- resolve_lags(
    historic_df = data.frame(),
    covariates  = c("rainfall", "mean_temperature"),
    user_options = list(n_lags = 99),  # should be ignored — cache wins
    lags_path = cache
  )
  expect_equal(unname(result), c(5L, 9L))
  expect_named(result, c("rainfall", "mean_temperature"))
})

test_that("resolve_lags falls back to manual override when cache is missing or incomplete", {
  missing <- tempfile(fileext = ".rds")  # never created
  result <- resolve_lags(
    historic_df = data.frame(),
    covariates  = c("rainfall", "mean_temperature"),
    user_options = list(n_lags = 4),
    lags_path = missing
  )
  expect_equal(unname(result), c(4L, 4L))

  partial <- tempfile(fileext = ".rds")
  saveRDS(list(rainfall = 5L), partial)  # missing mean_temperature
  result2 <- resolve_lags(
    historic_df = data.frame(),
    covariates  = c("rainfall", "mean_temperature"),
    user_options = list(n_lags = 3),
    lags_path = partial
  )
  expect_equal(unname(result2), c(3L, 3L))
})

test_that("resolve_lags falls back to CV when neither cache nor override is set", {
  # Stub the scorer so we don't need INLA in tests.
  fake_score <- function(df, location, covariate, lag, train_idx, test_idx) {
    if (covariate == "rainfall" && lag == 2)  return(0)
    if (covariate == "mean_temperature" && lag == 4) return(0)
    -1
  }
  with_stub <- function() {
    local_score_fn <- fake_score
    select_lags_per_district <<- function(df, covariates, candidate_lags,
                                          n_folds = 3,
                                          score_fn = score_inla_holdout) {
      # Delegate to the real function but with our stub baked in.
      rows <- list()
      for (loc in unique(df$location)) {
        for (cov in covariates) for (lag in candidate_lags) {
          rows[[length(rows) + 1L]] <- data.frame(
            location = loc, covariate = cov, lag = lag,
            score = local_score_fn(df, loc, cov, lag, 1L, 2L)
          )
        }
      }
      do.call(rbind, rows)
    }
  }
  saved <- select_lags_per_district
  on.exit(select_lags_per_district <<- saved, add = TRUE)
  with_stub()

  result <- resolve_lags(
    historic_df = data.frame(location = "A"),
    covariates = c("rainfall", "mean_temperature"),
    user_options = list(candidate_lags = c(2, 3, 4), lag_selection_cv_folds = 1),
    lags_path = NULL
  )
  expect_equal(unname(result), c(2L, 4L))
})
