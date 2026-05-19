# Tests for resolve_lags's cache + override + CV-fallback priority chain.
# Run with: Rscript -e "testthat::test_file('tests/test_resolve_lags.R')"

library(testthat)
source(test_path("..", "lib.R"))

test_that("resolve_lags reads cached per-location lag_map when the companion file exists", {
  cache <- tempfile(fileext = ".rds")
  saveRDS(data.frame(
    location  = c("A", "A", "B", "B"),
    covariate = rep(c("rainfall", "mean_temperature"), 2),
    lag       = c(5L, 9L, 3L, 11L)
  ), cache)
  result <- resolve_lags(
    historic_df = data.frame(location = c("A", "B")),
    covariates  = c("rainfall", "mean_temperature"),
    user_options = list(n_lags = 99),  # should be ignored — cache wins
    lags_path = cache
  )
  expect_setequal(names(result), c("location", "covariate", "lag"))
  expect_equal(nrow(result), 4)
  expect_equal(result$lag[result$location == "A" & result$covariate == "rainfall"], 5L)
  expect_equal(result$lag[result$location == "B" & result$covariate == "mean_temperature"], 11L)
})

test_that("resolve_lags rejects cached files of the wrong shape and falls back", {
  bad <- tempfile(fileext = ".rds")
  saveRDS(list(rainfall = 5L), bad)  # legacy named-list shape
  result <- resolve_lags(
    historic_df = data.frame(location = c("A")),
    covariates  = "rainfall",
    user_options = list(n_lags = 4),
    lags_path = bad
  )
  expect_equal(nrow(result), 1)
  expect_equal(result$lag, 4L)
})

test_that("resolve_lags expands manual n_lags uniformly across the historic_df's locations", {
  result <- resolve_lags(
    historic_df = data.frame(location = c("A", "A", "B", "B", "C")),
    covariates  = c("rainfall", "mean_temperature"),
    user_options = list(n_lags = c(4, 6)),
    lags_path = tempfile(fileext = ".rds")  # never created
  )
  expect_equal(nrow(result), 6)  # 3 locations × 2 covariates
  expect_setequal(unique(result$location), c("A", "B", "C"))
  expect_equal(result$lag[result$covariate == "rainfall"], rep(4L, 3))
  expect_equal(result$lag[result$covariate == "mean_temperature"], rep(6L, 3))
})

test_that("resolve_lags falls back to CV (per-location) when neither cache nor override is set", {
  fake_score <- function(df, location, covariate, lag, train_idx, test_idx) {
    if (location == "A" && lag == 2) return(0)
    if (location == "B" && lag == 4) return(0)
    -1
  }
  saved <- select_lags_per_district
  select_lags_per_district <<- function(df, covariates, candidate_lags,
                                        n_folds = 3,
                                        score_fn = score_inla_holdout) {
    rows <- list()
    for (loc in unique(df$location)) {
      for (cov in covariates) for (lag in candidate_lags) {
        rows[[length(rows) + 1L]] <- data.frame(
          location = loc, covariate = cov, lag = lag,
          score = fake_score(df, loc, cov, lag, 1L, 2L)
        )
      }
    }
    do.call(rbind, rows)
  }
  on.exit(select_lags_per_district <<- saved, add = TRUE)

  result <- resolve_lags(
    historic_df = data.frame(location = c("A", "B")),
    covariates  = "rainfall",
    user_options = list(candidate_lags = c(2, 3, 4),
                        lag_selection_cv_folds = 1),
    lags_path = NULL
  )
  expect_equal(nrow(result), 2)
  expect_equal(result$lag[result$location == "A"], 2L)
  expect_equal(result$lag[result$location == "B"], 4L)
})
