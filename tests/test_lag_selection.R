# Tests for the lag-selection helpers in lib.R.
# Run with: Rscript -e "testthat::test_file('tests/test_lag_selection.R')"
# from the ewars_plus_template directory.

library(testthat)
source(test_path("..", "lib.R"))

test_that("make_expanding_window_folds produces non-overlapping growing windows", {
  folds <- make_expanding_window_folds(40, 3)
  expect_length(folds, 3)
  expect_equal(length(folds[[1]]$train_idx), 10)
  expect_equal(length(folds[[2]]$train_idx), 20)
  expect_equal(length(folds[[3]]$train_idx), 30)
  for (f in folds) {
    expect_false(any(f$test_idx %in% f$train_idx))
    expect_true(min(f$test_idx) > max(f$train_idx))
  }
})

test_that("select_lags_per_district returns one row per (location, covariate, lag) and uses injected score_fn", {
  df <- data.frame(
    location = rep(c("A", "B"), each = 30),
    ID_year  = 1, week = rep(1:30, 2),
    Cases    = rep(10, 60),
    rainfall = runif(60),
    E        = 1000
  )
  fake_score <- function(df, location, covariate, lag, train_idx, test_idx) {
    if (location == "A" && lag == 7)  return(0)
    if (location == "B" && lag == 12) return(0)
    -1
  }
  scores <- select_lags_per_district(df, "rainfall", c(7, 10, 12),
                                     n_folds = 3, score_fn = fake_score)
  expect_equal(nrow(scores), 6)  # 2 locations x 1 covariate x 3 lags
  expect_setequal(names(scores), c("location", "covariate", "lag", "score"))
  expect_equal(scores$score[scores$location == "A" & scores$lag == 7], 0)
  expect_equal(scores$score[scores$location == "B" & scores$lag == 12], 0)
})

test_that("pick_best_lag_per_covariate aggregates by mean across districts and breaks ties to the smallest lag", {
  score_df <- data.frame(
    location  = rep(c("A", "B"), each = 3),
    covariate = "rainfall",
    lag       = rep(c(7, 10, 12), 2),
    score     = c(-1, -2,  0,
                   0, -2, -1)
  )
  best <- pick_best_lag_per_covariate(score_df)
  expect_named(best, "rainfall")
  # Mean: lag7=-0.5, lag10=-2, lag12=-0.5 -> tie between 7 and 12, pick smallest.
  expect_equal(best$rainfall, 7)
})

test_that("add_lagged_columns shifts within each location and produces NA prefix of length k", {
  df <- data.frame(
    location = rep(c("A", "B"), each = 5),
    ID_year  = 1,
    week     = c(1:5, 1:5),
    rainfall = c(10, 20, 30, 40, 50,    # location A
                 100, 200, 300, 400, 500)  # location B
  )
  out <- add_lagged_columns(df, "rainfall", 2)
  expect_true("rainfall_LAG2" %in% names(out))
  out_a <- out[out$location == "A", ]
  out_a <- out_a[order(out_a$week), ]
  expect_equal(out_a$rainfall_LAG2, c(NA, NA, 10, 20, 30))
  out_b <- out[out$location == "B", ]
  out_b <- out_b[order(out_b$week), ]
  expect_equal(out_b$rainfall_LAG2, c(NA, NA, 100, 200, 300))
})

test_that("add_lagged_columns supports multiple covariates with distinct lags", {
  df <- data.frame(
    location = rep("A", 6),
    ID_year  = 1,
    week     = 1:6,
    rainfall = 1:6,
    mean_temperature = c(10, 20, 30, 40, 50, 60)
  )
  out <- add_lagged_columns(df, c("rainfall", "mean_temperature"), c(1, 3))
  out <- out[order(out$week), ]
  expect_equal(out$rainfall_LAG1, c(NA, 1, 2, 3, 4, 5))
  expect_equal(out$mean_temperature_LAG3, c(NA, NA, NA, 10, 20, 30))
})

test_that("pick_best_lag_per_covariate handles two covariates independently", {
  score_df <- data.frame(
    location  = "A",
    covariate = c("rainfall", "rainfall", "rainfall",
                  "mean_temperature", "mean_temperature", "mean_temperature"),
    lag       = c(7, 10, 12, 7, 10, 12),
    score     = c(-2, -1,  0,
                   0, -1, -2)
  )
  best <- pick_best_lag_per_covariate(score_df)
  expect_equal(best$rainfall, 12)
  expect_equal(best$mean_temperature, 7)
})
