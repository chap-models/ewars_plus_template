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

test_that("pick_best_lag_per_location_covariate picks argmax per (location, covariate) with smallest-lag tie-break", {
  score_df <- data.frame(
    location  = rep(c("A", "B"), each = 3),
    covariate = "rainfall",
    lag       = rep(c(7, 10, 12), 2),
    score     = c(-1, -2,  0,    # A: lag 12 best
                   0, -2, -1)    # B: lag 7 best
  )
  best <- pick_best_lag_per_location_covariate(score_df)
  expect_setequal(names(best), c("location", "covariate", "lag"))
  expect_equal(nrow(best), 2)
  expect_equal(best$lag[best$location == "A"], 12L)
  expect_equal(best$lag[best$location == "B"], 7L)
})

test_that("pick_best_lag_per_location_covariate breaks ties to the smallest lag per (location, covariate)", {
  score_df <- data.frame(
    location  = c("A", "A", "A"),
    covariate = "rainfall",
    lag       = c(7, 10, 12),
    score     = c(0, -1, 0)  # tie between 7 and 12
  )
  best <- pick_best_lag_per_location_covariate(score_df)
  expect_equal(best$lag[best$location == "A"], 7L)
})

test_that("pick_best_lag_per_location_covariate handles multiple covariates and locations independently", {
  score_df <- expand.grid(
    location  = c("A", "B"),
    covariate = c("rainfall", "mean_temperature"),
    lag       = c(7, 10, 12),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  score_df$score <- 0
  # A/rainfall best at 12, A/temp best at 7, B/rainfall best at 7, B/temp best at 10.
  score_df$score[score_df$location == "A" & score_df$covariate == "rainfall" & score_df$lag == 12] <- 1
  score_df$score[score_df$location == "A" & score_df$covariate == "mean_temperature" & score_df$lag == 7] <- 1
  score_df$score[score_df$location == "B" & score_df$covariate == "rainfall" & score_df$lag == 7] <- 1
  score_df$score[score_df$location == "B" & score_df$covariate == "mean_temperature" & score_df$lag == 10] <- 1
  best <- pick_best_lag_per_location_covariate(score_df)
  expect_equal(nrow(best), 4)
  expect_equal(best$lag[best$location == "A" & best$covariate == "rainfall"], 12L)
  expect_equal(best$lag[best$location == "A" & best$covariate == "mean_temperature"], 7L)
  expect_equal(best$lag[best$location == "B" & best$covariate == "rainfall"], 7L)
  expect_equal(best$lag[best$location == "B" & best$covariate == "mean_temperature"], 10L)
})

test_that("add_lagged_columns shifts within each location using that location's lag and produces NA prefix", {
  df <- data.frame(
    location = rep(c("A", "B"), each = 5),
    ID_year  = 1,
    week     = c(1:5, 1:5),
    rainfall = c(10, 20, 30, 40, 50,        # location A
                 100, 200, 300, 400, 500)   # location B
  )
  lag_map <- data.frame(
    location  = c("A", "B"),
    covariate = "rainfall",
    lag       = c(2L, 3L)
  )
  out <- add_lagged_columns(df, "rainfall", lag_map)
  expect_true("rainfall_lag" %in% names(out))

  out_a <- out[out$location == "A", ]
  out_a <- out_a[order(out_a$week), ]
  expect_equal(out_a$rainfall_lag, c(NA, NA, 10, 20, 30))  # lag = 2

  out_b <- out[out$location == "B", ]
  out_b <- out_b[order(out_b$week), ]
  expect_equal(out_b$rainfall_lag, c(NA, NA, NA, 100, 200))  # lag = 3
})

test_that("add_lagged_columns supports multiple covariates with per-(location, covariate) lags", {
  df <- data.frame(
    location = rep("A", 6),
    ID_year  = 1,
    week     = 1:6,
    rainfall = 1:6,
    mean_temperature = c(10, 20, 30, 40, 50, 60)
  )
  lag_map <- data.frame(
    location  = "A",
    covariate = c("rainfall", "mean_temperature"),
    lag       = c(1L, 3L)
  )
  out <- add_lagged_columns(df, c("rainfall", "mean_temperature"), lag_map)
  out <- out[order(out$week), ]
  expect_equal(out$rainfall_lag, c(NA, 1, 2, 3, 4, 5))
  expect_equal(out$mean_temperature_lag, c(NA, NA, NA, 10, 20, 30))
})

test_that("add_lagged_columns errors clearly when lag_map is missing a (location, covariate)", {
  df <- data.frame(
    location = rep("A", 3),
    ID_year  = 1, week = 1:3,
    rainfall = 1:3
  )
  lag_map <- data.frame(location = "B", covariate = "rainfall", lag = 1L)
  expect_error(
    add_lagged_columns(df, "rainfall", lag_map),
    "lag_map must have exactly one row"
  )
})
