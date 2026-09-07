# Tests for parent-org-unit lag grouping (lag_grouping: parent).
# Run with: Rscript -e "testthat::test_file('tests/test_lag_grouping.R')"

library(testthat)
source(test_path("..", "lib.R"))

# Two provinces, two districts each, 12 weeks of data per district.
grouped_df <- function() {
  districts <- c("d1", "d2", "d3", "d4")
  parents   <- c("p1", "p1", "p2", "p2")
  do.call(rbind, lapply(seq_along(districts), function(i) {
    data.frame(
      location = districts[i],
      parent   = parents[i],
      ID_year  = 1,
      week     = 1:12,
      Cases    = i,
      E        = i * 1000,
      rainfall = i * 10,
      stringsAsFactors = FALSE
    )
  }))
}

test_that("has_usable_parents rejects a missing column and the '-' placeholder", {
  df <- grouped_df()
  expect_true(has_usable_parents(df))
  expect_false(has_usable_parents(df[, setdiff(names(df), "parent")]))

  placeholder <- df
  placeholder$parent <- "-"
  expect_false(has_usable_parents(placeholder))
})

test_that("location_group_map returns one row per location", {
  map <- location_group_map(grouped_df())
  expect_equal(nrow(map), 4)
  expect_setequal(names(map), c("location", "group"))
  expect_equal(map$group[map$location == "d1"], "p1")
  expect_equal(map$group[map$location == "d4"], "p2")
})

test_that("aggregate_to_parent sums cases and population and keys on the parent id", {
  agg <- aggregate_to_parent(grouped_df(), "rainfall")
  expect_equal(nrow(agg), 24)  # 2 parents x 12 weeks
  expect_setequal(unique(agg$location), c("p1", "p2"))

  p1_w1 <- agg[agg$location == "p1" & agg$week == 1, ]
  expect_equal(p1_w1$Cases, 3)     # d1 (1) + d2 (2)
  expect_equal(p1_w1$E, 3000)      # 1000 + 2000
})

test_that("aggregate_to_parent weights covariates by population, not by district count", {
  agg <- aggregate_to_parent(grouped_df(), "rainfall")
  p1_w1 <- agg[agg$location == "p1" & agg$week == 1, ]
  # d1: rainfall 10 at E 1000, d2: rainfall 20 at E 2000 -> (10*1000 + 20*2000)/3000
  expect_equal(p1_w1$rainfall, 50 / 3)
  # A plain unweighted mean would give 15.
  expect_false(isTRUE(all.equal(p1_w1$rainfall, 15)))
})

test_that("aggregate_to_parent falls back to equal weights when population is zero", {
  df <- grouped_df()
  df$E <- 0
  agg <- aggregate_to_parent(df, "rainfall")
  p1_w1 <- agg[agg$location == "p1" & agg$week == 1, ]
  expect_equal(p1_w1$rainfall, 15)  # unweighted mean of 10 and 20
})

test_that("expand_group_lags fans a group lag out to every location in the group", {
  group_lag_map <- data.frame(
    location  = c("p1", "p2"),  # group ids sit in `location` after selection
    covariate = "rainfall",
    lag       = c(7L, 12L),
    stringsAsFactors = FALSE
  )
  out <- expand_group_lags(group_lag_map, location_group_map(grouped_df()))
  expect_setequal(names(out), c("location", "covariate", "lag"))
  expect_equal(nrow(out), 4)
  expect_equal(out$lag[out$location == "d1"], 7L)
  expect_equal(out$lag[out$location == "d2"], 7L)
  expect_equal(out$lag[out$location == "d3"], 12L)
  expect_equal(out$lag[out$location == "d4"], 12L)
})

test_that("expand_group_lags output satisfies add_lagged_columns' completeness requirement", {
  df <- grouped_df()
  group_lag_map <- data.frame(
    location = c("p1", "p2"), covariate = "rainfall", lag = c(2L, 3L),
    stringsAsFactors = FALSE
  )
  lag_map <- expand_group_lags(group_lag_map, location_group_map(df))
  expect_silent(add_lagged_columns(df, "rainfall", lag_map))
})

test_that("resolve_lag_grouping defaults to location and rejects unknown values", {
  expect_equal(resolve_lag_grouping(grouped_df(), NULL), "location")
  expect_equal(resolve_lag_grouping(grouped_df(), "parent"), "parent")
  expect_error(resolve_lag_grouping(grouped_df(), "province"),
               "Unknown lag_grouping")
})

test_that("resolve_lag_grouping warns and falls back when parents are placeholders", {
  df <- grouped_df()
  df$parent <- "-"
  expect_warning(result <- resolve_lag_grouping(df, "parent"),
                 "no usable")
  expect_equal(result, "location")
})

test_that("resolve_lags with lag_grouping='parent' selects per group and returns per location", {
  df <- grouped_df()
  # Score by group id, so a per-location selection could not produce this answer.
  fake_score <- function(df, location, covariate, lag, train_idx, test_idx) {
    if (location == "p1" && lag == 7)  return(0)
    if (location == "p2" && lag == 12) return(0)
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
          score = fake_score(df, loc, cov, lag, 1L, 2L),
          stringsAsFactors = FALSE
        )
      }
    }
    do.call(rbind, rows)
  }
  on.exit(select_lags_per_district <<- saved, add = TRUE)

  result <- resolve_lags(
    historic_df = df,
    covariates  = "rainfall",
    user_options = list(candidate_lags = c(7, 10, 12),
                        lag_selection_cv_folds = 1,
                        lag_grouping = "parent"),
    lags_path = NULL
  )
  expect_equal(nrow(result), 4)  # one row per district, not per province
  expect_setequal(result$location, c("d1", "d2", "d3", "d4"))
  expect_equal(result$lag[result$location == "d1"], 7L)
  expect_equal(result$lag[result$location == "d3"], 12L)
})
