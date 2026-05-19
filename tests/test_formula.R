# Tests for the INLA formula construction in lib.R.
# Stubs inla.group so the test runs without INLA installed.

library(testthat)
source(test_path("..", "lib.R"))

# Stub inla.group with identity for test purposes (the real function bins
# numeric values into discrete groups for RW1 priors; identity preserves
# shape, which is enough to assert on formula structure / column presence).
# Inject into the global env so unqualified lookup from generate_lagged_model
# resolves to the stub.
assign("inla.group", function(x, ...) x, envir = globalenv())

minimal_df <- function() {
  data.frame(
    location         = rep(c("A", "B"), each = 6),
    ID_year          = 1,
    ID_spat          = rep(c("A", "B"), each = 6),
    ID_time_cyclic   = 1,
    week             = c(1:6, 1:6),
    Cases            = 1:12,
    E                = 1000,
    rainfall         = 1:12,
    mean_temperature = 13:24
  )
}

minimal_lag_map <- function(covariates = c("rainfall", "mean_temperature")) {
  data.frame(
    location  = rep(c("A", "B"), times = length(covariates)),
    covariate = rep(covariates, each = 2),
    lag       = 1L,
    stringsAsFactors = FALSE
  )
}

test_that("generate_lagged_model produces one shared RW1 smooth per covariate by default", {
  out <- generate_lagged_model(
    df = minimal_df(),
    covariates = c("rainfall", "mean_temperature"),
    lag_map = minimal_lag_map(),
    region_seasonal = FALSE
  )
  rhs <- paste(as.character(out$formula)[3], collapse = " ")
  expect_match(rhs, "f\\(rainfall_lag_grp, *model")
  expect_match(rhs, "f\\(mean_temperature_lag_grp, *model")
  expect_false(grepl("_grp_loc", rhs))
  expect_true("rainfall_lag_grp" %in% names(out$data))
  expect_true("mean_temperature_lag_grp" %in% names(out$data))
  expect_false("rainfall_lag_grp_loc" %in% names(out$data))
})

test_that("generate_lagged_model adds a per-location RW1 deviation when location_specific_effects=TRUE", {
  out <- generate_lagged_model(
    df = minimal_df(),
    covariates = c("rainfall", "mean_temperature"),
    lag_map = minimal_lag_map(),
    region_seasonal = FALSE,
    location_specific_effects = TRUE
  )
  rhs <- paste(as.character(out$formula)[3], collapse = " ")
  # Shared global term is still present.
  expect_match(rhs, "f\\(rainfall_lag_grp, *model")
  expect_match(rhs, "f\\(mean_temperature_lag_grp, *model")
  # Plus a per-location deviation on a distinct column name with replicate=ID_spat.
  expect_match(rhs, "f\\(rainfall_lag_grp_loc,.*replicate *= *ID_spat")
  expect_match(rhs, "f\\(mean_temperature_lag_grp_loc,.*replicate *= *ID_spat")
  expect_true("rainfall_lag_grp_loc" %in% names(out$data))
  expect_true("mean_temperature_lag_grp_loc" %in% names(out$data))
  # The local column holds the same values as the global one.
  expect_equal(out$data$rainfall_lag_grp_loc, out$data$rainfall_lag_grp)
})

test_that("get_nonlinearity_backend returns the registered backends and errors on unknown names", {
  expect_identical(get_nonlinearity_backend("rw1_inla_group"), backend_rw1_inla_group)
  expect_identical(get_nonlinearity_backend("linear"), backend_linear)
  expect_error(get_nonlinearity_backend("does_not_exist"),
               "Unknown nonlinearity backend")
})

test_that("backend_linear emits a single standardised linear term and no inla.group column by default", {
  df <- data.frame(
    location = rep("A", 4), ID_year = 1, ID_spat = "A",
    week = 1:4, Cases = 1:4, E = 1000,
    rainfall = c(10, 20, 30, 40),
    rainfall_lag = c(NA, 10, 20, 30)
  )
  res <- backend_linear(df, "rainfall", location_specific_effects = FALSE)
  expect_equal(res$terms, "rainfall_lag_z")
  expect_false("rainfall_lag_grp" %in% names(res$df))
  expect_true("rainfall_lag_z" %in% names(res$df))
  # NA in source column gets imputed to 0 (standardised mean).
  expect_equal(res$df$rainfall_lag_z[1], 0)
  # Standardised column has approximately zero mean and unit sd over non-NA values.
  z_nonna <- res$df$rainfall_lag_z[-1]
  expect_lt(abs(mean(z_nonna)), 1e-9)
  expect_lt(abs(stats::sd(z_nonna) - 1), 1e-9)
})

test_that("backend_linear emits a per-location random-slope term when location_specific_effects=TRUE", {
  df <- data.frame(
    location = rep(c("A", "B"), each = 2), ID_year = 1, ID_spat = c("A","A","B","B"),
    week = c(1,2,1,2), Cases = 1:4, E = 1000,
    rainfall = c(10, 20, 30, 40),
    rainfall_lag = c(NA, 10, NA, 30)
  )
  res <- backend_linear(df, "rainfall", location_specific_effects = TRUE)
  expect_length(res$terms, 2)
  expect_equal(res$terms[1], "rainfall_lag_z")
  expect_match(res$terms[2], "^f\\(ID_spat_rainfall, rainfall_lag_z, model='iid'\\)$")
  expect_true("ID_spat_rainfall" %in% names(res$df))
})

test_that("generate_lagged_model accepts an alternative nonlinearity backend", {
  out <- generate_lagged_model(
    df = minimal_df(),
    covariates = c("rainfall", "mean_temperature"),
    lag_map = minimal_lag_map(),
    region_seasonal = FALSE,
    nonlinearity = backend_linear
  )
  rhs <- paste(as.character(out$formula)[3], collapse = " ")
  expect_match(rhs, "\\+ *rainfall_lag_z")
  expect_match(rhs, "\\+ *mean_temperature_lag_z")
  expect_false(grepl("rainfall_lag_grp", rhs))
  expect_false(grepl("inla.group", rhs))
})

test_that("generate_lagged_model still adds the region_seasonal term when requested", {
  out <- generate_lagged_model(
    df = minimal_df(),
    covariates = "rainfall",
    lag_map = minimal_lag_map("rainfall"),
    region_seasonal = TRUE,
    location_specific_effects = TRUE
  )
  rhs <- paste(as.character(out$formula)[3], collapse = " ")
  expect_match(rhs, "f\\(ID_time_cyclic2, *model")
})
