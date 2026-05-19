# Cases = number of cases | E = population | week/month period component
# ID_year = year | ID_spat = location

library(yaml)
library(jsonlite)
library(INLA)
library(dplyr)
source("lib.R")

# Build the production formula matching ewars_Plus's `selected_Model_form_rw`:
# a per-covariate RW1 smooth on the inla.group()'d shifted column, no separate
# linear term (the RW1 captures the exposure-response shape).
generate_lagged_model <- function(df, covariates, lag_map, region_seasonal) {
  df <- add_lagged_columns(df, covariates, lag_map)

  smooth_terms <- character()
  for (cov in covariates) {
    col <- lagged_col_name(cov)
    grp <- paste0(col, "_grp")
    df[[grp]] <- inla.group(df[[col]])
    smooth_terms <- c(
      smooth_terms,
      sprintf("f(%s, model='rw1', scale.model=TRUE)", grp)
    )
  }

  formula_str <- paste(
    "Cases ~ 1 +",
    "f(ID_spat, model='iid', replicate=ID_year) +",
    "f(ID_time_cyclic, model='rw1', cyclic=TRUE, scale.model=TRUE) +",
    paste(smooth_terms, collapse = " + ")
  )
  if (region_seasonal) {
    formula_str <- paste(formula_str,
      "+ f(ID_time_cyclic2, model='rw1', cyclic=TRUE, scale.model=TRUE, replicate=ID_spat)")
  }
  list(formula = as.formula(formula_str), data = df)
}

generate_basic_model <- function(df, region_seasonal) {
  formula_str <- paste(
    "Cases ~ 1 +",
    "f(ID_spat, model='iid', replicate=ID_year) +",
    "f(ID_time_cyclic, model='rw1', cyclic=TRUE, scale.model=TRUE)"
  )
  if (region_seasonal) {
    formula_str <- paste(formula_str,
      "+ f(ID_time_cyclic2, model='rw1', cyclic=TRUE, scale.model=TRUE, replicate=ID_spat)")
  }
  list(formula = as.formula(formula_str), data = df)
}

predict_chap <- function(model_fn, hist_fn, future_fn, preds_fn, config_fn = "") {
  if (config_fn != "") {
    config <- parse_model_configuration(config_fn)
    covariate_names <- config$additional_continuous_covariates
    user_options    <- config$user_option_values
  } else {
    covariate_names <- c("rainfall", "mean_temperature")
    user_options    <- list()
  }
  precision       <- user_options$precision %||% 0.01
  region_seasonal <- user_options$region_seasonal %||% FALSE

  historic_df <- read.csv(hist_fn)
  future_df   <- read.csv(future_fn)
  future_df$Cases <- rep(NA, nrow(future_df))
  future_df$disease_cases <- rep(NA, nrow(future_df))
  df <- rbind(historic_df, future_df)

  if ("week" %in% colnames(df)) {
    df <- mutate(df, ID_time_cyclic = week)
    df <- offset_years_and_weeks(df)
  } else {
    df <- mutate(df, ID_time_cyclic = month)
    df <- offset_years_and_months(df)
  }
  df$ID_time_cyclic2 <- df$ID_time_cyclic
  df$ID_year <- df$ID_year - min(df$ID_year) + 1

  if (length(covariate_names) == 0) {
    generated <- generate_basic_model(df, region_seasonal)
  } else {
    lag_map <- resolve_lags(historic_df, covariate_names, user_options,
                            lags_path = lags_companion_path(model_fn))
    generated <- generate_lagged_model(df, covariate_names, lag_map, region_seasonal)
  }
  formula_used <- generated$formula
  df <- generated$data
  df$ID_spat <- as.integer(as.factor(df$ID_spat))

  model <- inla(
    formula = formula_used, data = df, family = "nbinomial", offset = log(E),
    control.inla     = list(strategy = "adaptive"),
    control.compute  = list(dic = TRUE, config = TRUE, cpo = TRUE, return.marginals = FALSE),
    control.fixed    = list(correlation.matrix = TRUE, prec.intercept = 1e-4, prec = precision),
    control.predictor = list(link = 1, compute = TRUE),
    verbose = FALSE, safe = FALSE
  )

  casestopred <- df$Cases
  idx.pred <- which(is.na(casestopred))
  mpred <- length(idx.pred)
  s <- 1000
  y.pred <- matrix(NA, mpred, s)
  xx <- inla.posterior.sample(s, model)
  xx.s <- inla.posterior.sample.eval(
    function(idx.pred) c(theta[1], Predictor[idx.pred]), xx, idx.pred = idx.pred
  )
  for (s.idx in seq_len(s)) {
    xx.sample <- xx.s[, s.idx]
    y.pred[, s.idx] <- rnbinom(mpred, mu = exp(xx.sample[-1]), size = xx.sample[1])
  }

  new.df <- data.frame(
    time_period = df$time_period[idx.pred],
    location    = df$location[idx.pred],
    y.pred
  )
  colnames(new.df) <- c("time_period", "location", paste0("sample_", 0:(s - 1)))

  write.csv(new.df, preds_fn, row.names = FALSE)
  saveRDS(model, file = model_fn)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  model_fn  <- args[1]
  hist_fn   <- args[2]
  future_fn <- args[3]
  preds_fn  <- args[4]
  config_fn <- if (length(args) == 5) args[5] else ""
  predict_chap(model_fn, hist_fn, future_fn, preds_fn, config_fn)
}
