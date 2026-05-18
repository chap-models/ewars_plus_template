# Cases = number of cases | E = population | week/month period component
# ID_year = year | ID_spat = location
#
# Train computes the per-covariate lag map once (manual override, or CV on the
# training slice) and writes it to `<model>_lags.rds`. predict.R reads it back,
# so the CV loop runs once per backtest split rather than once per predict call.

library(yaml)
library(jsonlite)
library(INLA)
source("lib.R")

train_chap <- function(train_fn, model_fn, config_fn = "") {
  if (config_fn == "") {
    message("No config supplied; nothing to do.")
    return(invisible(NULL))
  }
  config <- parse_model_configuration(config_fn)
  covariate_names <- config$additional_continuous_covariates
  user_options    <- config$user_option_values

  if (length(covariate_names) == 0) {
    message("No additional_continuous_covariates; nothing to select.")
    return(invisible(NULL))
  }

  train_df <- read.csv(train_fn)
  selected_lags <- resolve_lags(train_df, covariate_names, user_options,
                                lags_path = NULL)

  lags_path <- lags_companion_path(model_fn)
  saveRDS(as.list(selected_lags), file = lags_path)
  message("Wrote selected lags to ", lags_path, ".")
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2) {
  train_fn  <- args[1]
  model_fn  <- args[2]
  config_fn <- if (length(args) >= 3) args[3] else ""
  train_chap(train_fn, model_fn, config_fn)
}
