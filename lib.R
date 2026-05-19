### Period-offset helpers (shared with ewars_template) ########################

get_last_month <- function(df) {
  df <- df[!is.na(df$Cases), ]
  df$month[length(df$month)]
}

get_last_week <- function(df) {
  df <- df[!is.na(df$Cases), ]
  df$week[length(df$week)]
}

get_week_diff <- function(df) {
  last_week <- get_last_week(df)
  if (last_week <= 26) 26 - last_week else 78 - last_week
}

offset_years_and_weeks <- function(df) {
  week_diff <- get_week_diff(df)
  new_week <- df$week + week_diff
  df$week <- ((new_week - 1) %% 52) + 1
  df$ID_year <- ifelse(new_week > 52, df$ID_year + 1, df$ID_year)
  df
}

get_month_diff <- function(df) {
  last_month <- get_last_month(df)
  if (last_month <= 6) 6 - last_month else 18 - last_month
}

offset_years_and_months <- function(df) {
  month_diff <- get_month_diff(df)
  new_month <- df$month + month_diff
  df$month <- ((new_month - 1) %% 12) + 1
  df$ID_year <- ifelse(new_month > 12, df$ID_year + 1, df$ID_year)
  df
}

`%||%` <- function(a, b) if (is.null(a)) b else a

### Single shifted-lag column construction ####################################
#
# Match ewars_Plus's production formula: a single shifted column per covariate,
# rather than a dlnm crossbasis spread over 1..K. Each location's rows are
# shifted by *its own* selected lag, so the column carries values from
# different lags across locations. The shared INLA-grouped RW1 smooth on this
# column models a common exposure-response shape; the per-location lag map
# determines which past observation enters each row.

lagged_col_name <- function(covariate) paste0(covariate, "_lag")

# `lag_map`: data.frame(location, covariate, lag) — one row per
# (location, covariate) pair.
add_lagged_columns <- function(df, covariates, lag_map) {
  stopifnot(all(c("location", "covariate", "lag") %in% names(lag_map)))

  time_col <- if ("week" %in% names(df)) "week" else "month"
  df <- df[order(df$location, df$ID_year, df[[time_col]]), , drop = FALSE]

  for (cov in covariates) {
    new_col <- lagged_col_name(cov)
    df[[new_col]] <- NA_real_
    for (loc in unique(df$location)) {
      k_row <- lag_map[lag_map$location == loc & lag_map$covariate == cov, ]
      if (nrow(k_row) != 1L) {
        stop(sprintf("lag_map must have exactly one row for (location=%s, covariate=%s); got %d.",
                     loc, cov, nrow(k_row)))
      }
      k <- as.integer(k_row$lag)
      idx <- which(df$location == loc)
      vals <- df[[cov]][idx]
      shifted <- c(rep(NA_real_, k), utils::head(vals, length(vals) - k))
      df[[new_col]][idx] <- shifted
    }
  }
  df
}

### Lag selection #############################################################
#
# We pick one lag per (location, covariate) using expanding-window CV scores.
# Each location can end up with a different lag for each covariate — closer
# in spirit to ewars_Plus's per-district selection. The design matrix still
# has consistent column shape across locations because each row's lagged
# value enters a shared `<cov>_lag` column.

make_expanding_window_folds <- function(n, n_folds) {
  stopifnot(n_folds >= 1, n >= n_folds + 1)
  fold_size <- n %/% (n_folds + 1)
  lapply(seq_len(n_folds), function(k) {
    list(
      train_idx = seq_len(k * fold_size),
      test_idx  = (k * fold_size + 1):((k + 1) * fold_size)
    )
  })
}

# Default scoring backend: fit a single-covariate NB+INLA on the train slice
# using the single shifted column at lag k, return the sum of NB log-densities
# of held-out rows. Higher is better.
score_inla_holdout <- function(df, location, covariate, lag,
                               train_idx, test_idx) {
  loc_df <- df[df$location == location, , drop = FALSE]
  if (nrow(loc_df) < max(test_idx)) return(NA_real_)
  time_col <- if ("week" %in% names(loc_df)) "week" else "month"
  loc_df <- loc_df[order(loc_df$ID_year, loc_df[[time_col]]), , drop = FALSE]

  vals <- loc_df[[covariate]]
  lagged_x <- c(rep(NA_real_, lag), utils::head(vals, length(vals) - lag))

  fit_df <- data.frame(
    Cases    = loc_df$Cases,
    E        = loc_df$E,
    lagged_x = lagged_x
  )
  fit_df$Cases_obs <- fit_df$Cases
  fit_df$Cases[test_idx] <- NA

  model <- tryCatch(
    INLA::inla(
      Cases ~ 1 + lagged_x, data = fit_df, family = "nbinomial",
      offset = log(E),
      control.compute   = list(config = TRUE),
      control.predictor = list(link = 1, compute = TRUE)
    ),
    error = function(e) NULL
  )
  if (is.null(model)) return(NA_real_)

  mu   <- model$summary.fitted.values$mean[test_idx]
  size <- model$summary.hyperpar["size for the nbinomial observations (1/overdispersion)", "mean"]
  obs  <- fit_df$Cases_obs[test_idx]
  sum(stats::dnbinom(obs, size = size, mu = mu, log = TRUE), na.rm = TRUE)
}

select_lags_per_district <- function(df, covariates, candidate_lags,
                                     n_folds = 3, score_fn = score_inla_holdout) {
  stopifnot(length(covariates) >= 1, length(candidate_lags) >= 1)
  locations <- unique(df$location)

  rows <- list()
  for (loc in locations) {
    loc_df <- df[df$location == loc, , drop = FALSE]
    folds <- make_expanding_window_folds(nrow(loc_df), n_folds)
    for (cov in covariates) {
      for (lag in candidate_lags) {
        fold_scores <- vapply(folds, function(f) {
          score_fn(df, location = loc, covariate = cov, lag = lag,
                   train_idx = f$train_idx, test_idx = f$test_idx)
        }, numeric(1))
        rows[[length(rows) + 1L]] <- data.frame(
          location  = loc,
          covariate = cov,
          lag       = lag,
          score     = mean(fold_scores, na.rm = TRUE)
        )
      }
    }
  }
  do.call(rbind, rows)
}

# Picks one lag per (location, covariate) by argmax of CV log-score. Ties
# resolved to the smallest lag (more parsimonious).
# Returns data.frame(location, covariate, lag).
pick_best_lag_per_location_covariate <- function(score_df) {
  out <- list()
  for (loc in unique(score_df$location)) {
    for (cov in unique(score_df$covariate)) {
      sub <- score_df[score_df$location == loc & score_df$covariate == cov, ]
      sub <- sub[order(-sub$score, sub$lag), ]
      out[[length(out) + 1L]] <- data.frame(
        location  = loc,
        covariate = cov,
        lag       = as.integer(sub$lag[1])
      )
    }
  }
  do.call(rbind, out)
}

### Config + lag resolution (shared by train.R and predict.R) ##################

parse_model_configuration <- function(file_path) {
  config <- yaml::yaml.load_file(file_path)
  user_option_values <- if (!is.null(config$user_option_values))
    jsonlite::fromJSON(jsonlite::toJSON(config$user_option_values)) else list()
  additional_continuous_covariates <- if (!is.null(config$additional_continuous_covariates))
    config$additional_continuous_covariates else character()
  list(user_option_values = user_option_values,
       additional_continuous_covariates = additional_continuous_covariates)
}

# Companion file written by train.R and read by predict.R so the CV runs once.
lags_companion_path <- function(model_fn) paste0(model_fn, "_lags.rds")

# Renders a (location, covariate, lag) lag_map compactly for log messages.
format_lag_map <- function(lag_map) {
  parts <- vapply(unique(lag_map$covariate), function(cov) {
    sub <- lag_map[lag_map$covariate == cov, ]
    sub <- sub[order(sub$location), ]
    sprintf("%s={%s}", cov,
            paste(sub$location, sub$lag, sep = ":", collapse = ", "))
  }, character(1))
  paste(parts, collapse = " | ")
}

# Expands a scalar/vector manual `n_lags` into a per-(location, covariate)
# uniform lag map across the locations seen in `historic_df`.
expand_manual_lags <- function(historic_df, covariates, manual) {
  if (length(manual) == 1) manual <- rep(manual, length(covariates))
  stopifnot(length(manual) == length(covariates))
  locations <- unique(historic_df$location)
  cov_lag <- data.frame(covariate = covariates, lag = as.integer(manual),
                        stringsAsFactors = FALSE)
  out <- merge(data.frame(location = locations, stringsAsFactors = FALSE),
               cov_lag, by = character())
  out[order(out$location, out$covariate), c("location", "covariate", "lag")]
}

# Resolves lags in priority order:
#   1. Cached file at `lags_path` (written by train.R).
#   2. Manual override `user_options$n_lags` (fanned out uniformly per location).
#   3. CV selection on `historic_df` (the in-predict fallback).
# Returns data.frame(location, covariate, lag).
resolve_lags <- function(historic_df, covariates, user_options,
                         lags_path = NULL) {
  if (!is.null(lags_path) && file.exists(lags_path)) {
    cached <- readRDS(lags_path)
    has_required_cols <- is.data.frame(cached) &&
      all(c("location", "covariate", "lag") %in% names(cached))
    if (has_required_cols && all(covariates %in% cached$covariate)) {
      result <- cached[cached$covariate %in% covariates, , drop = FALSE]
      message("Loaded selected lags from ", lags_path, ": ",
              format_lag_map(result))
      return(result)
    }
    message("Cached lags at ", lags_path,
            " missing one of [", paste(covariates, collapse = ", "),
            "] or wrong shape; falling back to in-predict resolution.")
  }

  manual <- user_options$n_lags
  if (!is.null(manual) && length(manual) > 0) {
    result <- expand_manual_lags(historic_df, covariates, manual)
    message("Using manual n_lags override (uniform across locations): ",
            paste(covariates,
                  if (length(manual) == 1) rep(manual, length(covariates)) else manual,
                  sep = "=", collapse = ", "))
    return(result)
  }

  candidate_lags <- user_options$candidate_lags %||% c(7, 10, 12)
  n_folds <- user_options$lag_selection_cv_folds %||% 3
  message("Selecting lags from candidates [",
          paste(candidate_lags, collapse = ", "),
          "] with ", n_folds, "-fold expanding-window CV (per-location)...")
  scores <- select_lags_per_district(
    historic_df, covariates, candidate_lags, n_folds = n_folds
  )
  result <- pick_best_lag_per_location_covariate(scores)
  message("Selected lags: ", format_lag_map(result))
  result
}
