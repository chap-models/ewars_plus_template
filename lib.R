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

### Parent-org-unit lag grouping ###############################################
#
# chap-core writes a `parent` column — the parent org unit id taken from the
# dataset's geojson feature properties — into the CSVs handed to the model.
# With `lag_grouping: parent`, selection runs once per parent group on an
# aggregated series instead of once per location, and the winning lag is fanned
# back out to every location in the group. For 146 districts across 18
# provinces that is 18 CV sweeps rather than 146.
#
# Beyond the cost, this is also a shrinkage choice: in sparse district-week
# data a per-district argmax is largely fitting noise, so sharing one lag
# across a province is often the better estimate.

parent_placeholder <- "-"

# chap-core fills the `parent` column with "-" when the dataset has no geojson,
# or when the features carry no `parent` property. Grouping on that would
# silently collapse every location into one group, so callers check this first.
has_usable_parents <- function(df) {
  if (!"parent" %in% names(df)) return(FALSE)
  parents <- df$parent
  !all(is.na(parents) | parents == parent_placeholder)
}

# data.frame(location, group) — one row per location, taking the first parent
# seen for it.
location_group_map <- function(df) {
  map <- data.frame(location = df$location, group = df$parent,
                    stringsAsFactors = FALSE)
  map[!duplicated(map$location), , drop = FALSE]
}

# Collapses the per-location frame to one series per parent group. Cases and
# population sum; covariates are population-weighted means, so a province's
# climate series is not dominated by its smallest district. The group id is
# written into `location` so the existing selection code runs against the
# aggregate unchanged.
aggregate_to_parent <- function(df, covariates) {
  time_col <- if ("week" %in% names(df)) "week" else "month"
  weights <- df$E
  weights[!is.finite(weights) | weights < 0] <- 0

  keys <- paste(df$parent, df$ID_year, df[[time_col]], sep = "\r")
  rows <- lapply(split(seq_len(nrow(df)), keys), function(idx) {
    w <- weights[idx]
    if (sum(w) <= 0) w <- rep(1, length(idx))
    row <- data.frame(
      location = df$parent[idx[1]],
      ID_year  = df$ID_year[idx[1]],
      Cases    = sum(df$Cases[idx], na.rm = TRUE),
      E        = sum(df$E[idx], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    row[[time_col]] <- df[[time_col]][idx[1]]
    for (cov in covariates) {
      row[[cov]] <- stats::weighted.mean(df[[cov]][idx], w, na.rm = TRUE)
    }
    row
  })
  agg <- do.call(rbind, rows)
  rownames(agg) <- NULL
  agg[order(agg$location, agg$ID_year, agg[[time_col]]), , drop = FALSE]
}

# Fans a group-level lag map back out to data.frame(location, covariate, lag),
# which is the shape add_lagged_columns requires — it errors on any missing
# (location, covariate) pair.
expand_group_lags <- function(group_lag_map, group_map) {
  merged <- merge(group_map, group_lag_map,
                  by.x = "group", by.y = "location")
  merged <- merged[order(merged$location, merged$covariate), , drop = FALSE]
  rownames(merged) <- NULL
  merged[, c("location", "covariate", "lag")]
}

### Nonlinearity backends ######################################################
#
# A backend is a function with signature
#   function(df, covariate, location_specific_effects) -> list(df, terms)
# It receives a data frame that already has the shifted `<cov>_lag` column,
# mutates it with whatever columns its parameterisation needs, and returns the
# formula fragments it wants the linear predictor to contain. The orchestrator
# `generate_lagged_model` paste()s the fragments together.
#
# Adding a backend: define the function, register it in `nonlinearity_backends`.
# Existing callers pick a backend via `user_options$nonlinearity` (config).

# Default: a single RW1 smooth on inla.group'd shifted column per covariate.
# Matches ewars_Plus's `selected_Model_form_rw`. With location_specific_effects
# additionally adds a per-location RW1 deviation on the same grouped column
# (distinct name so INLA accepts both terms) — global shape + per-location
# partial-pooled deviation.
backend_rw1_inla_group <- function(df, covariate,
                                   location_specific_effects = FALSE) {
  col <- lagged_col_name(covariate)
  grp <- paste0(col, "_grp")
  df[[grp]] <- inla.group(df[[col]])
  terms <- sprintf("f(%s, model='rw1', scale.model=TRUE)", grp)
  if (location_specific_effects) {
    grp_loc <- paste0(col, "_grp_loc")
    df[[grp_loc]] <- df[[grp]]  # INLA needs a distinct column name per f()
    terms <- c(
      terms,
      sprintf("f(%s, model='rw1', scale.model=TRUE, replicate=ID_spat)",
              grp_loc)
    )
  }
  list(df = df, terms = terms)
}

# Simpler baseline: linear in the *standardised* shifted column. Standardising
# is necessary so INLA's Newton-Raphson optimiser converges on raw covariate
# scales (rainfall in mm can run into the hundreds). NA values in the shifted
# column are imputed to 0 (the standardised mean) so prediction rows with NA
# covariates still produce a finite linear predictor. With
# `location_specific_effects = TRUE`, adds a per-location random slope on the
# standardised column via `f(ID_spat_<cov>, <col>_z, model='iid')`.
backend_linear <- function(df, covariate,
                           location_specific_effects = FALSE) {
  col <- lagged_col_name(covariate)
  z_col <- paste0(col, "_z")
  mu_x <- mean(df[[col]], na.rm = TRUE)
  sd_x <- stats::sd(df[[col]], na.rm = TRUE)
  if (!is.finite(sd_x) || sd_x == 0) sd_x <- 1
  z <- (df[[col]] - mu_x) / sd_x
  z[is.na(z)] <- 0
  df[[z_col]] <- z
  terms <- z_col
  if (location_specific_effects) {
    spat_col <- paste0("ID_spat_", covariate)
    df[[spat_col]] <- as.integer(as.factor(df$location))
    terms <- c(terms,
      sprintf("f(%s, %s, model='iid')", spat_col, z_col))
  }
  list(df = df, terms = terms)
}

nonlinearity_backends <- list(
  rw1_inla_group = backend_rw1_inla_group,
  linear         = backend_linear
)

get_nonlinearity_backend <- function(name) {
  if (!name %in% names(nonlinearity_backends)) {
    stop("Unknown nonlinearity backend: '", name, "'. Available: ",
         paste(names(nonlinearity_backends), collapse = ", "), ".",
         call. = FALSE)
  }
  nonlinearity_backends[[name]]
}

### Model formula builder ######################################################

generate_lagged_model <- function(df, covariates, lag_map, region_seasonal,
                                  location_specific_effects = FALSE,
                                  nonlinearity = backend_rw1_inla_group) {
  df <- add_lagged_columns(df, covariates, lag_map)

  smooth_terms <- character()
  for (cov in covariates) {
    res <- nonlinearity(df, cov, location_specific_effects)
    df          <- res$df
    smooth_terms <- c(smooth_terms, res$terms)
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
  list(formula = stats::as.formula(formula_str), data = df)
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
  list(formula = stats::as.formula(formula_str), data = df)
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

# Validates `lag_grouping` and downgrades "parent" to "location" when the frame
# carries no usable parent ids — grouping on the "-" placeholder would put every
# location in one group, which is silently wrong rather than loudly wrong.
resolve_lag_grouping <- function(historic_df, lag_grouping) {
  lag_grouping <- lag_grouping %||% "location"
  if (!lag_grouping %in% c("location", "parent")) {
    stop("Unknown lag_grouping: '", lag_grouping,
         "'. Available: location, parent.", call. = FALSE)
  }
  if (lag_grouping == "parent" && !has_usable_parents(historic_df)) {
    warning("lag_grouping='parent' requested but the data has no usable ",
            "`parent` column (chap-core fills '", parent_placeholder,
            "' when the dataset has no geojson). Falling back to per-location ",
            "selection.", call. = FALSE)
    lag_grouping <- "location"
  }
  lag_grouping
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
  lag_grouping <- resolve_lag_grouping(historic_df, user_options$lag_grouping)

  selection_df <- historic_df
  group_map <- NULL
  if (lag_grouping == "parent") {
    group_map <- location_group_map(historic_df)
    selection_df <- aggregate_to_parent(historic_df, covariates)
    message("Grouping lag selection by parent org unit: ",
            length(unique(group_map$group)), " groups for ",
            nrow(group_map), " locations.")
  }

  message("Selecting lags from candidates [",
          paste(candidate_lags, collapse = ", "),
          "] with ", n_folds, "-fold expanding-window CV (per-",
          lag_grouping, ")...")
  scores <- select_lags_per_district(
    selection_df, covariates, candidate_lags, n_folds = n_folds
  )
  result <- pick_best_lag_per_location_covariate(scores)
  if (lag_grouping == "parent") {
    result <- expand_group_lags(result, group_map)
  }
  message("Selected lags: ", format_lag_map(result))
  result
}
