#' Optimise conductivity correction parameters
#'
#' Two-stage approach:
#' 1. **tau** (decay constant) is fitted globally from well-mixed profiles,
#'    or can be overridden manually. Supports single or double exponential.
#' 2. **A** (amplitude) is optimised per era.
#'
#' The double exponential model captures two processes:
#'   - Fast equilibration (tau1 ~ 2-3m): cell flushing
#'   - Slow thermal mass (tau2 ~ 10-20m): temperature compensation lag
#'   correction = A1 * exp(-d/tau1) + A2 * exp(-d/tau2)
#'   A1 and A2 are optimised per era; tau1 and tau2 are global.
#'
#' @param pro Data frame with profile data
#' @param analysis Data frame from [cond_analyse()] with `era` column
#' @param col_datetime,col_depth,col_cond,col_temp Column names
#' @param profile_id Integer vector of profile IDs matching `pro` rows
#' @param model "single" or "double" exponential (default "double")
#' @param tau Override for tau (single exp) or c(tau1, tau2) (double exp).
#'   If NULL, estimated from deep-mixed profiles.
#' @param A_range Range for amplitude search (default c(0, 40))
#' @param A_n Grid points for amplitude search (default 20)
#' @param mld_threshold Min mixed layer depth for tau estimation (default 15)
#' @param delta_t Temperature threshold for mixed layer (default 0.2)
#' @return A list with:
#'   - `model`: "single" or "double"
#'   - `tau`: global tau (single) or c(tau1, tau2) (double)
#'   - `tau_n`: number of profiles used for tau estimation
#'   - `params`: data frame with era, A (or A1/A2), tau (or tau1/tau2)
#' @export
cond_optimise <- function(pro,
                          analysis,
                          col_datetime = "DateTime",
                          col_depth = "DptSns",
                          col_cond = "Cond",
                          col_temp = "TmpWtr",
                          profile_id = NULL,
                          model = "double",
                          tau = NULL,
                          A_range = c(0, 40),
                          A_n = 20,
                          mld_threshold = 15,
                          delta_t = 0.2) {

  depth <- pro[[col_depth]]
  cond  <- pro[[col_cond]]
  temp  <- if (!is.null(col_temp) && col_temp %in% names(pro)) pro[[col_temp]] else NULL

  if (is.null(profile_id)) {
    profile_id <- detect_profiles(pro[[col_datetime]], depth)
  }

  if (!"era" %in% names(analysis)) {
    analysis$era <- 1L
  }

  # ===================================================================
  # Stage 1: Estimate tau globally
  # ===================================================================
  tau_n <- 0L

  if (!is.null(tau)) {
    # Manual override
    if (model == "double") {
      stopifnot(length(tau) == 2)
      tau1 <- tau[1]; tau2 <- tau[2]
    } else {
      tau1 <- tau[1]; tau2 <- NULL
    }
    message(sprintf("Using manual tau: %s", paste(round(tau, 1), collapse = ", ")))
  } else {
    # Estimate from deep-mixed profiles
    deep_mixed <- analysis[analysis$fit_ok &
                            !is.na(analysis$mld) &
                            analysis$mld >= mld_threshold, ]
    tau_n <- nrow(deep_mixed)

    if (model == "double") {
      # Fit double exponential to deep-mixed profiles
      tau_fit <- .estimate_double_tau(pro, deep_mixed, profile_id,
                                      col_depth, col_cond, col_temp, delta_t)
      tau1 <- tau_fit$tau1; tau2 <- tau_fit$tau2
      message(sprintf("Double exp: tau1=%.1f m (fast), tau2=%.1f m (slow) from %d profiles",
                      tau1, tau2, tau_n))
    } else {
      valid_tau <- deep_mixed$tau[!is.na(deep_mixed$tau) &
                                   deep_mixed$tau > 1 & deep_mixed$tau < 40]
      tau1 <- if (length(valid_tau) >= 5) stats::median(valid_tau) else 5
      tau2 <- NULL
      message(sprintf("Single exp: tau=%.1f m from %d profiles", tau1, tau_n))
    }
  }

  # ===================================================================
  # Stage 2: Optimise amplitudes per era
  # ===================================================================
  eras <- sort(unique(stats::na.omit(analysis$era)))
  results <- vector("list", length(eras))

  for (i in seq_along(eras)) {
    e <- eras[i]
    era_pids <- analysis$profile_id[analysis$era == e & analysis$fit_ok]

    if (length(era_pids) < 3) {
      era_rows <- analysis[analysis$era == e & analysis$fit_ok, ]
      if (model == "double") {
        results[[i]] <- data.frame(era = e, A1 = NA_real_, A2 = NA_real_,
                                    tau1 = tau1, tau2 = tau2,
                                    n_profiles = length(era_pids), objective = NA_real_)
      } else {
        results[[i]] <- data.frame(era = e, A = NA_real_, tau = tau1,
                                    n_profiles = length(era_pids), objective = NA_real_)
      }
      next
    }

    # Subsample
    use_pids <- if (length(era_pids) > 50) {
      era_pids[round(seq(1, length(era_pids), length.out = 50))]
    } else era_pids

    # Pre-extract data for speed
    era_data <- lapply(use_pids, function(pid) {
      idx <- which(profile_id == pid)
      d <- depth[idx]; c_ <- cond[idx]
      t_ <- if (!is.null(temp)) temp[idx] else NULL
      valid <- !is.na(d) & !is.na(c_) & if (!is.null(t_)) !is.na(t_) else TRUE
      if (sum(valid) < 5) return(NULL)
      d_v <- d[valid]; c_v <- c_[valid]
      t_v <- if (!is.null(t_)) t_[valid] else NULL
      ord <- order(d_v); d_v <- d_v[ord]; c_v <- c_v[ord]
      if (!is.null(t_v)) t_v <- t_v[ord]
      mask <- if (!is.null(t_v)) {
        mixed_layer_mask(d_v, t_v, delta_t = delta_t, min_depth = 3)
      } else rep(TRUE, length(d_v))
      if (sum(mask) < 4) return(NULL)
      list(d = d_v[mask], c = c_v[mask])
    })
    era_data <- Filter(Negate(is.null), era_data)

    if (length(era_data) < 3) {
      if (model == "double") {
        results[[i]] <- data.frame(era = e, A1 = 0, A2 = 0,
                                    tau1 = tau1, tau2 = tau2,
                                    n_profiles = length(era_pids), objective = NA_real_)
      } else {
        results[[i]] <- data.frame(era = e, A = 0, tau = tau1,
                                    n_profiles = length(era_pids), objective = NA_real_)
      }
      next
    }

    if (model == "double") {
      opt_result <- .optimise_double_A(era_data, tau1, tau2, A_range, A_n)
      results[[i]] <- data.frame(era = e, A1 = opt_result$A1, A2 = opt_result$A2,
                                  tau1 = tau1, tau2 = tau2,
                                  n_profiles = length(era_pids),
                                  objective = opt_result$objective)
    } else {
      opt_result <- .optimise_single_A(era_data, tau1, A_range, A_n)
      results[[i]] <- data.frame(era = e, A = opt_result$A, tau = tau1,
                                  n_profiles = length(era_pids),
                                  objective = opt_result$objective)
    }
  }

  result_df <- do.call(rbind, results)
  tau_out <- if (model == "double") c(tau1, tau2) else tau1
  list(model = model, tau = tau_out, tau_n = tau_n, params = result_df)
}


#' Estimate tau1 and tau2 for double exponential from deep-mixed profiles
#' @keywords internal
.estimate_double_tau <- function(pro, deep_mixed, profile_id,
                                 col_depth, col_cond, col_temp, delta_t) {
  depth <- pro[[col_depth]]
  cond  <- pro[[col_cond]]
  temp  <- pro[[col_temp]]

  # Sample up to 100 deep-mixed profiles
  pids <- deep_mixed$profile_id
  if (length(pids) > 100) pids <- pids[round(seq(1, length(pids), length.out = 100))]

  # Collect excess curves
  all_excess <- list()
  for (pid in pids) {
    idx <- which(profile_id == pid)
    d <- depth[idx]; c_ <- cond[idx]; t_ <- temp[idx]
    valid <- !is.na(d) & !is.na(c_) & !is.na(t_)
    if (sum(valid) < 8) next
    d_v <- d[valid]; c_v <- c_[valid]; t_v <- t_[valid]
    ord <- order(d_v); d_v <- d_v[ord]; c_v <- c_v[ord]; t_v <- t_v[ord]
    mask <- mixed_layer_mask(d_v, t_v, delta_t = delta_t, min_depth = 10)
    if (sum(mask) < 8) next
    d_ml <- d_v[mask]; c_ml <- c_v[mask]
    baseline <- stats::median(tail(c_ml, 3))
    excess <- c_ml - baseline
    if (max(excess[1:min(3, length(excess))]) < 1) next
    all_excess[[length(all_excess) + 1]] <- data.frame(depth = d_ml, excess = excess)
  }

  if (length(all_excess) < 10) {
    return(list(tau1 = 2.5, tau2 = 12))
  }

  # Pool all excess curves and fit double exponential
  pooled <- do.call(rbind, all_excess)
  pooled <- pooled[pooled$excess > 0, ]

  fit <- tryCatch({
    stats::nls(excess ~ A1 * exp(-depth / tau1) + A2 * exp(-depth / tau2),
               data = pooled,
               start = list(A1 = max(pooled$excess) * 0.6, tau1 = 2,
                            A2 = max(pooled$excess) * 0.4, tau2 = 12),
               control = stats::nls.control(maxiter = 100, warnOnly = TRUE))
  }, error = function(e) NULL)

  if (!is.null(fit)) {
    cf <- stats::coef(fit)
    # Ensure tau1 < tau2 (tau1 = fast, tau2 = slow)
    if (cf["tau1"] > cf["tau2"]) {
      return(list(tau1 = max(1, cf["tau2"]), tau2 = min(40, cf["tau1"])))
    }
    return(list(tau1 = max(1, cf["tau1"]), tau2 = min(40, cf["tau2"])))
  }

  list(tau1 = 2.5, tau2 = 12)
}


#' Optimise A for single exponential (1D)
#' @keywords internal
.optimise_single_A <- function(era_data, tau, A_range, A_n) {
  obj_fn <- function(A) {
    if (A < 0) return(1e6)
    total <- 0; n <- 0
    for (pd in era_data) {
      corrected <- pd$c - A * exp(-pd$d / tau)
      slope <- abs(stats::coef(stats::lm(corrected ~ pd$d))[2])
      total <- total + slope
      # Overcorrection penalty
      sh <- corrected[1:min(3, length(corrected))]
      dp <- tail(corrected, min(3, length(corrected)))
      if (mean(sh) < mean(dp) - 0.5) total <- total + (mean(dp) - mean(sh))^2
      n <- n + 1
    }
    if (n == 0) return(1e6)
    total / n
  }

  A_grid <- seq(A_range[1], A_range[2], length.out = A_n)
  grid_obj <- vapply(A_grid, obj_fn, numeric(1))
  best_A <- A_grid[which.min(grid_obj)]
  opt <- tryCatch(
    stats::optimize(obj_fn, c(max(0, best_A - 5), best_A + 5)),
    error = function(e) list(minimum = best_A, objective = min(grid_obj))
  )
  list(A = max(0, opt$minimum), objective = opt$objective)
}


#' Optimise A1 and A2 for double exponential (2D)
#' @keywords internal
.optimise_double_A <- function(era_data, tau1, tau2, A_range, A_n) {
  obj_fn <- function(params) {
    A1 <- params[1]; A2 <- params[2]
    if (A1 < 0 || A2 < 0) return(1e6)
    total <- 0; n <- 0
    for (pd in era_data) {
      corrected <- pd$c - A1 * exp(-pd$d / tau1) - A2 * exp(-pd$d / tau2)
      slope <- abs(stats::coef(stats::lm(corrected ~ pd$d))[2])
      total <- total + slope
      sh <- corrected[1:min(3, length(corrected))]
      dp <- tail(corrected, min(3, length(corrected)))
      if (mean(sh) < mean(dp) - 0.5) total <- total + (mean(dp) - mean(sh))^2
      n <- n + 1
    }
    if (n == 0) return(1e6)
    total / n
  }

  # 2D grid search (A_n/2 x A_n/2 for speed)
  n2 <- max(5, A_n %/% 2)
  A1_grid <- seq(A_range[1], A_range[2], length.out = n2)
  A2_grid <- seq(A_range[1], A_range[2] / 2, length.out = n2)
  grid <- expand.grid(A1 = A1_grid, A2 = A2_grid)
  grid_obj <- vapply(seq_len(nrow(grid)), function(j) {
    obj_fn(c(grid$A1[j], grid$A2[j]))
  }, numeric(1))
  best <- grid[which.min(grid_obj), ]

  opt <- tryCatch(
    stats::optim(c(best$A1, best$A2), obj_fn, method = "Nelder-Mead",
                 control = list(maxit = 200)),
    error = function(e) list(par = c(best$A1, best$A2), value = min(grid_obj))
  )

  list(A1 = max(0, opt$par[1]), A2 = max(0, opt$par[2]), objective = opt$value)
}
