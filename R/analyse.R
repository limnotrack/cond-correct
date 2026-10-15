#' Analyse conductivity equilibration artifacts across profiles
#'
#' Fits an exponential decay model (`excess = A * exp(-depth / tau)`) to the
#' **mixed surface layer only** of each profile. Within the mixed layer,
#' conductivity should be uniform with depth — any depth-dependent signal is
#' purely sensor artifact (cell equilibration). This avoids contamination
#' from real stratification-driven conductivity structure.
#'
#' The mixed layer is detected from the temperature profile using a tight
#' threshold (default 0.05°C) to ensure genuine uniformity.
#'
#' @param pro Data frame with columns for datetime, depth, conductivity, and
#'   temperature.
#' @param col_datetime Name of the datetime column (default "DateTime")
#' @param col_depth Name of the depth column (default "DptSns")
#' @param col_cond Name of the conductivity column (default "Cond")
#' @param col_temp Name of the temperature column (default "TmpWtr")
#' @param profile_id Optional integer vector of profile IDs. If NULL, profiles
#'   are detected automatically from depth resets and time gaps.
#' @param min_points Minimum points within the mixed layer to attempt fitting
#'   (default 5)
#' @param delta_t Temperature threshold for mixed layer detection (°C,
#'   default 0.05)
#' @param min_mld Minimum mixed layer depth in metres (default 3)
#' @return A data frame with one row per profile:
#'   - `profile_id`: profile identifier
#'   - `datetime`: representative datetime for the profile
#'   - `A`: fitted amplitude
#'   - `tau`: fitted decay constant (metres)
#'   - `baseline`: conductivity baseline (deepest mixed-layer reading)
#'   - `residual_sd`: fit residual standard deviation
#'   - `n_points`: number of mixed-layer points used in fit
#'   - `mld`: mixed layer depth (m)
#'   - `fit_ok`: logical, whether the fit converged with sensible params
#' @export
cond_analyse <- function(pro,
                         col_datetime = "DateTime",
                         col_depth = "DptSns",
                         col_cond = "Cond",
                         col_temp = "TmpWtr",
                         profile_id = NULL,
                         min_points = 5,
                         delta_t = 0.05,
                         min_mld = 3) {

  datetime <- pro[[col_datetime]]
  depth    <- pro[[col_depth]]
  cond     <- pro[[col_cond]]
  temp     <- pro[[col_temp]]

  if (is.null(profile_id)) {
    profile_id <- detect_profiles(datetime, depth)
  }

  pids <- sort(unique(profile_id))
  results <- vector("list", length(pids))

  for (i in seq_along(pids)) {
    pid <- pids[i]
    idx <- which(profile_id == pid)

    d <- depth[idx]; c_ <- cond[idx]; t_ <- temp[idx]; dt <- datetime[idx]

    valid <- !is.na(d) & !is.na(c_) & !is.na(t_)
    if (sum(valid) < min_points) {
      results[[i]] <- .empty_result(pid, dt)
      next
    }

    d_v <- d[valid]; c_v <- c_[valid]; t_v <- t_[valid]; dt_v <- dt[valid]
    ord <- order(d_v)
    d_v <- d_v[ord]; c_v <- c_v[ord]; t_v <- t_v[ord]

    # Find mixed layer
    mld <- find_mixed_layer_depth(d_v, t_v, delta_t = delta_t, min_depth = min_mld)

    if (is.na(mld)) {
      # No usable mixed layer — can't fit reliably
      results[[i]] <- data.frame(
        profile_id = pid, datetime = dt_v[1],
        A = NA_real_, tau = NA_real_, baseline = NA_real_,
        residual_sd = NA_real_, n_points = 0L,
        mld = NA_real_, fit_ok = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }

    # Subset to mixed layer only
    ml_mask <- d_v <= mld
    d_ml <- d_v[ml_mask]; c_ml <- c_v[ml_mask]

    if (length(d_ml) < min_points) {
      results[[i]] <- data.frame(
        profile_id = pid, datetime = dt_v[1],
        A = NA_real_, tau = NA_real_, baseline = NA_real_,
        residual_sd = NA_real_, n_points = length(d_ml),
        mld = mld, fit_ok = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }

    # Baseline: deepest readings in the mixed layer
    n_base <- min(3, length(d_ml))
    baseline <- stats::median(tail(c_ml, n_base))
    excess <- c_ml - baseline

    # Only fit if there's a clear surface excess (> 1 µS/cm)
    if (max(excess[1:min(3, length(excess))]) < 1) {
      results[[i]] <- data.frame(
        profile_id = pid, datetime = dt_v[1],
        A = 0, tau = NA_real_, baseline = baseline,
        residual_sd = 0, n_points = length(d_ml),
        mld = mld, fit_ok = TRUE,  # valid result: no artifact detected
        stringsAsFactors = FALSE
      )
      next
    }

    fit_result <- .fit_exponential(d_ml, excess)

    results[[i]] <- data.frame(
      profile_id = pid, datetime = dt_v[1],
      A = fit_result$A, tau = fit_result$tau,
      baseline = baseline, residual_sd = fit_result$residual_sd,
      n_points = length(d_ml), mld = mld,
      fit_ok = fit_result$ok,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}


#' Fit exponential decay: excess = A * exp(-depth / tau)
#' @keywords internal
.fit_exponential <- function(depth, excess) {
  # Log-linear initial estimate (fast)
  pos <- excess > 0.1
  if (sum(pos) < 3) {
    return(list(A = 0, tau = NA_real_, residual_sd = NA_real_, ok = FALSE))
  }

  log_fit <- tryCatch({
    stats::lm(log(excess[pos]) ~ depth[pos])
  }, error = function(e) NULL)

  if (is.null(log_fit)) {
    A_init <- max(excess, na.rm = TRUE)
    tau_init <- 10
  } else {
    A_init <- exp(stats::coef(log_fit)[1])
    tau_init <- -1 / stats::coef(log_fit)[2]
    if (!is.finite(tau_init) || tau_init < 1 || tau_init > 50) tau_init <- 10
    if (!is.finite(A_init) || A_init < 0) A_init <- max(excess, na.rm = TRUE)
  }

  # Refine with nls
  df <- data.frame(depth = depth, excess = excess)
  nls_fit <- tryCatch({
    fit <- stats::nls(excess ~ A * exp(-depth / tau), data = df,
               start = list(A = A_init, tau = tau_init),
               control = stats::nls.control(maxiter = 30, warnOnly = TRUE))
    cf <- stats::coef(fit)
    ok <- cf["A"] > 0 & cf["tau"] > 0.5 & cf["tau"] < 50
    res_sd <- tryCatch(summary(fit)$sigma, error = function(e) NA_real_)
    list(A = unname(cf["A"]), tau = unname(cf["tau"]),
         residual_sd = res_sd, ok = ok)
  }, error = function(e) NULL)

  if (!is.null(nls_fit)) return(nls_fit)

  # Fall back to log-linear estimate
  list(A = A_init, tau = tau_init, residual_sd = NA_real_, ok = FALSE)
}


#' @keywords internal
.empty_result <- function(pid, dt) {
  data.frame(
    profile_id = pid, datetime = dt[1],
    A = NA_real_, tau = NA_real_, baseline = NA_real_,
    residual_sd = NA_real_, n_points = 0L,
    mld = NA_real_, fit_ok = FALSE,
    stringsAsFactors = FALSE
  )
}
