#' Identify profiles from DateTime and DptSns
#'
#' Groups consecutive readings into profiles based on depth resets
#' (large negative depth change) or time gaps.
#'
#' @param datetime POSIXct vector
#' @param depth Numeric depth vector
#' @param gap_hours Time gap (hours) that triggers a new profile
#' @param reset_m Depth decrease (metres) that triggers a new profile
#' @return Integer vector of profile IDs
#' @keywords internal
detect_profiles <- function(datetime, depth, gap_hours = 2, reset_m = 5) {
  n <- length(datetime)
  if (n == 0L) return(integer(0))
  dt_num <- as.numeric(datetime)
  gap_secs <- gap_hours * 3600
  pid <- integer(n)
  current <- 1L
  pid[1] <- current
  for (i in 2:n) {
    time_gap <- is.na(dt_num[i]) || is.na(dt_num[i - 1L]) ||
      (dt_num[i] - dt_num[i - 1L]) > gap_secs
    depth_reset <- !is.na(depth[i]) && !is.na(depth[i - 1L]) &&
      (depth[i] - depth[i - 1L]) < -reset_m
    if (time_gap || depth_reset) current <- current + 1L
    pid[i] <- current
  }
  pid
}


#' Find the mixed layer depth from a temperature profile
#'
#' Descends from the shallowest reading until temperature deviates from
#' the surface value by more than `delta_t`. Uses a tight threshold to
#' ensure the mixed layer is genuinely uniform — important for
#' conductivity correction where small temperature gradients drive real
#' conductivity changes.
#'
#' @param depth Numeric depth vector (need not be sorted)
#' @param temp Numeric temperature vector
#' @param delta_t Temperature deviation threshold (°C). Default 0.05 —
#'   intentionally tight to avoid including any stratified water.
#' @param min_depth Minimum mixed layer depth to return (m). If the mixed
#'   layer is shallower than this, returns NA (too thin to fit).
#' @return Mixed layer depth (numeric scalar), or NA if not found or too thin
#' @keywords internal
find_mixed_layer_depth <- function(depth, temp, delta_t = 0.05, min_depth = 3) {
  valid <- !is.na(depth) & !is.na(temp)
  if (sum(valid) < 4) return(NA_real_)
  d <- depth[valid]
  t <- temp[valid]
  ord <- order(d)
  d <- d[ord]; t <- t[ord]

  # Surface reference: median of shallowest 3 readings
  n_surf <- min(3, length(t))
  t_surf <- stats::median(t[1:n_surf])

  # Walk down until temperature deviates
  for (i in (n_surf + 1):length(t)) {
    if (abs(t[i] - t_surf) > delta_t) {
      mld <- d[i - 1]
      if (mld >= min_depth) return(mld)
      return(NA_real_)
    }
  }

  # Entire profile is mixed
  max(d)
}


#' Get the mixed layer mask for a profile
#'
#' Returns a logical vector: TRUE for readings within the well-mixed
#' surface layer, where conductivity should be depth-uniform and any
#' depth-dependent signal is purely sensor artifact.
#'
#' @param depth Numeric depth vector
#' @param temp Numeric temperature vector
#' @param delta_t Temperature threshold for mixed layer detection (°C)
#' @param min_depth Minimum mixed layer depth (m)
#' @return Logical vector: TRUE = within mixed layer
#' @keywords internal
mixed_layer_mask <- function(depth, temp, delta_t = 0.05, min_depth = 3) {
  mld <- find_mixed_layer_depth(depth, temp, delta_t = delta_t,
                                 min_depth = min_depth)
  if (is.na(mld)) return(rep(FALSE, length(depth)))
  !is.na(depth) & depth <= mld
}
