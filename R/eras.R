#' Detect eras with different conductivity decay characteristics
#'
#' Splits the analysis results into eras based on step changes in the
#' deep-water baseline or decay amplitude A. Designed to find sensor swaps,
#' fouling, or calibration changes — typically produces 1-6 eras per site.
#'
#' The detection uses a monthly-smoothed A time series and finds step changes
#' that exceed `change_threshold` relative to the running median.
#'
#' @param analysis Data frame from [cond_analyse()]
#' @param min_era_days Minimum era duration in days (default 60)
#' @param change_threshold Relative change in monthly A that triggers a new
#'   era (default 0.4 = 40% change)
#' @param max_eras Maximum number of eras (default 6). If more are detected,
#'   only the largest step changes are kept.
#' @return The input data frame with an additional `era` column (integer)
#' @export
cond_detect_eras <- function(analysis, min_era_days = 60,
                             change_threshold = 0.4, max_eras = 6) {
  df <- analysis[analysis$fit_ok, ]
  if (nrow(df) < 10) {
    analysis$era <- 1L
    return(analysis)
  }

  df <- df[order(df$datetime), ]
  dt_num <- as.numeric(df$datetime)

  # Monthly median of A (smooths out per-profile noise)
  month_breaks <- seq(min(dt_num), max(dt_num), by = 30 * 86400)
  month_idx <- findInterval(dt_num, month_breaks)
  monthly_A <- tapply(df$A, month_idx, stats::median, na.rm = TRUE)
  monthly_dt <- tapply(dt_num, month_idx, stats::median)
  monthly_A <- as.numeric(monthly_A)
  monthly_dt <- as.numeric(monthly_dt)

  if (length(monthly_A) < 3) {
    analysis$era <- 1L
    return(analysis)
  }

  # Find step changes in monthly A
  diffs <- abs(diff(monthly_A))
  rel_diffs <- diffs / pmax(monthly_A[-length(monthly_A)], 1)

  # Score each potential boundary
  candidates <- which(rel_diffs > change_threshold)
  if (length(candidates) == 0) {
    analysis$era <- 1L
    return(analysis)
  }

  # Enforce minimum era duration and max eras
  # Sort by magnitude, greedily keep the biggest that are far enough apart
  candidate_scores <- rel_diffs[candidates]
  ord <- order(candidate_scores, decreasing = TRUE)
  candidates <- candidates[ord]
  candidate_scores <- candidate_scores[ord]

  kept <- integer(0)
  for (c in candidates) {
    if (length(kept) >= (max_eras - 1)) break
    c_dt <- monthly_dt[c]
    # Check distance from existing boundaries
    if (length(kept) == 0 ||
        all(abs(monthly_dt[kept] - c_dt) > min_era_days * 86400)) {
      kept <- c(kept, c)
    }
  }

  if (length(kept) == 0) {
    analysis$era <- 1L
    return(analysis)
  }

  # Convert monthly boundaries to datetime boundaries
  boundary_dts <- sort(monthly_dt[kept])

  # Assign eras to profiles
  era <- findInterval(dt_num, boundary_dts) + 1L
  df$era <- era

  # Map back to full analysis
  analysis$era <- NA_integer_
  analysis$era[match(df$profile_id, analysis$profile_id)] <- df$era

  # Fill NA eras from nearest neighbour
  for (i in which(is.na(analysis$era))) {
    dt_i <- as.numeric(analysis$datetime[i])
    valid <- which(!is.na(analysis$era))
    if (length(valid) > 0) {
      nearest <- valid[which.min(abs(as.numeric(analysis$datetime[valid]) - dt_i))]
      analysis$era[i] <- analysis$era[nearest]
    }
  }

  analysis
}
