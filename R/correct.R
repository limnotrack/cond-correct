#' Apply conductivity equilibration correction
#'
#' Fast vectorised correction for use in QC pipelines. Supports single
#' or double exponential models.
#'
#' @param cond Numeric vector of conductivity values
#' @param depth Numeric vector of depth values
#' @param A Amplitude (single exp). Ignored if `params` is provided.
#' @param tau Decay constant in metres (single exp). Ignored if `params`.
#' @param A1,A2,tau1,tau2 Double exponential parameters. Ignored if `params`.
#' @param era Optional integer vector for era-specific correction.
#' @param params Data frame from [cond_optimise()]$params. If provided,
#'   overrides individual A/tau arguments. Requires `era`.
#' @param clamp_min Minimum corrected value (default NULL).
#' @return Numeric vector of corrected conductivity values
#' @export
#'
#' @examples
#' # Simple single exponential
#' corrected <- cond_correct(df$Cond, df$DptSns, A = 18, tau = 13.5)
#'
#' # Double exponential
#' corrected <- cond_correct(df$Cond, df$DptSns,
#'                           A1 = 10, tau1 = 2.5, A2 = 8, tau2 = 12)
#'
#' # From optimise results (era-specific)
#' corrected <- cond_correct(df$Cond, df$DptSns,
#'                           params = opt$params, era = df$era)
cond_correct <- function(cond, depth,
                         A = NULL, tau = NULL,
                         A1 = NULL, A2 = NULL, tau1 = NULL, tau2 = NULL,
                         era = NULL, params = NULL, clamp_min = NULL) {

  n <- length(cond)

  if (!is.null(params)) {
    # Era-specific correction from params data frame
    if (is.null(era)) stop("'era' must be provided when using 'params'")
    correction <- numeric(n)
    is_double <- "A1" %in% names(params)

    for (i in seq_len(nrow(params))) {
      mask <- era == params$era[i] & !is.na(era)
      if (!any(mask)) next
      d <- depth[mask]
      if (is_double) {
        correction[mask] <- params$A1[i] * exp(-d / params$tau1[i]) +
                            params$A2[i] * exp(-d / params$tau2[i])
      } else {
        correction[mask] <- params$A[i] * exp(-d / params$tau[i])
      }
    }

  } else if (!is.null(A1) && !is.null(tau1) && !is.null(A2) && !is.null(tau2)) {
    # Double exponential (global)
    correction <- A1 * exp(-depth / tau1) + A2 * exp(-depth / tau2)

  } else if (!is.null(A) && !is.null(tau)) {
    # Single exponential (global)
    correction <- A * exp(-depth / tau)

  } else {
    stop("Provide A+tau, A1+tau1+A2+tau2, or params+era")
  }

  result <- cond - correction
  if (!is.null(clamp_min)) result <- pmax(result, clamp_min)
  result
}
