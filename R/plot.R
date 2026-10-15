#' Plot example profiles showing raw vs corrected conductivity
#'
#' @param pro Data frame with profile data
#' @param analysis Data frame from [cond_analyse()]
#' @param params Data frame from [cond_optimise()] or list with A, tau
#' @param col_datetime,col_depth,col_cond,col_temp Column names
#' @param profile_id Integer vector of profile IDs
#' @param n_profiles Number of example profiles to plot (default 12)
#' @param era Optional integer vector for era-specific correction
#' @return A ggplot object
#' @export
cond_plot_profiles <- function(pro, analysis, params,
                               col_datetime = "DateTime",
                               col_depth = "DptSns",
                               col_cond = "Cond",
                               col_temp = "TmpWtr",
                               profile_id = NULL,
                               n_profiles = 12,
                               era = NULL) {

  if (is.null(profile_id)) {
    profile_id <- detect_profiles(pro[[col_datetime]], pro[[col_depth]])
  }

  # Pick profiles spread across the record
  good <- analysis[analysis$fit_ok, ]
  if (nrow(good) > n_profiles) {
    idx <- round(seq(1, nrow(good), length.out = n_profiles))
    good <- good[idx, ]
  }

  plot_data <- lapply(seq_len(nrow(good)), function(i) {
    pid <- good$profile_id[i]
    mask <- profile_id == pid
    d <- pro[[col_depth]][mask]
    c_ <- pro[[col_cond]][mask]
    valid <- !is.na(d) & !is.na(c_)
    if (sum(valid) < 3) return(NULL)

    # Handle list (from cond_optimise) or data frame
    p_df <- if (is.list(params) && "params" %in% names(params)) params$params else params

    if (is.data.frame(p_df) && "era" %in% names(p_df) && !is.null(era)) {
      e <- era[mask][valid][1]
      p_row <- p_df[p_df$era == e, ]
      if (nrow(p_row) == 0) return(NULL)
      if ("A1" %in% names(p_row)) {
        corrected <- c_[valid] - p_row$A1[1] * exp(-d[valid] / p_row$tau1[1]) -
                                  p_row$A2[1] * exp(-d[valid] / p_row$tau2[1])
      } else {
        corrected <- c_[valid] - p_row$A[1] * exp(-d[valid] / p_row$tau[1])
      }
    } else if (is.data.frame(p_df)) {
      if ("A1" %in% names(p_df)) {
        corrected <- c_[valid] - p_df$A1[1] * exp(-d[valid] / p_df$tau1[1]) -
                                  p_df$A2[1] * exp(-d[valid] / p_df$tau2[1])
      } else {
        corrected <- c_[valid] - p_df$A[1] * exp(-d[valid] / p_df$tau[1])
      }
    } else {
      corrected <- c_[valid] - p_df$A * exp(-d[valid] / p_df$tau)
    }

    data.frame(
      profile = format(good$datetime[i], "%Y-%m-%d %H:%M"),
      depth = rep(d[valid], 2),
      value = c(c_[valid], corrected),
      type = rep(c("Raw", "Corrected"), each = sum(valid))
    )
  })

  plot_data <- do.call(rbind, plot_data)
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    message("No plottable profiles")
    return(ggplot2::ggplot())
  }

  ggplot2::ggplot(plot_data, ggplot2::aes(x = value, y = depth, colour = type)) +
    ggplot2::geom_path(linewidth = 0.5) +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::scale_y_reverse() +
    ggplot2::scale_colour_manual(values = c(Raw = "grey60", Corrected = "#E41A1C")) +
    ggplot2::facet_wrap(~ profile, scales = "free_x", ncol = 4) +
    ggplot2::labs(x = "Conductivity", y = "Depth (m)", colour = NULL) +
    ggplot2::theme_bw(base_size = 9)
}


#' Plot diagnostic summary of conductivity analysis
#'
#' Four-panel plot: (1) A over time with era boundaries, (2) tau over time,
#' (3) baseline over time, (4) histogram of residuals.
#'
#' @param analysis Data frame from [cond_analyse()] with `era` column
#' @param params Optional data frame from [cond_optimise()]
#' @return A patchwork plot object
#' @export
cond_plot_diagnostics <- function(analysis, params = NULL) {
  df <- analysis[analysis$fit_ok, ]

  a_col <- if ("A" %in% names(df)) "A" else "A1"
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = datetime, y = .data[[a_col]])) +
    ggplot2::geom_point(size = 0.5, alpha = 0.3) +
    ggplot2::labs(x = NULL, y = "A (amplitude)", title = "Decay amplitude") +
    ggplot2::theme_bw(base_size = 9)

  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = datetime, y = tau)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.3) +
    ggplot2::labs(x = NULL, y = "tau (m)", title = "Decay constant") +
    ggplot2::theme_bw(base_size = 9)

  # Clamp A axis to remove outliers
  A_ylim <- c(0, stats::quantile(df$A, 0.99, na.rm = TRUE) * 1.2)
  p1 <- p1 + ggplot2::coord_cartesian(ylim = A_ylim)

  p3 <- ggplot2::ggplot(df, ggplot2::aes(x = datetime, y = baseline)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.3) +
    ggplot2::labs(x = NULL, y = "Baseline (µS/cm)", title = "Mixed-layer baseline") +
    ggplot2::theme_bw(base_size = 9)

  # Mixed layer depth panel
  if ("mld" %in% names(df)) {
    p_mld <- ggplot2::ggplot(df[!is.na(df$mld), ],
                              ggplot2::aes(x = datetime, y = mld)) +
      ggplot2::geom_point(size = 0.5, alpha = 0.3) +
      ggplot2::scale_y_reverse() +
      ggplot2::labs(x = NULL, y = "MLD (m)", title = "Mixed layer depth used for fitting") +
      ggplot2::theme_bw(base_size = 9)
  }

  p4 <- ggplot2::ggplot(df[!is.na(df$residual_sd), ],
                         ggplot2::aes(x = residual_sd)) +
    ggplot2::geom_histogram(bins = 40, fill = "steelblue", colour = "white") +
    ggplot2::labs(x = "Residual SD", y = "Count", title = "Fit quality") +
    ggplot2::theme_bw(base_size = 9)

  # Add era boundaries and optimised values
  if ("era" %in% names(df)) {
    era_bounds <- df |>
      dplyr::group_by(era) |>
      dplyr::summarise(dt_start = min(datetime), .groups = "drop")
    p1 <- p1 + ggplot2::geom_vline(data = era_bounds,
                                    ggplot2::aes(xintercept = dt_start),
                                    linetype = "dashed", colour = "red", linewidth = 0.3)
    p2 <- p2 + ggplot2::geom_vline(data = era_bounds,
                                    ggplot2::aes(xintercept = dt_start),
                                    linetype = "dashed", colour = "red", linewidth = 0.3)

    if (!is.null(params)) {
      # Handle both single (A) and double (A1) exp params
      if (is.list(params) && "params" %in% names(params)) params <- params$params
      era_params <- dplyr::left_join(era_bounds, params, by = "era")
      a_param <- if ("A" %in% names(era_params)) "A" else "A1"
      t_param <- if ("tau" %in% names(era_params)) "tau" else "tau1"
      p1 <- p1 + ggplot2::geom_hline(data = era_params,
                                       ggplot2::aes(yintercept = .data[[a_param]]),
                                       colour = "red", linewidth = 0.5)
      p2 <- p2 + ggplot2::geom_hline(data = era_params,
                                       ggplot2::aes(yintercept = .data[[t_param]]),
                                       colour = "red", linewidth = 0.5)
    }
  }

  if ("mld" %in% names(df)) {
    (p1 | p2) / (p3 | p_mld) / p4 +
      patchwork::plot_annotation(title = "Conductivity equilibration analysis (mixed-layer fits)")
  } else {
    (p1 | p2) / (p3 | p4) +
      patchwork::plot_annotation(title = "Conductivity equilibration analysis")
  }
}


#' Plot era timeline
#'
#' Horizontal bar chart showing eras with their A and tau values.
#'
#' @param analysis Data frame with `era` column
#' @param params Data frame from [cond_optimise()]
#' @return A ggplot object
#' @export
cond_plot_eras <- function(analysis, params) {
  df <- analysis[analysis$fit_ok & !is.na(analysis$era), ]

  era_summary <- df |>
    dplyr::group_by(era) |>
    dplyr::summarise(
      dt_start = min(datetime),
      dt_end = max(datetime),
      n = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::left_join(params, by = "era")

  ggplot2::ggplot(era_summary) +
    ggplot2::geom_segment(ggplot2::aes(x = dt_start, xend = dt_end,
                                        y = era, yend = era,
                                        colour = A),
                           linewidth = 4) +
    ggplot2::scale_colour_viridis_c(name = "A") +
    ggplot2::geom_text(ggplot2::aes(x = dt_start, y = era,
                                     label = sprintf("A=%.1f tau=%.1f", A, tau)),
                        hjust = 0, vjust = -1, size = 3) +
    ggplot2::labs(x = NULL, y = "Era", title = "Conductivity correction eras") +
    ggplot2::theme_bw(base_size = 10)
}
