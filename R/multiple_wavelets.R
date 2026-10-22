# ---------------------------------------------------------------------------
# Multiple wavelets: MODWT-based wavelet correlation (WC), cross-correlation
# (WCC), wavelet multiple correlation (WMC), wavelet multiple cross-correlation
# (WMCC) and wavelet local multiple correlation (WLMC).
# Computations follow Fernandez-Macho (2012) / Polanco-Martinez (W2CWM2C);
# plotting is re-implemented with ggplot2.
# ---------------------------------------------------------------------------

FILTER_LENGTH <- c(haar = 2, d4 = 4, d6 = 6, d8 = 8, d16 = 16, la8 = 8, la16 = 16,
                   la20 = 20, bl14 = 14, bl20 = 20, c6 = 6, c12 = 12)

FILTER_CHOICES <- c("Least asymmetric LA(8)" = "la8", "LA(16)" = "la16", "LA(20)" = "la20",
                    "Haar" = "haar", "Daubechies D(4)" = "d4", "D(6)" = "d6", "D(8)" = "d8",
                    "D(16)" = "d16", "Best localised BL(14)" = "bl14", "BL(20)" = "bl20",
                    "Coiflet C(6)" = "c6", "C(12)" = "c12")

WINDOW_CHOICES <- c("Gaussian" = "gauss", "Uniform" = "unif", "Bartlett (triangular)" = "bart",
                    "Epanechnikov" = "epan", "Tricube (Cleveland)" = "clev", "Wendland" = "wend")

# Largest J such that brick-walled coefficients remain and CI formula is defined
max_levels <- function(N, wf = "la8", min_coefs = 10) {
  L  <- FILTER_LENGTH[[wf]]
  js <- 1:20
  ok <- (N - (2^js - 1) * (L - 1) >= min_coefs) & (N / 2^js - 3 > 0)
  if (!any(ok)) return(1L)
  as.integer(max(which(ok)))
}

level_labels <- function(J, unit = "obs") {
  sprintf("d%d  [%d-%d %s]", seq_len(J), 2^seq_len(J), 2^(seq_len(J) + 1), unit)
}

modwt_list <- function(vals, wf, J) {
  L <- lapply(vals, function(x) waveslim::brick.wall(waveslim::modwt(as.numeric(x), wf, n.levels = J), wf))
  names(L) <- names(vals)
  L
}

# Plain MODWT (not brick-walled) of one series as a long data.frame for plotting
modwt_decomposition <- function(x, dates, wf, J) {
  w <- waveslim::modwt(as.numeric(x), wf, n.levels = J)
  idx <- if (is.null(dates)) seq_along(x) else dates
  comps <- c("Series", paste0("d", seq_len(J)), paste0("s", J))
  vals  <- c(list(as.numeric(x)), w[seq_len(J)], w[J + 1])
  long  <- do.call(rbind, lapply(seq_along(comps), function(i)
    data.frame(time = idx, component = comps[i], value = vals[[i]])))
  long$component <- factor(long$component, levels = comps)
  long
}

plot_modwt <- function(long, title = NULL) {
  ggplot2::ggplot(long, ggplot2::aes(time, value)) +
    ggplot2::geom_line(colour = WA_COLORS$series, linewidth = 0.35) +
    ggplot2::facet_wrap(~component, ncol = 1, scales = "free_y", strip.position = "right") +
    wa_labs(x = NULL, y = NULL, title = title) +
    app_theme()
}

# ---- Pairwise wavelet correlation (WC) --------------------------------------
pairwise_wc <- function(L, N, J) {
  nm  <- names(L)
  cmb <- utils::combn(length(L), 2)
  out <- do.call(rbind, lapply(seq_len(ncol(cmb)), function(k) {
    i <- cmb[1, k]; j <- cmb[2, k]
    w <- as.data.frame(as.matrix(waveslim::wave.correlation(L[[i]], L[[j]], N = N)))[seq_len(J), ]
    names(w) <- c("wavecor", "lower", "upper")
    data.frame(pair = paste(nm[i], "~", nm[j]), var1 = nm[i], var2 = nm[j],
               level = seq_len(J), scale = 2^(seq_len(J) - 1),
               wavecor = w$wavecor, lower = w$lower, upper = w$upper,
               stringsAsFactors = FALSE)
  }))
  out$significant <- !is.na(out$lower) & (out$lower > 0 | out$upper < 0)
  out
}

plot_pairwise_wc <- function(df, J, unit = "obs", title = "Pairwise wavelet correlation") {
  sums <- tapply(df$wavecor, df$pair, sum, na.rm = TRUE)
  df$pair  <- factor(df$pair, levels = names(sort(sums, decreasing = TRUE)))
  df$level_lab <- factor(level_labels(J, unit)[df$level], levels = level_labels(J, unit))
  df$label <- ifelse(is.na(df$wavecor), "", sprintf("%.2f%s", df$wavecor, ifelse(df$significant, "*", "")))
  npairs <- nlevels(df$pair)
  if (npairs > 30) df$label <- ""          # too crowded: colours only (use the matrix layout instead)
  ggplot2::ggplot(df, ggplot2::aes(pair, level_lab, fill = wavecor)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = label, colour = abs(wavecor) > 0.55),
                       size = if (npairs > 15) 2.6 else 3.4, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = WA_COLORS$ink)) +
    ggplot2::scale_fill_gradient2(low = WA_COLORS$div_low, mid = WA_COLORS$div_mid, high = WA_COLORS$div_high,
                                  limits = c(-1, 1), name = "Correlation", na.value = "#e1e0d9") +
    wa_labs(x = NULL, y = "Wavelet level  [period band]", title = title,
                  subtitle = "* 95% confidence interval excludes zero; pairs sorted by overall correlation") +
    app_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = if (npairs > 6) 45 else 0,
                                                       hjust = if (npairs > 6) 1 else 0.5),
                   panel.grid = ggplot2::element_blank())
}

# Alternative layout: one correlation matrix (lower triangle) per wavelet level.
# Readable for any number of variables.
plot_wc_matrix <- function(df, vars, J, unit = "obs", title = "Pairwise wavelet correlation by level") {
  df$level_lab <- factor(level_labels(J, unit)[df$level], levels = level_labels(J, unit))
  df$row <- factor(df$var2, levels = rev(vars))
  df$col <- factor(df$var1, levels = vars)
  df$label <- ifelse(is.na(df$wavecor), "", sprintf("%.2f%s", df$wavecor, ifelse(df$significant, "*", "")))
  nv <- length(vars)
  txt <- if (nv <= 8) 3 else if (nv <= 12) 2.3 else 0
  p <- ggplot2::ggplot(df, ggplot2::aes(col, row, fill = wavecor)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_gradient2(low = WA_COLORS$div_low, mid = WA_COLORS$div_mid, high = WA_COLORS$div_high,
                                  limits = c(-1, 1), name = "Correlation", na.value = "#e1e0d9") +
    ggplot2::facet_wrap(~level_lab, ncol = if (J <= 4) J else 3) +
    wa_labs(x = NULL, y = NULL, title = title,
                  subtitle = "* 95% confidence interval excludes zero") +
    app_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = if (nv > 10) 8 else 10),
                   axis.text.y = ggplot2::element_text(size = if (nv > 10) 8 else 10),
                   panel.grid = ggplot2::element_blank())
  if (txt > 0)
    p <- p + ggplot2::geom_text(ggplot2::aes(label = label, colour = abs(wavecor) > 0.55), size = txt, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = WA_COLORS$ink))
  p
}

# ---- Pairwise wavelet cross-correlation (WCC) -------------------------------
pair_wcc <- function(Lx, Ly, lmax, N, J, name_x = "x", name_y = "y") {
  lags <- -lmax:lmax
  m <- sapply(seq_len(J), function(j) waveslim::spin.correlation(Lx[[j]], Ly[[j]], lmax))
  m <- matrix(m, nrow = length(lags), ncol = J)
  n_j <- matrix(trunc(N / 2^(seq_len(J))), nrow = length(lags), ncol = J, byrow = TRUE)
  se  <- 1 / sqrt(pmax(n_j - 3, 1))
  lower <- tanh(atanh(m) - stats::qnorm(0.975) * se)
  upper <- tanh(atanh(m) + stats::qnorm(0.975) * se)
  df <- data.frame(lag = rep(lags, J), level = rep(seq_len(J), each = length(lags)),
                   scale = rep(2^(seq_len(J) - 1), each = length(lags)),
                   wcc = as.vector(m), lower = as.vector(lower), upper = as.vector(upper))
  df$significant <- !is.na(df$lower) & (df$lower > 0 | df$upper < 0)
  attr(df, "pair") <- c(name_x, name_y)
  df
}

# shared heatmap for WCC / WMCC (lag x level)
plot_lag_heatmap <- function(df, value, J, unit = "obs", title = "", side_labels = NULL,
                             only_sig = TRUE, lag_unit = "obs") {
  df$val <- df[[value]]
  df$fill <- if (only_sig) ifelse(df$significant, df$val, NA) else df$val
  df$level_lab <- factor(level_labels(J, unit)[df$level], levels = level_labels(J, unit))
  maxlag <- do.call(rbind, lapply(split(df, df$level), function(d) {
    d <- d[!is.na(d$val), ]
    if (!nrow(d)) return(NULL)
    d[which.max(d$val), c("lag", "level_lab")]
  }))
  lmax <- max(abs(df$lag))
  p <- ggplot2::ggplot(df, ggplot2::aes(lag, level_lab)) +
    ggplot2::geom_tile(ggplot2::aes(fill = fill)) +
    ggplot2::scale_fill_gradient2(low = WA_COLORS$div_low, mid = WA_COLORS$div_mid, high = WA_COLORS$div_high,
                                  limits = c(-1, 1), name = "Correlation", na.value = "#e1e0d9") +
    ggplot2::geom_vline(xintercept = 0, colour = WA_COLORS$axis, linewidth = 0.5) +
    wa_labs(x = sprintf("Lag (%s)", lag_unit), y = "Wavelet level  [period band]",
                  title = title,
                  subtitle = if (only_sig) "Only values whose 95% CI excludes zero are coloured; black tick = lag of maximum correlation" else "Black tick = lag of maximum correlation") +
    app_theme() + ggplot2::theme(panel.grid = ggplot2::element_blank())
  if (!is.null(maxlag) && nrow(maxlag))
    p <- p + ggplot2::geom_segment(data = maxlag,
                                   ggplot2::aes(x = lag, xend = lag,
                                                y = as.numeric(level_lab) - 0.5, yend = as.numeric(level_lab) + 0.5),
                                   linewidth = 1.1, colour = WA_COLORS$ink)
  if (!is.null(side_labels))
    p <- p + ggplot2::annotate("label", x = lmax, y = seq_len(J), label = side_labels,
                               hjust = 1, size = 3.2, fill = "white", alpha = 0.85)
  p
}

plot_lag_lines <- function(df, value, J, unit = "obs", title = "", lag_unit = "obs") {
  df$val <- df[[value]]
  df$level_lab <- factor(level_labels(J, unit)[df$level], levels = level_labels(J, unit))
  ggplot2::ggplot(df, ggplot2::aes(lag, val)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = WA_COLORS$series, alpha = 0.18) +
    ggplot2::geom_hline(yintercept = 0, colour = WA_COLORS$axis, linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0, colour = WA_COLORS$axis, linewidth = 0.5) +
    ggplot2::geom_line(colour = WA_COLORS$series, linewidth = 0.7) +
    ggplot2::facet_wrap(~level_lab, ncol = 3) +
    wa_labs(x = sprintf("Lag (%s)", lag_unit), y = "Correlation", title = title,
                  subtitle = "Shaded band: 95% confidence interval") +
    app_theme()
}

# ---- Wavelet multiple correlation (WMC) -------------------------------------
wmc_calc <- function(L, N, J) {
  r <- wavemulcor::wave.multiple.correlation(L, N = N)
  m <- r$xy.mulcor[seq_len(J), , drop = FALSE]
  data.frame(level = seq_len(J), scale = 2^(seq_len(J) - 1),
             period_from = 2^seq_len(J), period_to = 2^(seq_len(J) + 1),
             wmc = m[, 1], lower = m[, 2], upper = m[, 3],
             max_var = names(L)[r$YmaxR[seq_len(J)]], stringsAsFactors = FALSE)
}

plot_wmc <- function(df, unit = "obs", title = "Wavelet multiple correlation") {
  df$level_lab <- factor(level_labels(nrow(df), unit), levels = level_labels(nrow(df), unit))
  ggplot2::ggplot(df, ggplot2::aes(level_lab, wmc, group = 1)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = WA_COLORS$series, alpha = 0.18) +
    ggplot2::geom_line(colour = WA_COLORS$series, linewidth = 0.9) +
    ggplot2::geom_point(colour = WA_COLORS$series, size = 2.6) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.15, colour = WA_COLORS$series) +
    ggplot2::geom_text(ggplot2::aes(y = pmin(lower, 0.98) - 0.06, label = max_var), size = 3.3, colour = "grey20") +
    ggplot2::geom_hline(yintercept = 0, colour = WA_COLORS$axis, linewidth = 0.5) +
    ggplot2::scale_y_continuous(limits = c(-0.1, 1.05), breaks = seq(0, 1, 0.2)) +
    wa_labs(x = "Wavelet level  [period band]", y = "Wavelet multiple correlation", title = title,
                  subtitle = "Band: 95% CI. Label below each point: variable best explained by the others at that scale") +
    app_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

# ---- Wavelet multiple cross-correlation (WMCC) ------------------------------
wmcc_calc <- function(L, lmax, J) {
  r <- wavemulcor::wave.multiple.cross.correlation(L, lag.max = lmax)
  lags <- -lmax:lmax
  val <- r$xy.mulcor[seq_len(J), , drop = FALSE]
  lo  <- r$ci.mulcor$lower[seq_len(J), , drop = FALSE]
  up  <- r$ci.mulcor$upper[seq_len(J), , drop = FALSE]
  df <- data.frame(lag = rep(lags, each = J), level = rep(seq_len(J), length(lags)),
                   scale = rep(2^(seq_len(J) - 1), length(lags)),
                   wmcc = as.vector(val), lower = as.vector(lo), upper = as.vector(up))
  df$significant <- !is.na(df$lower) & (df$lower > 0 | df$upper < 0)
  attr(df, "max_var") <- names(L)[r$YmaxR[seq_len(J)]]
  df
}

# ---- Wavelet local multiple correlation (WLMC) ------------------------------
wlmc_calc <- function(L, M, window, J, dates) {
  r <- wavemulcor::wave.local.multiple.correlation(L, M = M, window = window)
  N <- nrow(r$val)
  idx <- if (is.null(dates)) seq_len(N) else dates
  cols <- seq_len(J)
  df <- data.frame(time = rep(idx, J), level = rep(seq_len(J), each = N),
                   wlmc = as.vector(as.matrix(r$val[, cols])),
                   lower = as.vector(as.matrix(r$lo[, cols])),
                   upper = as.vector(as.matrix(r$up[, cols])))
  df
}

plot_wlmc <- function(df, J, unit = "obs", title = "Wavelet local multiple correlation") {
  df$level_lab <- factor(level_labels(J, unit)[df$level], levels = level_labels(J, unit))
  # tile width = typical spacing of observations, so weekly/monthly data draw without gaps
  tt <- sort(unique(df$time))
  tile_w <- if (length(tt) > 1) stats::median(diff(as.numeric(tt))) else 1
  ggplot2::ggplot(df, ggplot2::aes(time, level_lab, fill = wlmc)) +
    ggplot2::geom_tile(width = tile_w) +
    ggplot2::scale_fill_gradientn(colours = WA_COLORS$seq, limits = c(0, 1), name = "WLMC", na.value = "#e1e0d9") +
    wa_labs(x = NULL, y = "Wavelet level  [period band]", title = title,
                  subtitle = "Time-varying multiple correlation (moving window); each row is one wavelet scale") +
    app_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

plot_wlmc_lines <- function(df, J, unit = "obs", title = "Wavelet local multiple correlation by scale") {
  df$level_lab <- factor(level_labels(J, unit)[df$level], levels = level_labels(J, unit))
  ggplot2::ggplot(df, ggplot2::aes(time, wlmc)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = WA_COLORS$series, alpha = 0.18) +
    ggplot2::geom_line(colour = WA_COLORS$series, linewidth = 0.6) +
    ggplot2::facet_wrap(~level_lab, ncol = 2) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    wa_labs(x = NULL, y = "WLMC", title = title, subtitle = "Shaded band: 95% confidence interval") +
    app_theme()
}
