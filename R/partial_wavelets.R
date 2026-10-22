# ---------------------------------------------------------------------------
# Partial wavelet coherence (biwavelet): coherence between y and x1 after
# controlling for a third series x2, following Mihanovic et al. (2009) and
# Ng & Chan (2012); biwavelet::pwtc() does the computation.  The workflow
# mirrors the author's Bi_Partial_Wavelets scripts (wtc of the pair, then
# pwtc(y, x1 | x2) and pwtc(y, x2 | x1)).
# ---------------------------------------------------------------------------

PW_DJ_CHOICES <- c("1/12 (biwavelet default)" = "1/12", "1/8 (faster)" = "1/8", "1/20 (finer)" = "1/20")

run_partial <- function(dates, vals, y, x1, x2, dj = 1/12, nrands = 50, sig.level = 0.95,
                        max.scale = NULL, compute_wtc = TRUE) {
  vars <- c(y, x1, x2)
  if (length(unique(vars)) < 3) stop("y, x1 and x2 must be three different variables.")
  if (!all(vars %in% names(vals))) stop("Selected variables not found in the data.")
  N   <- nrow(vals)
  idx <- seq_len(N)
  Y  <- cbind(idx, as.numeric(vals[[y]]))
  X1 <- cbind(idx, as.numeric(vals[[x1]]))
  X2 <- cbind(idx, as.numeric(vals[[x2]]))
  if (!is.null(max.scale) && (is.na(max.scale) || max.scale <= 2)) max.scale <- NULL
  nrands <- max(0L, as.integer(nrands))
  args <- list(dj = dj, nrands = nrands, sig.level = sig.level, quiet = TRUE)
  if (!is.null(max.scale)) args$max.scale <- max.scale
  out <- list(y = y, x1 = x1, x2 = x2, N = N, dates = dates, nrands = nrands, sig.level = sig.level,
              show_date = !is.null(dates))
  # the AR(1) fits inside biwavelet's Monte-Carlo test often emit optim warnings; they are harmless
  if (compute_wtc) {
    out$wtc_y_x1 <- suppressWarnings(do.call(biwavelet::wtc, c(list(Y, X1), args)))
    out$wtc_y_x2 <- suppressWarnings(do.call(biwavelet::wtc, c(list(Y, X2), args)))
  }
  out$pwtc_y_x1 <- suppressWarnings(do.call(biwavelet::pwtc, c(list(Y, X1, X2), args)))
  out$pwtc_y_x2 <- suppressWarnings(do.call(biwavelet::pwtc, c(list(Y, X2, X1), args)))
  out
}

# date axis for biwavelet plots drawn on an observation index
.pw_time_axis <- function(dates, N) {
  if (is.null(dates)) { graphics::axis(1); return(invisible()) }
  yrs <- format(dates, "%Y")
  starts <- match(unique(yrs), yrs)
  if (length(starts) > 12) starts <- starts[seq(1, length(starts), by = ceiling(length(starts) / 12))]
  if (length(starts) < 3) {        # short samples: label months
    ym <- format(dates, "%Y-%m"); starts <- match(unique(ym), ym)
    if (length(starts) > 12) starts <- starts[seq(1, length(starts), by = ceiling(length(starts) / 12))]
    graphics::axis(1, at = starts, labels = format(dates[starts], "%b %Y"), las = 2, cex.axis = 0.8)
  } else {
    graphics::axis(1, at = starts, labels = yrs[starts])
  }
  invisible()
}

plot_biwavelet <- function(obj, dates, unit = "observation", title = "", plot.phase = TRUE, plot.sig = TRUE, palette = "blue") {
  op <- graphics::par(oma = c(0, 0, 0, 1), mar = c(5, 4.5, 4, 5) + 0.1)
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(obj, plot.phase = plot.phase, plot.sig = plot.sig, lty.coi = 1, col.coi = "grey", lwd.coi = 2,
                 lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.12, plot.cb = TRUE, fill.cols = seq_ramp(64, palette),
                 col.sig = WA_COLORS$ink,
                 ylab = paste0("Period (", unit_plural(unit), ")"), xlab = "Time", main = title, xaxt = "n")
  .pw_time_axis(dates, length(obj$t))
  invisible()
}

plot_partial_pair <- function(p, which = 1, unit = "observation", palette = "blue") {
  if (which == 1) plot_biwavelet(p$pwtc_y_x1, p$dates, unit,
                                 sprintf("Partial wavelet coherence: %s ~ %s | %s", p$y, p$x1, p$x2), palette = palette)
  else            plot_biwavelet(p$pwtc_y_x2, p$dates, unit,
                                 sprintf("Partial wavelet coherence: %s ~ %s | %s", p$y, p$x2, p$x1), palette = palette)
}

plot_wtc_pair <- function(p, which = 1, unit = "observation", palette = "blue") {
  obj <- if (which == 1) p$wtc_y_x1 else p$wtc_y_x2
  other <- if (which == 1) p$x1 else p$x2
  if (is.null(obj)) stop("Wavelet coherence was not computed.")
  plot_biwavelet(obj, p$dates, unit, sprintf("Wavelet coherence: %s and %s", p$y, other), palette = palette)
}

# WTC above, PWTC below, for the same pair
plot_partial_compare <- function(p, which = 1, unit = "observation", palette = "blue") {
  obj  <- if (which == 1) p$wtc_y_x1 else p$wtc_y_x2
  pobj <- if (which == 1) p$pwtc_y_x1 else p$pwtc_y_x2
  other <- if (which == 1) p$x1 else p$x2
  ctrl  <- if (which == 1) p$x2 else p$x1
  op <- graphics::par(mfrow = c(2, 1), oma = c(0, 0, 0, 1), mar = c(4, 4.5, 3, 5) + 0.1)
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(obj, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, lwd.sig = 2,
                 arrow.lwd = 0.03, arrow.len = 0.1, plot.cb = TRUE, xaxt = "n", fill.cols = seq_ramp(64, palette), col.sig = WA_COLORS$ink,
                 ylab = paste0("Period (", unit_plural(unit), ")"), xlab = "",
                 main = sprintf("Wavelet coherence: %s and %s", p$y, other))
  .pw_time_axis(p$dates, length(obj$t))
  graphics::plot(pobj, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, lwd.sig = 2,
                 arrow.lwd = 0.03, arrow.len = 0.1, plot.cb = TRUE, xaxt = "n", fill.cols = seq_ramp(64, palette), col.sig = WA_COLORS$ink,
                 ylab = paste0("Period (", unit_plural(unit), ")"), xlab = "Time",
                 main = sprintf("Partial wavelet coherence: %s ~ %s | %s", p$y, other, ctrl))
  .pw_time_axis(p$dates, length(pobj$t))
  invisible()
}

# Period-by-period summary: mean coherence inside the cone of influence and
# share of time significant, for WTC and PWTC of both pairs
.pw_period_stats <- function(obj) {
  if (is.null(obj)) return(NULL)
  rsq <- obj$rsq
  inside <- outer(obj$period, obj$coi, FUN = "<=")      # TRUE where period <= coi (inside cone)
  rsq_in <- ifelse(inside, rsq, NA)
  avg <- rowMeans(rsq_in, na.rm = TRUE)
  # biwavelet stores the significance threshold per cell; a cell is significant
  # where rsq / signif >= 1 (that is the contour plot.biwavelet draws)
  sig <- if (!is.null(obj$signif) && is.matrix(obj$signif) && all(dim(obj$signif) == dim(rsq))) {
    s <- ifelse(inside, (rsq / obj$signif) >= 1, NA); rowMeans(s, na.rm = TRUE)
  } else rep(NA_real_, length(obj$period))
  data.frame(Period = obj$period, avg = avg, share_sig = sig)
}

partial_table <- function(p) {
  a <- .pw_period_stats(p$wtc_y_x1); b <- .pw_period_stats(p$pwtc_y_x1)
  c2 <- .pw_period_stats(p$wtc_y_x2); d <- .pw_period_stats(p$pwtc_y_x2)
  df <- data.frame(Period = round(b$Period, 2))
  if (!is.null(a)) df[[sprintf("WTC %s~%s", p$y, p$x1)]] <- round(a$avg, 3)
  df[[sprintf("PWTC %s~%s|%s", p$y, p$x1, p$x2)]] <- round(b$avg, 3)
  if (!is.null(c2)) df[[sprintf("WTC %s~%s", p$y, p$x2)]] <- round(c2$avg, 3)
  df[[sprintf("PWTC %s~%s|%s", p$y, p$x2, p$x1)]] <- round(d$avg, 3)
  if (p$nrands > 0) {
    df[[sprintf("Sig. share PWTC %s~%s|%s", p$y, p$x1, p$x2)]] <- round(b$share_sig, 3)
    df[[sprintf("Sig. share PWTC %s~%s|%s", p$y, p$x2, p$x1)]] <- round(d$share_sig, 3)
  }
  df
}

# Long table of the partial coherence matrix (period x time) for export
partial_matrix <- function(obj, dates) {
  m <- obj$rsq
  dimnames(m) <- list(round(obj$period, 3), if (!is.null(dates)) format(dates) else seq_len(ncol(m)))
  m
}
