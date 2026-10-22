# ---------------------------------------------------------------------------
# Bivariate wavelets: continuous Morlet wavelet power, cross-wavelet power,
# wavelet coherence and phase differences via WaveletComp.
# ---------------------------------------------------------------------------

SIM_METHODS <- c("White noise" = "white.noise", "Shuffle" = "shuffle",
                 "Fourier randomisation" = "Fourier.rand", "AR(1)" = "AR", "ARIMA" = "ARIMA")

DJ_CHOICES <- c("1/20 (default)" = "1/20", "1/50" = "1/50", "1/100" = "1/100", "1/250 (slow)" = "1/250")

WINDOW_TYPES <- c("Bartlett" = 1, "Triangular (non-Bartlett)" = 2, "Box car (uniform)" = 3,
                  "Epanechnikov" = 4, "Gauss" = 5, "Hanning" = 6)

run_coherency <- function(dates, vals, pair, dj = 1/20, lowerPeriod = 2, upperPeriod = NULL,
                          loess.span = 0, make.pval = TRUE, method = "white.noise", n.sim = 100,
                          window.type.t = 1, window.type.s = 1, window.size.t = 5, window.size.s = 1/4) {
  stopifnot(length(pair) == 2, all(pair %in% names(vals)))
  if (pair[1] == pair[2]) stop("Please select two different variables.")
  N <- nrow(vals)
  if (is.null(upperPeriod) || is.na(upperPeriod) || upperPeriod <= lowerPeriod) upperPeriod <- floor(N / 3)
  upperPeriod <- min(upperPeriod, N)
  my.data <- data.frame(x = vals[[pair[1]]], y = vals[[pair[2]]])
  names(my.data) <- pair
  show_date <- !is.null(dates)
  if (show_date) my.data$date <- as.POSIXct(dates, tz = "UTC")
  wc <- WaveletComp::analyze.coherency(my.data, my.pair = pair, loess.span = loess.span, dt = 1, dj = dj,
                                       lowerPeriod = lowerPeriod, upperPeriod = upperPeriod,
                                       window.type.t = window.type.t, window.type.s = window.type.s,
                                       window.size.t = window.size.t, window.size.s = window.size.s,
                                       make.pval = make.pval, method = method, n.sim = n.sim,
                                       verbose = FALSE)
  list(wc = wc, pair = pair, show_date = show_date, make_pval = make.pval)
}

.leg <- function(lab) list(lab = lab, mar = 5.1, label.digits = 2, n.ticks = 6)

plot_wt_power <- function(b, series, unit, siglvl = 0.05, palette = "blue") {
  op <- title_par(); on.exit(graphics::par(op), add = TRUE)
  WaveletComp::wt.image(b$wc, my.series = series, color.key = "interval", siglvl = siglvl,
                        color.palette = wc_palette_string(palette), col.contour = WA_COLORS$ink,
                        legend.params = .leg("wavelet power levels"),
                        main = fit_main(paste("Wavelet power spectrum:", series)),
                        show.date = b$show_date, date.format = "%F",
                        timelab = "Time", periodlab = paste0("Period (", unit_plural(unit), ")"),
                        graphics.reset = TRUE)
}

plot_cross_power <- function(b, unit, siglvl = 0.05, palette = "blue") {
  op <- title_par(); on.exit(graphics::par(op), add = TRUE)
  WaveletComp::wc.image(b$wc, which.image = "wp", color.key = "interval", siglvl.arrow = siglvl,
                        siglvl.contour = siglvl, legend.params = .leg("cross-wavelet power levels"),
                        color.palette = wc_palette_string(palette), col.contour = WA_COLORS$ink,
                        main = fit_main(paste("Cross-wavelet power:", b$pair[1], "and", b$pair[2])),
                        show.date = b$show_date, date.format = "%F",
                        timelab = "Time", periodlab = paste0("Period (", unit_plural(unit), ")"),
                        graphics.reset = TRUE)
}

plot_coherence <- function(b, unit, siglvl = 0.05, palette = "blue") {
  op <- title_par(); on.exit(graphics::par(op), add = TRUE)
  WaveletComp::wc.image(b$wc, which.image = "wc", color.key = "interval", siglvl.arrow = siglvl,
                        siglvl.contour = siglvl, legend.params = .leg("wavelet coherence levels"),
                        color.palette = wc_palette_string(palette), col.contour = WA_COLORS$ink,
                        main = fit_main(paste("Wavelet coherence:", b$pair[1], "and", b$pair[2])),
                        show.date = b$show_date, date.format = "%F",
                        timelab = "Time", periodlab = paste0("Period (", unit_plural(unit), ")"),
                        graphics.reset = TRUE)
}

plot_avg <- function(b, which = "wc", unit) {
  lab <- if (which == "wc") "average coherence" else "average cross-wavelet power"
  op <- title_par(); on.exit(graphics::par(op), add = TRUE)
  WaveletComp::wc.avg(b$wc, which.avg = which, siglvl = c(0.05, 0.1), sigpch = 20,
                      sigcol = c(WA_COLORS$div_high, WA_COLORS$series), col = WA_COLORS$ink2, lwd = 1.5,
                      main = fit_main(paste0(toupper(substr(lab, 1, 1)), substr(lab, 2, nchar(lab)), ": ",
                                             b$pair[1], " and ", b$pair[2]), frac = 0.9),
                      averagelab = lab, periodlab = paste0("Period (", unit_plural(unit), ")"))
}

plot_phase <- function(b, lower, upper, siglvl = 0.05, which.sig = "wc") {
  op <- title_par(); on.exit(graphics::par(op), add = TRUE)
  WaveletComp::wc.sel.phases(b$wc, sel.lower = lower, sel.upper = upper, only.sig = b$make_pval,
                             which.sig = which.sig, siglvl = siglvl, legend.coords = "topright",
                             phase.cols = c(WA_COLORS$div_high, WA_COLORS$series), Angle.col = WA_COLORS$ink,
                             legend.horiz = FALSE, phaselim = c(-pi, pi),
                             main = fit_main(sprintf("Phase difference %s vs %s, periods %s-%s", b$pair[1], b$pair[2],
                                                     format(round(lower, 1)), format(round(upper, 1))), frac = 0.9),
                             sub = "", show.date = b$show_date, date.format = "%F", timelab = "Time")
}

plot_phasediff_image <- function(b, unit, siglvl = 0.05, palette = "blue") {
  op <- title_par(); on.exit(graphics::par(op), add = TRUE)
  WaveletComp::wc.phasediff.image(b$wc, which.contour = "wc", siglvl = siglvl,
                                  color.palette = wc_palette_string(palette, diverging = TRUE), col.contour = WA_COLORS$ink,
                                  legend.params = list(lab = "phase difference levels", mar = 5.1),
                                  main = fit_main(paste("Phase difference:", b$pair[1], "and", b$pair[2])),
                                  show.date = b$show_date, date.format = "%F",
                                  timelab = "Time", periodlab = paste0("Period (", unit_plural(unit), ")"),
                                  graphics.reset = TRUE)
}

plot_reconstruct <- function(b, series, siglvl = 0.05, waves = FALSE) {
  op <- title_par(); on.exit(graphics::par(op), add = TRUE)
  WaveletComp::reconstruct(b$wc, my.series = series, only.sig = b$make_pval, siglvl = siglvl,
                           plot.waves = waves, plot.rec = !waves, lwd = c(1, 2), col = c(WA_COLORS$muted, WA_COLORS$series),
                           legend.coords = "topleft",
                           main = fit_main(paste(if (waves) "Significant waves:" else "Reconstruction:", series), frac = 0.9),
                           show.date = b$show_date, date.format = "%F", timelab = "Time", verbose = FALSE)
}

coherence_table <- function(b) {
  wc <- b$wc
  df <- data.frame(Period = wc$Period,
                   `Avg cross power` = wc$Power.xy.avg,
                   `Avg coherence` = wc$Coherence.avg, check.names = FALSE)
  if (!is.null(wc$Power.xy.avg.pval)) df$`Cross power p` <- wc$Power.xy.avg.pval
  if (!is.null(wc$Coherence.avg.pval)) df$`Coherence p` <- wc$Coherence.avg.pval
  df
}

# Mean phase difference (in radians) inside the significant coherence region, by period
phase_summary <- function(b, siglvl = 0.05) {
  wc <- b$wc
  ang <- wc$Angle
  mask <- if (!is.null(wc$Coherence.pval)) wc$Coherence.pval < siglvl else matrix(TRUE, nrow(ang), ncol(ang))
  mean_ang <- sapply(seq_len(nrow(ang)), function(i) {
    a <- ang[i, mask[i, ]]
    if (!length(a)) return(NA_real_)
    atan2(mean(sin(a)), mean(cos(a)))
  })
  share <- rowMeans(mask)
  data.frame(Period = wc$Period, `Share of time significant` = round(share, 3),
             `Mean phase diff (rad)` = round(mean_ang, 3),
             Regime = ifelse(is.na(mean_ang), "",
                             ifelse(abs(mean_ang) < pi / 2,
                                    ifelse(mean_ang > 0, "in phase, x leads", "in phase, y leads"),
                                    ifelse(mean_ang > 0, "anti-phase, y leads", "anti-phase, x leads"))),
             check.names = FALSE)
}
