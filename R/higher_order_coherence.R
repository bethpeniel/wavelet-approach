# ---------------------------------------------------------------------------
# Higher-order wavelet coherencies: multiple wavelet coherence (MWC) and
# partial wavelet coherence (PWC) with one or several controls.
#
# Theory and conventions follow
#   * Aguiar-Conraria & Soares (2014, J. Econ. Surveys), Appendix A: the p x p
#     Hermitian matrix S of smoothed (cross-)wavelet spectra; squared multiple
#     coherency R^2 = 1 - det(S) / (S_11 det(S^11))  (A1), and complex partial
#     coherency gamma_1j.q = -S^d_j1 / sqrt(S^d_11 S^d_jj) (A2), i.e. partial
#     correlation applied to S; the ASToolbox (MPAWCOG.m) implements these.
#   * Hu & Si (2016, HESS): MWC = S_yX S_XX^-1 S_Xy / S_yy with Monte-Carlo
#     AR(1) significance (per-scale 95th percentile outside the cone).
#   * Hu & Si (2021, HESS): complex PWC with several excluding variables and
#     its phase; shows the "classical" one-control formula of Mihanovic et al.
#     (2009) / Ng & Chan (2012) uses real instead of complex coherence.
#   * Ng & Chan (2012, JTECH): PWC / MWC for one control / two predictors.
#   * Grinsted, Moore & Jevrejeva (2004): smoothing operator and significance
#     test design; wavelet transform and smoothing come from biwavelet.
#
# Both quantities reduce to the partial spectral matrix (Koopmans 1995):
#   S_AA.Z = S_AA - S_AZ S_ZZ^-1 S_ZA
#   MWC(y | X):      R^2 = 1 - S_yy.X / S_yy
#   PWC(y, x | Z):   gamma = S_yx.Z / sqrt(S_yy.Z S_xx.Z);  rsq = |gamma|^2,
#                    phase = Arg(gamma)
# ---------------------------------------------------------------------------

HO_DJ_CHOICES <- c("1/12 (biwavelet default)" = "1/12", "1/8 (faster)" = "1/8", "1/20 (finer)" = "1/20")
HO_PWC_METHODS <- c("Complex coherence, any number of controls (Hu & Si 2021; Aguiar-Conraria & Soares 2014)" = "complex",
                    "Classical real-coherence form, one control (Mihanovic et al. 2009; Ng & Chan 2012)" = "classical")

# ---- smoothed spectral matrix ------------------------------------------------
# Returns list(S = S[[i]][[j]] complex matrices (scales x times), wt = list of
# wt objects, coi, period, scale, dt, t, vars, ar1)
ho_spectra <- function(vals, vars, dj = 1/12, max.scale = NULL, pad = TRUE) {
  N   <- nrow(vals)
  idx <- seq_len(N)
  dt  <- 1
  s0  <- 2 * dt
  if (is.null(max.scale) || is.na(max.scale) || max.scale <= 2) max.scale <- (N * 0.17) * 2 * dt
  J1 <- round(log2(max.scale / s0) / dj)
  wts <- lapply(vars, function(v)
    biwavelet::wt(cbind(idx, as.numeric(vals[[v]])), pad = pad, dj = dj, s0 = s0, J1 = J1,
                  max.scale = max.scale, mother = "morlet", do.sig = FALSE))
  names(wts) <- vars
  ar1 <- vapply(vars, function(v) {
    a <- tryCatch(stats::arima(as.numeric(vals[[v]]), order = c(1, 0, 0), method = "CSS")$coef[[1]], error = function(e) NA_real_)
    if (is.na(a)) a <- stats::cor(vals[[v]][-1], vals[[v]][-N])
    max(min(a, 0.98), -0.98)
  }, numeric(1))
  out <- ho_smooth_matrix(lapply(wts, `[[`, "wave"), wts[[1]]$scale, dt, dj)
  out$coi    <- Reduce(pmin, lapply(wts, `[[`, "coi"))
  out$period <- wts[[1]]$period
  out$scale  <- wts[[1]]$scale
  out$dt <- dt; out$dj <- dj; out$s0 <- s0; out$J1 <- J1; out$max.scale <- max.scale; out$pad <- pad
  out$t <- idx; out$vars <- vars; out$ar1 <- ar1; out$N <- N
  out
}

# smoothed auto/cross spectra of a list of wave matrices (biwavelet convention:
# scale-normalised, then smoothed in time and scale)
ho_smooth_matrix <- function(waves, scale, dt, dj) {
  p <- length(waves)
  n <- ncol(waves[[1]])
  s.inv <- matrix(rep(1 / scale, n), nrow = length(scale))
  S <- vector("list", p)
  for (i in seq_len(p)) {
    S[[i]] <- vector("list", p)
    for (j in seq_len(p)) {
      if (j < i) { S[[i]][[j]] <- Conj(S[[j]][[i]]); next }
      cw <- if (i == j) abs(waves[[i]])^2 else waves[[i]] * Conj(waves[[j]])
      S[[i]][[j]] <- biwavelet::smooth.wavelet(s.inv * cw, dt, dj, scale)
    }
  }
  list(S = S, p = p)
}

# ---- partial spectral matrix S_AA.Z (elementwise over the time-frequency plane)
# A, Z: integer index vectors into S. Returns list-of-lists over A.
# Closed-form inverses for |Z| <= 3; generic per-cell solve otherwise.
ho_partial_spec <- function(S, A, Z) {
  q <- length(Z)
  getS <- function(i, j) S[[i]][[j]]
  if (q == 0) {
    return(lapply(A, function(i) lapply(A, function(j) getS(i, j))))
  }
  # elements of S_ZZ^-1 as a q x q list of matrices
  if (q == 1) {
    inv <- list(list(1 / Re(getS(Z[1], Z[1]))))
  } else if (q == 2) {
    a <- Re(getS(Z[1], Z[1])); c <- Re(getS(Z[2], Z[2])); b <- getS(Z[1], Z[2])
    det <- a * c - Mod(b)^2
    inv <- list(list(c / det, -b / det), list(-Conj(b) / det, a / det))
  } else if (q == 3) {
    a <- Re(getS(Z[1], Z[1])); d <- Re(getS(Z[2], Z[2])); f <- Re(getS(Z[3], Z[3]))
    b <- getS(Z[1], Z[2]); c <- getS(Z[1], Z[3]); e <- getS(Z[2], Z[3])
    # Hermitian [[a, b, c], [b*, d, e], [c*, e*, f]]; adjugate / determinant
    det <- Re(a * d * f + 2 * Re(b * e * Conj(c)) - a * Mod(e)^2 - d * Mod(c)^2 - f * Mod(b)^2)
    i11 <- (d * f - Mod(e)^2) / det
    i12 <- (c * Conj(e) - b * f) / det
    i13 <- (b * e - c * d) / det
    i22 <- (a * f - Mod(c)^2) / det
    i23 <- (Conj(b) * c - a * e) / det
    i33 <- (a * d - Mod(b)^2) / det
    inv <- list(list(i11, i12, i13), list(Conj(i12), i22, i23), list(Conj(i13), Conj(i23), i33))
  } else {
    inv <- ho_inverse_generic(S, Z)
  }
  # S_AA.Z = S_AA - S_AZ inv S_ZA
  out <- vector("list", length(A))
  for (ia in seq_along(A)) {
    out[[ia]] <- vector("list", length(A))
    for (ib in seq_along(A)) {
      if (ib < ia) { out[[ia]][[ib]] <- Conj(out[[ib]][[ia]]); next }
      acc <- 0
      for (k in seq_len(q)) for (l in seq_len(q))
        acc <- acc + getS(A[ia], Z[k]) * inv[[k]][[l]] * getS(Z[l], A[ib])
      out[[ia]][[ib]] <- getS(A[ia], A[ib]) - acc
    }
  }
  out
}

# generic elementwise inverse of the Hermitian q x q sub-matrix S[Z, Z] (q >= 4)
ho_inverse_generic <- function(S, Z) {
  q <- length(Z)
  dims <- dim(S[[Z[1]]][[Z[1]]])
  ncell <- prod(dims)
  arr <- array(NA_complex_, c(ncell, q, q))
  for (k in seq_len(q)) for (l in seq_len(q)) arr[, k, l] <- as.vector(S[[Z[k]]][[Z[l]]])
  invarr <- array(NA_complex_, c(ncell, q, q))
  for (cell in seq_len(ncell)) {
    m <- matrix(arr[cell, , ], q, q)
    invarr[cell, , ] <- tryCatch(solve(m), error = function(e) matrix(NA_complex_, q, q))
  }
  lapply(seq_len(q), function(k) lapply(seq_len(q), function(l) matrix(invarr[, k, l], dims[1], dims[2])))
}

# ---- the two statistics -----------------------------------------------------
# multiple wavelet coherence of series y_idx on predictors X_idx (squared, real)
ho_mwc_stat <- function(S, y_idx, X_idx) {
  ps <- ho_partial_spec(S, A = y_idx, Z = X_idx)
  rsq <- 1 - Re(ps[[1]][[1]]) / Re(S[[y_idx]][[y_idx]])
  rsq[!is.finite(rsq)] <- NA
  pmin(pmax(rsq, 0), 1)
}

# complex partial wavelet coherence of y_idx and x_idx given controls Z_idx
ho_pwc_stat <- function(S, y_idx, x_idx, Z_idx, method = "complex") {
  if (method == "classical") {
    # Mihanovic et al. (2009) as implemented by Ng & Chan (2012) / biwavelet::pwtc:
    # real coherences R, one control z
    z <- Z_idx[1]
    R2 <- function(i, j) Mod(S[[i]][[j]])^2 / (Re(S[[i]][[i]]) * Re(S[[j]][[j]]))
    ryx <- sqrt(R2(y_idx, x_idx)); ryz <- sqrt(R2(y_idx, z)); rxz <- sqrt(R2(x_idx, z))
    rsq <- (ryx - ryz * rxz)^2 / ((1 - ryz^2) * (1 - rxz^2))
    phase <- Arg(S[[y_idx]][[x_idx]])          # plain cross-spectrum phase (no partial phase in this form)
  } else {
    ps <- ho_partial_spec(S, A = c(y_idx, x_idx), Z = Z_idx)
    gamma <- ps[[1]][[2]] / sqrt(Re(ps[[1]][[1]]) * Re(ps[[2]][[2]]))
    rsq   <- Mod(gamma)^2
    phase <- Arg(gamma)
  }
  rsq[!is.finite(rsq)] <- NA
  list(rsq = pmin(pmax(rsq, 0), 1), phase = phase)
}

# ---- Monte-Carlo significance -----------------------------------------------
# AR(1) surrogates for every series (Grinsted et al. 2004; Hu & Si 2016, 2021):
# the statistic is recomputed nrands times and, per scale, the sig.level
# quantile of all values outside the cone of influence is the threshold.
ho_signif <- function(spec, stat_fun, nrands = 100, sig.level = 0.95, progress = NULL) {
  if (nrands < 1) return(NULL)
  N <- spec$N; p <- spec$p
  outside <- outer(spec$period, spec$coi, FUN = "<=")   # TRUE where period <= coi
  nscale <- length(spec$scale)
  vals <- vector("list", nscale)
  sim_ar1 <- function(phi, n) {
    x <- tryCatch(biwavelet:::ar1_ma0_sim(0, phi, n), error = function(e) NULL)
    if (is.null(x) || length(x) != n) x <- as.numeric(stats::arima.sim(list(ar = phi), n = n))
    x
  }
  for (r in seq_len(nrands)) {
    sim <- as.data.frame(lapply(spec$ar1, function(phi) sim_ar1(phi, N)))
    names(sim) <- spec$vars
    idx <- seq_len(N)
    wts <- lapply(spec$vars, function(v)
      biwavelet::wt(cbind(idx, sim[[v]]), pad = spec$pad, dj = spec$dj, s0 = spec$s0, J1 = spec$J1,
                    max.scale = spec$max.scale, mother = "morlet", do.sig = FALSE))
    Ssim <- ho_smooth_matrix(lapply(wts, `[[`, "wave"), spec$scale, spec$dt, spec$dj)$S
    st <- stat_fun(Ssim)
    for (s in seq_len(nscale)) {
      v <- st[s, outside[s, ]]
      vals[[s]] <- c(vals[[s]], v[is.finite(v)])
    }
    if (is.function(progress)) progress(r / nrands)
  }
  thr <- vapply(vals, function(v) if (length(v)) stats::quantile(v, sig.level, names = FALSE) else NA_real_, numeric(1))
  matrix(thr, nrow = nscale, ncol = ncol(outside))
}

# ---- driver -----------------------------------------------------------------
# y: dependent; X: predictors (for MWC) / predictor-and-control set (for PWC)
# Computes MWC(y | X) and, for every x in X, PWC(y, x | X \ x). Optionally the
# plain bivariate coherence of each pair for comparison.
run_higher_order <- function(dates, vals, y, X, dj = 1/12, nrands = 50, sig.level = 0.95,
                             max.scale = NULL, pwc_method = "complex", compute_wtc = TRUE,
                             progress = NULL) {
  vars <- c(y, X)
  if (length(unique(vars)) < length(vars)) stop("The dependent variable and the predictors must all be different.")
  if (length(X) < 2) stop("Select at least two predictor / control variables.")
  if (!all(vars %in% names(vals))) stop("Selected variables not found in the data.")
  if (pwc_method == "classical" && length(X) != 2)
    stop("The classical partial coherence formula handles exactly one control: select two variables (x and z).")
  spec <- ho_spectra(vals, vars, dj = dj, max.scale = max.scale)
  S <- spec$S; p <- spec$p
  y_idx <- 1L; X_idx <- seq.int(2L, p)
  N <- spec$N
  show_date <- !is.null(dates)
  xaxis <- if (show_date) dates else seq_len(N)
  mk <- function(rsq, phase, signif, type = "pwtc")
    structure(list(coi = spec$coi, rsq = rsq, phase = phase, period = spec$period, scale = spec$scale,
                   dt = spec$dt, t = spec$t, xaxis = xaxis, s0 = spec$s0, dj = spec$dj,
                   mother = "morlet", type = type, signif = if (is.null(signif)) NA else signif),
              class = "biwavelet")

  out <- list(y = y, X = X, vars = vars, N = N, dates = dates, show_date = show_date,
              nrands = nrands, sig.level = sig.level, pwc_method = pwc_method, dj = dj)
  total_steps <- 1 + length(X) + (if (compute_wtc) length(X) else 0)
  step <- 0
  sub_progress <- function(k) function(f) if (is.function(progress)) progress((step + f) / total_steps)

  # multiple wavelet coherence
  mwc_rsq <- ho_mwc_stat(S, y_idx, X_idx)
  mwc_sig <- ho_signif(spec, function(Ss) ho_mwc_stat(Ss, y_idx, X_idx), nrands, sig.level, sub_progress(step))
  out$mwc <- mk(mwc_rsq, matrix(0, nrow(mwc_rsq), ncol(mwc_rsq)), mwc_sig)
  step <- step + 1

  # partial wavelet coherence of y with each predictor, controlling for the others
  out$pwc <- vector("list", length(X)); names(out$pwc) <- X
  for (k in seq_along(X)) {
    x_idx <- X_idx[k]; Z_idx <- setdiff(X_idx, x_idx)
    st <- ho_pwc_stat(S, y_idx, x_idx, Z_idx, pwc_method)
    sg <- ho_signif(spec, function(Ss) ho_pwc_stat(Ss, y_idx, x_idx, Z_idx, pwc_method)$rsq, nrands, sig.level, sub_progress(step))
    out$pwc[[k]] <- mk(st$rsq, st$phase, sg)
    attr(out$pwc[[k]], "controls") <- setdiff(X, X[k])
    step <- step + 1
  }

  # plain bivariate coherence of each pair (same transform and smoothing)
  if (compute_wtc) {
    out$wtc <- vector("list", length(X)); names(out$wtc) <- X
    for (k in seq_along(X)) {
      x_idx <- X_idx[k]
      rsq <- Mod(S[[y_idx]][[x_idx]])^2 / (Re(S[[y_idx]][[y_idx]]) * Re(S[[x_idx]][[x_idx]]))
      sg <- ho_signif(spec, function(Ss) Mod(Ss[[y_idx]][[x_idx]])^2 / (Re(Ss[[y_idx]][[y_idx]]) * Re(Ss[[x_idx]][[x_idx]])),
                      nrands, sig.level, sub_progress(step))
      out$wtc[[k]] <- mk(pmin(rsq, 1), Arg(S[[y_idx]][[x_idx]]), sg, type = "wtc")
      step <- step + 1
    }
  }
  out
}

# ---- plots (reuse the biwavelet-based renderer from partial_wavelets.R) --------
plot_ho_mwc <- function(h, unit = "observation", palette = "blue") {
  plot_biwavelet(h$mwc, h$dates, unit,
                 sprintf("Multiple wavelet coherence: %s | %s", h$y, paste(h$X, collapse = ", ")),
                 plot.phase = FALSE, plot.sig = h$nrands > 0, palette = palette)
}

plot_ho_pwc <- function(h, which = 1, unit = "observation", palette = "blue") {
  x <- h$X[which]; ctrl <- attr(h$pwc[[which]], "controls")
  plot_biwavelet(h$pwc[[which]], h$dates, unit,
                 sprintf("Partial wavelet coherence: %s ~ %s | %s", h$y, x, paste(ctrl, collapse = ", ")),
                 plot.phase = h$pwc_method == "complex", plot.sig = h$nrands > 0, palette = palette)
}

# WTC (top), PWC (middle) and MWC (bottom) for one predictor
plot_ho_compare <- function(h, which = 1, unit = "observation", palette = "blue") {
  x <- h$X[which]; ctrl <- attr(h$pwc[[which]], "controls")
  panels <- list()
  if (!is.null(h$wtc)) panels[[length(panels) + 1]] <- list(obj = h$wtc[[which]], title = sprintf("Wavelet coherence: %s and %s", h$y, x), phase = TRUE)
  panels[[length(panels) + 1]] <- list(obj = h$pwc[[which]], title = sprintf("Partial wavelet coherence: %s ~ %s | %s", h$y, x, paste(ctrl, collapse = ", ")), phase = h$pwc_method == "complex")
  panels[[length(panels) + 1]] <- list(obj = h$mwc, title = sprintf("Multiple wavelet coherence: %s | %s", h$y, paste(h$X, collapse = ", ")), phase = FALSE)
  op <- graphics::par(mfrow = c(length(panels), 1), oma = c(0, 0, 0, 1), mar = c(3.5, 4.5, 3, 5) + 0.1)
  on.exit(graphics::par(op), add = TRUE)
  for (pn in panels) {
    graphics::plot(pn$obj, plot.phase = pn$phase, plot.sig = h$nrands > 0, lty.coi = 1, col.coi = "grey", lwd.coi = 2, lwd.sig = 2,
                   arrow.lwd = 0.03, arrow.len = 0.1, plot.cb = TRUE, xaxt = "n", fill.cols = seq_ramp(64, palette), col.sig = WA_COLORS$ink,
                   ylab = paste0("Period (", unit_plural(unit), ")"), xlab = "", main = pn$title)
    .pw_time_axis(h$dates, length(pn$obj$t))
  }
  invisible()
}

# ---- tables -----------------------------------------------------------------
ho_table <- function(h) {
  st <- function(obj) .pw_period_stats(obj)
  m <- st(h$mwc)
  df <- data.frame(Period = round(m$Period, 2), check.names = FALSE)
  df[[sprintf("MWC %s | %s", h$y, paste(h$X, collapse = ","))]] <- round(m$avg, 3)
  for (k in seq_along(h$X)) {
    x <- h$X[k]; ctrl <- paste(attr(h$pwc[[k]], "controls"), collapse = ",")
    if (!is.null(h$wtc)) df[[sprintf("WTC %s~%s", h$y, x)]] <- round(st(h$wtc[[k]])$avg, 3)
    df[[sprintf("PWC %s~%s|%s", h$y, x, ctrl)]] <- round(st(h$pwc[[k]])$avg, 3)
  }
  if (h$nrands > 0) {
    df[["Sig. share MWC"]] <- round(m$share_sig, 3)
    for (k in seq_along(h$X)) df[[sprintf("Sig. share PWC %s~%s", h$y, h$X[k])]] <- round(st(h$pwc[[k]])$share_sig, 3)
  }
  df
}

# mean partial phase (radians) inside significant PWC regions, by period
ho_phase_table <- function(h) {
  if (h$pwc_method != "complex") return(NULL)
  per <- h$pwc[[1]]$period
  df <- data.frame(Period = round(per, 2))
  for (k in seq_along(h$X)) {
    obj <- h$pwc[[k]]
    inside <- outer(obj$period, obj$coi, FUN = "<=")
    sig <- if (is.matrix(obj$signif)) (obj$rsq / obj$signif) >= 1 else matrix(TRUE, nrow(obj$rsq), ncol(obj$rsq))
    mask <- inside & sig
    ang <- vapply(seq_len(nrow(obj$phase)), function(i) {
      a <- obj$phase[i, mask[i, ]]
      if (!length(a)) NA_real_ else atan2(mean(sin(a)), mean(cos(a)))
    }, numeric(1))
    df[[sprintf("Mean partial phase %s~%s (rad)", h$y, h$X[k])]] <- round(ang, 3)
  }
  df
}

ho_matrix <- function(obj, dates) partial_matrix(obj, dates)
