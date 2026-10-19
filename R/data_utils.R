# ---------------------------------------------------------------------------
# Data loading, cleaning, transformation and summary statistics
# ---------------------------------------------------------------------------

REPLICATION <- list(
  resilience = list(label = "African financial resilience: 12 African markets, asset classes & risk indices (daily, 2017-2025)",
                    file = "data/Africa_Financial_resilience data.csv"),
  africa = list(label = "African stock market indices (weekly, 2011-2017)",
                file = "data/AFRICA_week.csv"),
  crypto = list(label = "Gold & cryptocurrency returns (daily, 2015-2019)",
                file = "data/gold_cryptos9.txt"),
  reit   = list(label = "International REIT returns (daily, 2017-2025)",
                file = "data/REIT.txt"),
  ghana  = list(label = "Ghana cocoa, gold & exchange rate (monthly, 2015-2023)",
                file = "data/Dataset3_Gold+Cocoa+Exrate.txt")
)

replication_choices <- function() {
  stats::setNames(names(REPLICATION), vapply(REPLICATION, `[[`, "", "label"))
}

# Read csv / txt (tab, comma, semicolon) / xlsx into a data.frame.
# Separator and decimal mark are detected automatically for text files.
read_user_file <- function(path, name = path, header = TRUE, sheet = 1) {
  ext <- tolower(tools::file_ext(name))
  if (ext %in% c("xlsx", "xls")) {
    df <- as.data.frame(readxl::read_excel(path, sheet = sheet, col_names = header), stringsAsFactors = FALSE)
  } else {
    lines <- readLines(path, n = 5, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    if (!length(lines)) stop("The file is empty.")
    first <- lines[1]
    sep <- if (grepl("\t", first)) "\t" else if (grepl(";", first)) ";" else if (grepl(",", first)) "," else if (grepl("\\|", first)) "|" else ""
    # decimal comma if the separator is ; (or tab) and numbers look like 1,23
    body <- paste(lines[-1], collapse = "\n")
    dec <- if (sep %in% c(";", "\t", "|") && grepl("[0-9],[0-9]", body) && !grepl("[0-9]\\.[0-9]", body)) "," else "."
    df <- utils::read.table(path, header = header, sep = sep, dec = dec, strip.white = TRUE,
                            stringsAsFactors = FALSE, comment.char = "", quote = "\"",
                            check.names = FALSE, fill = TRUE, na.strings = c("NA", "", "NaN", "#N/A", "n/a", "null", "."))
  }
  if (!header) names(df) <- paste0("V", seq_len(ncol(df)))
  names(df) <- ifelse(is.na(names(df)) | names(df) == "", paste0("V", seq_len(ncol(df))), names(df))
  df
}

DATE_FORMAT_CHOICES <- c("Auto-detect" = "auto", "YYYY-MM-DD" = "%Y-%m-%d", "YYYY/MM/DD" = "%Y/%m/%d",
                         "DD/MM/YYYY" = "%d/%m/%Y", "MM/DD/YYYY" = "%m/%d/%Y", "DD-MM-YYYY" = "%d-%m-%Y",
                         "DD.MM.YYYY" = "%d.%m.%Y", "YYYYMMDD" = "%Y%m%d", "YYYY-MM (monthly)" = "%Y-%m",
                         "Excel serial number" = "excel")

# Try a set of common date formats; return Date vector or NULL
parse_dates <- function(x, format = "auto") {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (format == "excel" || (is.numeric(x) && format == "auto")) {
    xn <- suppressWarnings(as.numeric(x))
    if (all(is.na(xn))) return(NULL)
    # Excel serial dates (1900 system) roughly between 1950 and 2100
    if (mean(xn >= 18000 & xn <= 73000, na.rm = TRUE) > 0.9) return(as.Date(xn, origin = "1899-12-30"))
    # plain year numbers
    if (format == "auto" && mean(xn >= 1800 & xn <= 2200 & xn == round(xn), na.rm = TRUE) > 0.9)
      return(as.Date(paste0(xn, "-01-01")))
    return(NULL)
  }
  x <- trimws(as.character(x))
  x <- sub("[T ]\\d{1,2}:\\d{2}.*$", "", x)      # drop a time part if any
  fmts <- if (format == "auto")
    c("%Y-%m-%d", "%Y/%m/%d", "%d/%m/%Y", "%d-%m-%Y", "%m/%d/%Y", "%d.%m.%Y", "%Y%m%d",
      "%b %d, %Y", "%d %b %Y", "%d-%b-%Y", "%b-%y", "%Y-%m", "%Y/%m", "%Y")
  else format
  for (f in fmts) {
    xf <- if (f %in% c("%Y-%m", "%Y/%m")) paste0(x, if (f == "%Y-%m") "-01" else "/01") else x
    ff <- if (f %in% c("%Y-%m", "%Y/%m")) paste0(f, if (f == "%Y-%m") "-%d" else "/%d") else f
    if (f == "%Y") { xf <- paste0(x, "-01-01"); ff <- "%Y-%m-%d" }
    if (f == "%b-%y") { xf <- paste0("01-", x); ff <- "%d-%b-%y" }
    d <- suppressWarnings(as.Date(xf, format = ff))
    ok <- !is.na(d)
    if (sum(ok) < 0.9 * length(x)) next
    # guard against lenient parsing (e.g. "05/01/2015" read as year 5 with %Y/%m/%d):
    # years must be plausible and numeric-only formats must round-trip exactly
    yr <- as.integer(format(d[ok], "%Y"))
    if (mean(yr >= 1800 & yr <= 2200) < 0.9) next
    if (!grepl("%b", ff) && mean(format(d[ok], ff) == xf[ok]) < 0.9) next
    return(d)
  }
  NULL
}

# Returns list(date = Date | NULL, data = numeric data.frame, note = character)
#   date_col    : name of the date column, "" / NULL for "use first column if it is a date",
#                 or "none" for no date column (index 1..N)
#   date_format : one of DATE_FORMAT_CHOICES
prepare_dataset <- function(df, date_col = NULL, date_format = "auto") {
  if (ncol(df) < 1) stop("The file has no columns.")
  orig_names <- trimws(names(df))
  names(df) <- make.names(orig_names, unique = TRUE)
  note  <- character(0)
  dates <- NULL
  if (is.null(date_col) || !nzchar(date_col)) date_col <- names(df)[1]
  if (identical(date_col, "none")) {
    vals <- df
  } else {
    date_col <- make.names(date_col)
    if (!(date_col %in% names(df))) date_col <- names(df)[1]
    dates <- parse_dates(df[[date_col]], date_format)
    if (!is.null(dates)) {
      vals <- df[setdiff(names(df), date_col)]
    } else {
      vals <- df
      note <- c(note, sprintf("Column '%s' could not be read as a date; observations are indexed 1..N.", date_col))
    }
  }
  if (ncol(vals) < 1) stop("The file needs at least one numeric column besides the date.")
  vals <- as.data.frame(lapply(vals, function(v) {
    if (is.numeric(v)) return(v)
    v <- gsub("[[:space:]]", "", as.character(v))
    v <- gsub("(?<=\\d),(?=\\d{3}(\\D|$))", "", v, perl = TRUE)   # thousands separators 1,234.5
    suppressWarnings(as.numeric(v))
  }), stringsAsFactors = FALSE)
  keep <- vapply(vals, function(v) mean(is.na(v)) < 0.5, logical(1))
  if (any(!keep)) note <- c(note, paste("Dropped non-numeric column(s):", paste(names(vals)[!keep], collapse = ", ")))
  vals <- vals[keep]
  if (ncol(vals) < 1) stop("No numeric columns found.")
  if (!is.null(dates)) {
    ok <- !is.na(dates)
    vals <- vals[ok, , drop = FALSE]; dates <- dates[ok]
    ord <- order(dates)
    if (is.unsorted(dates)) note <- c(note, "Observations were re-ordered chronologically.")
    vals <- vals[ord, , drop = FALSE]; dates <- dates[ord]
  }
  cc <- stats::complete.cases(vals)
  if (any(!cc)) {
    note <- c(note, sprintf("Removed %d row(s) with missing values.", sum(!cc)))
    vals <- vals[cc, , drop = FALSE]
    if (!is.null(dates)) dates <- dates[cc]
  }
  rownames(vals) <- NULL
  list(date = dates, data = vals, note = note)
}

TRANSFORM_CHOICES <- c("None" = "none",
                       "nth differences" = "diff",
                       "nth percentage changes" = "pct",
                       "nth log returns (x100)" = "logret",
                       "Standardized nth differences" = "std_diff")

VOLATILITY_CHOICES <- c("None" = "none", "Absolute values" = "abs", "Squared values" = "sq")

transform_label <- function(transformation, n, volatility) {
  tl <- names(TRANSFORM_CHOICES)[TRANSFORM_CHOICES == transformation]
  vl <- names(VOLATILITY_CHOICES)[VOLATILITY_CHOICES == volatility]
  if (transformation == "none" && volatility == "none") return("No transformation")
  out <- if (transformation == "none") "Levels" else sprintf("%s (n = %d)", tl, n)
  if (volatility != "none") out <- paste(out, "->", vl)
  out
}

# Apply transformation and volatility proxy. Returns list(date, data)
transform_data <- function(dates, vals, transformation = "none", n = 1, volatility = "none") {
  n <- max(1L, as.integer(n))
  N <- nrow(vals)
  if (transformation != "none" && N <= n + 10)
    stop("Not enough observations for the chosen lag operator.")
  out <- vals
  if (transformation != "none") {
    out <- as.data.frame(lapply(vals, function(x) {
      switch(transformation,
             diff     = diff(x, lag = n),
             pct      = 100 * (x[(n + 1):N] / x[1:(N - n)] - 1),
             logret   = { if (any(x <= 0)) stop("Log returns need strictly positive data."); 100 * diff(log(x), lag = n) },
             std_diff = as.numeric(scale(diff(x, lag = n))))
    }))
    names(out) <- names(vals)
    if (!is.null(dates)) dates <- dates[-(1:n)]
  }
  if (volatility == "abs") out <- as.data.frame(lapply(out, abs))
  if (volatility == "sq")  out <- as.data.frame(lapply(out, function(x) x^2))
  names(out) <- names(vals)
  note <- character(0)
  bad <- !is.finite(as.matrix(out))
  if (any(bad)) {
    keep <- rowSums(bad) == 0
    cols <- names(out)[colSums(bad) > 0]
    note <- sprintf("%d row(s) removed because the transformation produced non-finite values in: %s (e.g. percentage change of a zero).",
                    sum(!keep), paste(cols, collapse = ", "))
    out <- out[keep, , drop = FALSE]
    if (!is.null(dates)) dates <- dates[keep]
  }
  rownames(out) <- NULL
  list(date = dates, data = out, note = note)
}

# ---------------------------------------------------------------------------
# Summary statistics (in the spirit of the connectedness app)
# ---------------------------------------------------------------------------
skewness <- function(x) { x <- x - mean(x); mean(x^3) / (mean(x^2))^1.5 }
ex_kurtosis <- function(x) { x <- x - mean(x); mean(x^4) / (mean(x^2))^2 - 3 }

jarque_bera <- function(x) {
  n <- length(x); s <- skewness(x); k <- ex_kurtosis(x)
  stat <- n / 6 * (s^2 + k^2 / 4)
  c(stat = stat, p = 1 - stats::pchisq(stat, df = 2))
}

summary_stats <- function(vals) {
  rows <- lapply(names(vals), function(v) {
    x  <- vals[[v]]
    jb <- jarque_bera(x)
    q  <- stats::Box.test(x, lag = 20, type = "Ljung-Box")
    q2 <- stats::Box.test(x^2, lag = 20, type = "Ljung-Box")
    data.frame(Variable = v, N = length(x), Mean = mean(x), SD = stats::sd(x),
               Min = min(x), Max = max(x), Skewness = skewness(x),
               `Ex.Kurtosis` = ex_kurtosis(x), `JB` = jb[["stat"]], `JB p` = jb[["p"]],
               `Q(20)` = unname(q$statistic), `Q(20) p` = q$p.value,
               `Q2(20)` = unname(q2$statistic), `Q2(20) p` = q2$p.value,
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Long format helper for ggplot facets
to_long <- function(dates, vals) {
  idx <- if (is.null(dates)) seq_len(nrow(vals)) else dates
  do.call(rbind, lapply(names(vals), function(v)
    data.frame(time = idx, variable = v, value = vals[[v]], stringsAsFactors = FALSE)))
}

plot_series <- function(dates, vals, title = NULL) {
  long <- to_long(dates, vals)
  long$variable <- factor(long$variable, levels = names(vals))
  ggplot2::ggplot(long, ggplot2::aes(x = time, y = value)) +
    ggplot2::geom_line(colour = "#1f6f8b", linewidth = 0.4) +
    ggplot2::facet_wrap(~variable, scales = "free_y", ncol = 2) +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    app_theme()
}

app_theme <- function() {
  ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "right")
}

unit_plural <- function(unit) {
  switch(unit, day = "days", week = "weeks", month = "months", quarter = "quarters",
         year = "years", observation = "observations", unit)
}
