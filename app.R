# ===========================================================================
# Wavelet Approach - Bivariate & Multiple Wavelet Analysis
# Author : Peterson Owusu Junior (powusujunior@gmail.com)
# A Shiny app in the spirit of the Connectedness Approach app
# (https://davidgabauer.shinyapps.io/connectedness_approach/)
#
# Engines:
#   * Bivariate  : continuous Morlet wavelet transform, cross-wavelet power,
#                  wavelet coherence and phase difference  (WaveletComp)
#   * Multiple   : MODWT wavelet correlation, cross-correlation, wavelet
#                  multiple (cross-)correlation and local multiple correlation
#                  (waveslim / wavemulcor, following Fernandez-Macho 2012 and
#                  Polanco-Martinez's W2CWM2C)
# ===========================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(ggplot2)
  library(WaveletComp)
  library(waveslim)
  library(wavemulcor)
  library(shinycssloaders)
})

if (!exists("prepare_dataset")) {
  for (f in list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)) source(f, local = FALSE)
}

options(shiny.maxRequestSize = 30 * 1024^2, spinner.type = 6, spinner.color = "#1f6f8b")

APP_AUTHOR <- list(
  name        = "Peterson Owusu Junior",
  affiliation = "University of Cape Coast, Ghana",
  email       = "peterson.owusu@ucc.edu.gh",          # institutional (used for the footer link)
  emails      = c("peterson.owusu@ucc.edu.gh", "powusujunior@gmail.com"),
  links       = c("Google Scholar" = "https://scholar.google.com/citations?user=r6ZueU4AAAAJ&hl=en",
                  "Scopus"         = "https://www.scopus.com/authid/detail.uri?authorId=58833238800",
                  "LinkedIn"       = "https://www.linkedin.com/in/peterson-owusu-junior-phd-80858334/",
                  "ResearchGate"   = "https://www.researchgate.net/profile/Peterson-Owusu-Junior",
                  "IDEAS/RePEc"    = "https://ideas.repec.org/e/pow25.html")
)
APP_AUTHOR$links <- APP_AUTHOR$links[nzchar(APP_AUTHOR$links)]   # drop links not yet filled in

UNIT_CHOICES <- c("Daily" = "day", "Weekly" = "week", "Monthly" = "month",
                  "Quarterly" = "quarter", "Yearly" = "year", "Unspecified" = "observation")

BIV_TABS  <- c("tab_wt", "tab_xwt", "tab_coh", "tab_avg", "tab_phase", "tab_phimg", "tab_rec", "tab_bivtab")
MULT_TABS <- c("tab_modwt", "tab_wc", "tab_wcc", "tab_wmc", "tab_wmcc", "tab_wlmc", "tab_multtab")

plot_card <- function(id, height = "560px", ...) {
  card(full_screen = TRUE, card_body(..., withSpinner(plotOutput(id, height = height))))
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- page_sidebar(
  title = tags$span(tags$strong("Wavelet Approach"),
                    tags$small(class = "text-muted ms-2", "Bivariate & multiple wavelet analysis")),
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#1f6f8b",
                   base_font = font_google("Source Sans 3", local = FALSE),
                   heading_font = font_google("Source Sans 3", local = FALSE)),
  tags$head(tags$link(rel = "stylesheet", href = "custom.css")),
  fillable = FALSE,

  sidebar = sidebar(
    width = 370, open = "always",
    accordion(
      multiple = TRUE, open = c("Data", "Wavelet analysis"),
      accordion_panel(
        "Data", icon = icon("database"),
        radioButtons("source", NULL, c("Own dataset" = "own", "Example dataset" = "rep"),
                     selected = "own", inline = TRUE),
        conditionalPanel(
          "input.source == 'own'",
          fileInput("file", "Upload csv / txt / xlsx file", accept = c(".csv", ".txt", ".tsv", ".xlsx", ".xls")),
          helpText("One column per series, observations in rows. A date column (e.g. YYYY-MM-DD) is recommended; ",
                   "separators, decimal marks and most date formats are detected automatically."),
          checkboxInput("has_header", "First row contains column names", value = TRUE),
          layout_columns(col_widths = c(6, 6),
                         selectInput("date_col", "Date column", choices = c("First column" = "")),
                         selectInput("date_format", "Date format", DATE_FORMAT_CHOICES, selected = "auto"))
        ),
        conditionalPanel(
          "input.source == 'rep'",
          selectInput("rep_dataset", "Example dataset", choices = replication_choices(), selected = "africa"),
          helpText("Built-in examples to try the app; upload your own file under 'Own dataset'.")
        ),
        selectInput("unit", "Data frequency (labels only)", UNIT_CHOICES, selected = "week")
      ),
      accordion_panel(
        "Preprocessing", icon = icon("sliders"),
        selectInput("transformation", "Transformation", TRANSFORM_CHOICES, selected = "none"),
        numericInput("lag_n", "Lag operator (n)", value = 1, min = 1, step = 1),
        selectInput("volatility", "Volatility proxy", VOLATILITY_CHOICES, selected = "none")
      ),
      accordion_panel(
        "Wavelet analysis", icon = icon("wave-square"),
        radioButtons("mode", "Approach",
                     c("Bivariate wavelets (continuous Morlet)" = "biv",
                       "Multiple wavelets (MODWT)" = "mult"), selected = "biv"),
        conditionalPanel(
          "input.mode == 'biv'",
          layout_columns(col_widths = c(6, 6),
                         selectInput("var1", "Variable x", choices = NULL),
                         selectInput("var2", "Variable y", choices = NULL)),
          selectInput("dj", "Frequency resolution (dj)", DJ_CHOICES, selected = "1/20"),
          layout_columns(col_widths = c(6, 6),
                         numericInput("lower_period", "Lower period", value = 2, min = 2, step = 1),
                         numericInput("upper_period", "Upper period (blank = N/3)", value = NA, min = 4, step = 1)),
          numericInput("loess_span", "Loess span for detrending (0 = none)", value = 0, min = 0, max = 1, step = 0.05),
          checkboxInput("make_pval", "Significance test (Monte-Carlo)", value = TRUE),
          conditionalPanel(
            "input.make_pval",
            layout_columns(col_widths = c(7, 5),
                           selectInput("sim_method", "Null hypothesis", SIM_METHODS, selected = "white.noise"),
                           numericInput("n_sim", "Simulations", value = 100, min = 10, max = 1000, step = 10)),
            numericInput("siglvl", "Significance level", value = 0.05, min = 0.01, max = 0.2, step = 0.01)
          ),
          layout_columns(col_widths = c(6, 6),
                         selectInput("window_type_t", "Smoothing window (time)", WINDOW_TYPES, selected = 1),
                         selectInput("window_type_s", "Smoothing window (scale)", WINDOW_TYPES, selected = 1)),
          layout_columns(col_widths = c(6, 6),
                         numericInput("window_size_t", "Window size (time)", value = 5, min = 1, step = 1),
                         numericInput("window_size_s", "Window size (scale)", value = 0.25, min = 0.05, step = 0.05))
        ),
        conditionalPanel(
          "input.mode == 'mult'",
          selectizeInput("vars", "Variables (2 or more)", choices = NULL, multiple = TRUE,
                         options = list(plugins = list("remove_button"))),
          selectInput("wf", "Wavelet filter", FILTER_CHOICES, selected = "la8"),
          layout_columns(col_widths = c(6, 6),
                         numericInput("J", "Decomposition levels (J)", value = 5, min = 1, max = 12, step = 1),
                         numericInput("lmax", "Max lag (cross-corr.)", value = 10, min = 1, max = 200, step = 1)),
          uiOutput("j_hint"),
          layout_columns(col_widths = c(6, 6),
                         selectInput("wcc_a", "Cross-corr. pair: x", choices = NULL),
                         selectInput("wcc_b", "Cross-corr. pair: y", choices = NULL)),
          layout_columns(col_widths = c(6, 6),
                         numericInput("wlmc_M", "Local window M (blank = N/8)", value = NA, min = 5, step = 1),
                         selectInput("wlmc_window", "Local window type", WINDOW_CHOICES, selected = "gauss"))
        )
      )
    ),
    actionButton("run", "Start analysis!", class = "btn-primary w-100 mt-2", icon = icon("play")),
    downloadButton("download", "Download results", class = "btn-outline-secondary w-100 mt-2"),
    uiOutput("status")
  ),

  uiOutput("data_info"),

  navset_card_underline(
    id = "tabs",
    nav_panel("Data table", value = "tab_data", withSpinner(DTOutput("data_table"))),
    nav_panel("Raw data", value = "tab_raw", withSpinner(plotOutput("raw_plot", height = "620px"))),
    nav_panel("Transformed data", value = "tab_trans", withSpinner(plotOutput("trans_plot", height = "620px"))),
    nav_panel("Summary statistics", value = "tab_stats", withSpinner(DTOutput("stats_table"))),

    # ---- bivariate --------------------------------------------------------
    nav_panel("Wavelet power", value = "tab_wt",
              uiOutput("wt_series_ui"),
              withSpinner(plotOutput("wt_plot", height = "560px"))),
    nav_panel("Cross-wavelet power", value = "tab_xwt", withSpinner(plotOutput("xwt_plot", height = "560px"))),
    nav_panel("Wavelet coherence", value = "tab_coh", withSpinner(plotOutput("coh_plot", height = "560px"))),
    nav_panel("Average power & coherence", value = "tab_avg",
              layout_columns(col_widths = c(6, 6),
                             withSpinner(plotOutput("avg_wp_plot", height = "480px")),
                             withSpinner(plotOutput("avg_wc_plot", height = "480px")))),
    nav_panel("Phase difference", value = "tab_phase",
              uiOutput("phase_controls"),
              withSpinner(plotOutput("phase_plot", height = "520px")),
              tags$p(class = "text-muted small mt-2",
                     "Phase difference = phase(x) - phase(y). Values in (0, π/2): in phase, x leads; ",
                     "(-π/2, 0): in phase, y leads; (π/2, π): anti-phase, y leads; (-π, -π/2): anti-phase, x leads.")),
    nav_panel("Phase-difference image", value = "tab_phimg", withSpinner(plotOutput("phimg_plot", height = "560px"))),
    nav_panel("Reconstruction", value = "tab_rec",
              uiOutput("rec_controls"),
              withSpinner(plotOutput("rec_plot", height = "520px"))),
    nav_panel("Tables", value = "tab_bivtab",
              h5("Average cross-wavelet power and coherence by period"),
              withSpinner(DTOutput("coh_table")),
              h5(class = "mt-4", "Phase summary inside significant coherence regions"),
              withSpinner(DTOutput("phase_table"))),

    # ---- multiple ---------------------------------------------------------
    nav_panel("MODWT decomposition", value = "tab_modwt",
              uiOutput("modwt_controls"),
              withSpinner(plotOutput("modwt_plot", height = "760px"))),
    nav_panel("Wavelet correlation", value = "tab_wc",
              uiOutput("wc_controls"),
              withSpinner(plotOutput("wc_plot", height = "640px"))),
    nav_panel("Wavelet cross-correlation", value = "tab_wcc",
              withSpinner(plotOutput("wcc_plot", height = "520px")),
              withSpinner(plotOutput("wcc_lines", height = "520px"))),
    nav_panel("Multiple correlation", value = "tab_wmc", withSpinner(plotOutput("wmc_plot", height = "600px"))),
    nav_panel("Multiple cross-correlation", value = "tab_wmcc",
              withSpinner(plotOutput("wmcc_plot", height = "520px")),
              withSpinner(plotOutput("wmcc_lines", height = "520px"))),
    nav_panel("Local multiple correlation", value = "tab_wlmc",
              withSpinner(plotOutput("wlmc_plot", height = "520px")),
              withSpinner(plotOutput("wlmc_lines", height = "620px"))),
    nav_panel("Tables", value = "tab_multtab",
              h5("Wavelet multiple correlation"), withSpinner(DTOutput("wmc_table")),
              h5(class = "mt-4", "Pairwise wavelet correlation"), withSpinner(DTOutput("wc_table")),
              h5(class = "mt-4", "Wavelet multiple cross-correlation"), withSpinner(DTOutput("wmcc_table"))),

    nav_spacer(),
    nav_panel("About", value = "tab_about", icon = icon("circle-info"), uiOutput("about"))
  ),

  tags$footer(class = "text-muted small mt-3",
              HTML(sprintf('<strong>Wavelet Approach</strong> by <a href="mailto:%s">%s</a>, %s. ',
                           APP_AUTHOR$email, APP_AUTHOR$name, APP_AUTHOR$affiliation)),
              "Built with R Shiny, WaveletComp, waveslim and wavemulcor. ",
              "Layout inspired by the Connectedness Approach app by David Gabauer.")
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  # ---- data --------------------------------------------------------------
  # uploaded file, read once per file / header setting
  file_df <- reactive({
    req(input$file)
    tryCatch(read_user_file(input$file$datapath, input$file$name, header = isTRUE(input$has_header)),
             error = function(e) { showNotification(paste("Could not read file:", conditionMessage(e)), type = "error", duration = 10); NULL })
  })

  observeEvent(file_df(), {
    df <- file_df(); req(df)
    cols <- trimws(names(df))
    updateSelectInput(session, "date_col",
                      choices = c(stats::setNames(cols, cols), "No date column (use index)" = "none"),
                      selected = cols[1])
  })

  raw <- reactive({
    if (input$source == "rep") {
      f <- REPLICATION[[input$rep_dataset]]$file
      prepare_dataset(read_user_file(f, f))
    } else {
      df <- file_df(); req(df)
      tryCatch(prepare_dataset(df, input$date_col, input$date_format %||% "auto"),
               error = function(e) { showNotification(conditionMessage(e), type = "error", duration = 10); NULL })
    }
  })

  observeEvent(input$rep_dataset, {
    unit <- switch(input$rep_dataset, africa = "week", crypto = "day", reit = "day", ghana = "month", "observation")
    updateSelectInput(session, "unit", selected = unit)
  })

  observeEvent(raw(), {
    r <- raw(); req(r)
    vars <- names(r$data); N <- nrow(r$data)
    updateSelectInput(session, "var1", choices = vars, selected = vars[1])
    updateSelectInput(session, "var2", choices = vars, selected = vars[min(2, length(vars))])
    updateSelectizeInput(session, "vars", choices = vars, selected = vars[seq_len(min(7, length(vars)))])
    updateSelectInput(session, "wcc_a", choices = vars, selected = vars[1])
    updateSelectInput(session, "wcc_b", choices = vars, selected = vars[min(2, length(vars))])
    updateNumericInput(session, "J", value = min(isolate(input$J), max_levels(N, isolate(input$wf))))
  })

  output$j_hint <- renderUI({
    r <- raw(); req(r)
    helpText(sprintf("Max. levels for N = %d and %s filter: J = %d (boundary coefficients are removed).",
                     nrow(r$data), input$wf, max_levels(nrow(r$data), input$wf)))
  })

  # ---- tab visibility ----------------------------------------------------
  observeEvent(input$mode, {
    req(input$mode)
    for (t in BIV_TABS)  if (input$mode == "biv")  nav_show("tabs", t) else nav_hide("tabs", t)
    for (t in MULT_TABS) if (input$mode == "mult") nav_show("tabs", t) else nav_hide("tabs", t)
  }, ignoreNULL = FALSE)

  output$data_info <- renderUI({
    if (input$source == "own" && is.null(input$file))
      return(div(class = "alert alert-info small py-2 mb-3", icon("upload"),
                 " Upload a csv, txt or xlsx file in the sidebar (or switch to an example dataset), choose the settings, then click ",
                 tags$b("Start analysis!"), "."))
    r <- raw(); req(r)
    rng <- if (is.null(r$date)) "no dates" else paste(format(range(r$date)), collapse = " to ")
    div(class = "alert alert-light border small py-2 mb-3",
        tags$strong(sprintf("%d observations, %d variables", nrow(r$data), ncol(r$data))),
        " | ", rng, " | ", paste(names(r$data), collapse = ", "),
        if (length(r$note)) tags$div(class = "text-warning mt-1", paste(r$note, collapse = " ")))
  })

  # ---- analysis ----------------------------------------------------------
  res <- eventReactive(input$run, {
    r <- raw()
    validate(need(!is.null(r), "Please load a dataset first."))
    t0 <- Sys.time()
    withProgress(message = "Running wavelet analysis", value = 0.05, {
      prep <- tryCatch(transform_data(r$date, r$data, input$transformation, input$lag_n, input$volatility),
                       error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
      validate(need(!is.null(prep), "Transformation failed."))
      out <- list(mode = input$mode, unit = input$unit, raw = r, date = prep$date, data = prep$data,
                  transform = transform_label(input$transformation, input$lag_n, input$volatility),
                  stats = summary_stats(prep$data))
      N <- nrow(prep$data)

      if (input$mode == "biv") {
        validate(need(input$var1 != input$var2, "Please select two different variables."))
        incProgress(0.1, detail = sprintf("Wavelet coherence (%d simulations)", if (input$make_pval) input$n_sim else 0))
        dj <- eval(parse(text = input$dj))
        out$biv <- run_coherency(prep$date, prep$data, c(input$var1, input$var2), dj = dj,
                                 lowerPeriod = max(2, input$lower_period), upperPeriod = input$upper_period,
                                 loess.span = input$loess_span, make.pval = input$make_pval,
                                 method = input$sim_method, n.sim = max(10, input$n_sim),
                                 window.type.t = as.numeric(input$window_type_t),
                                 window.type.s = as.numeric(input$window_type_s),
                                 window.size.t = input$window_size_t, window.size.s = input$window_size_s)
        out$siglvl <- input$siglvl
      } else {
        vars <- intersect(input$vars, names(prep$data))
        validate(need(length(vars) >= 2, "Please select at least two variables."))
        Jmax <- max_levels(N, input$wf)
        J <- min(max(1, as.integer(input$J)), Jmax)
        if (J < input$J) showNotification(sprintf("J reduced to %d (maximum for this sample and filter).", J), type = "warning")
        lmax <- min(max(1, as.integer(input$lmax)), floor(N / 4))
        incProgress(0.1, detail = "MODWT decomposition")
        L <- modwt_list(prep$data[vars], input$wf, J)
        incProgress(0.15, detail = "Pairwise wavelet correlation")
        m <- list(vars = vars, J = J, wf = input$wf, lmax = lmax, L = L,
                  wc = pairwise_wc(L, N, J))
        incProgress(0.15, detail = "Wavelet multiple correlation")
        m$wmc  <- wmc_calc(L, N, J)
        incProgress(0.15, detail = "Wavelet multiple cross-correlation")
        m$wmcc <- wmcc_calc(L, lmax, J)
        a <- if (input$wcc_a %in% vars) input$wcc_a else vars[1]
        b <- if (input$wcc_b %in% vars && input$wcc_b != a) input$wcc_b else setdiff(vars, a)[1]
        m$wcc_pair <- c(a, b)
        m$wcc <- pair_wcc(L[[a]], L[[b]], lmax, N, J, a, b)
        incProgress(0.15, detail = "Local multiple correlation")
        M <- if (is.na(input$wlmc_M)) max(15, round(N / 8)) else max(5, as.integer(input$wlmc_M))
        m$M <- M
        m$wlmc <- tryCatch(wlmc_calc(L, M, input$wlmc_window, J, prep$date),
                           error = function(e) { showNotification(paste("WLMC failed:", conditionMessage(e)), type = "warning"); NULL })
        out$mult <- m
      }
      incProgress(0.2, detail = "Done")
    })
    out$elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    out
  })

  output$status <- renderUI({
    r <- res()
    div(class = "small text-success mt-2", icon("check"),
        sprintf(" Analysis finished in %s s (%s).", r$elapsed, r$transform))
  })

  # ---- common outputs ----------------------------------------------------
  output$data_table <- renderDT({
    r <- raw(); req(r)
    df <- if (!is.null(r$date)) cbind(Date = format(r$date), r$data) else r$data
    datatable(df, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "compact stripe") |>
      formatRound(names(r$data), 4)
  })
  output$raw_plot <- renderPlot({ r <- raw(); req(r); plot_series(r$date, r$data, "Raw data") }, res = 100)
  output$trans_plot <- renderPlot({ r <- res(); plot_series(r$date, r$data, r$transform) }, res = 100)
  output$stats_table <- renderDT({
    r <- res()
    datatable(r$stats, rownames = FALSE, options = list(pageLength = 25, scrollX = TRUE, dom = "t"), class = "compact stripe") |>
      formatRound(setdiff(names(r$stats), c("Variable", "N")), 3)
  })

  # ---- bivariate outputs -------------------------------------------------
  biv <- reactive({ r <- res(); validate(need(r$mode == "biv", "Run a bivariate analysis to see this tab.")); r })

  output$wt_series_ui <- renderUI({
    r <- biv()
    radioButtons("wt_series", NULL, choices = r$biv$pair, selected = r$biv$pair[1], inline = TRUE)
  })
  output$wt_plot <- renderPlot({
    r <- biv(); s <- input$wt_series; if (is.null(s) || !(s %in% r$biv$pair)) s <- r$biv$pair[1]
    plot_wt_power(r$biv, s, r$unit, r$siglvl)
  }, res = 100)
  output$xwt_plot <- renderPlot({ r <- biv(); plot_cross_power(r$biv, r$unit, r$siglvl) }, res = 100)
  output$coh_plot <- renderPlot({ r <- biv(); plot_coherence(r$biv, r$unit, r$siglvl) }, res = 100)
  output$avg_wp_plot <- renderPlot({ r <- biv(); plot_avg(r$biv, "wp", r$unit) }, res = 100)
  output$avg_wc_plot <- renderPlot({ r <- biv(); plot_avg(r$biv, "wc", r$unit) }, res = 100)

  output$phase_controls <- renderUI({
    r <- biv(); per <- r$biv$wc$Period
    lo <- min(per); hi <- max(per)
    layout_columns(col_widths = c(8, 4),
                   sliderInput("phase_range", sprintf("Period band (%s)", unit_plural(r$unit)),
                               min = floor(lo), max = ceiling(hi),
                               value = c(max(floor(lo), round(hi / 8)), min(ceiling(hi), round(hi / 2))), step = 1, width = "100%"),
                   selectInput("phase_sig", "Significance criterion",
                               c("Wavelet coherence" = "wc", "Cross-wavelet power" = "wp"), selected = "wc"))
  })
  output$phase_plot <- renderPlot({
    r <- biv(); req(input$phase_range)
    plot_phase(r$biv, input$phase_range[1], input$phase_range[2], r$siglvl, input$phase_sig %||% "wc")
  }, res = 100)
  output$phimg_plot <- renderPlot({ r <- biv(); plot_phasediff_image(r$biv, r$unit, r$siglvl) }, res = 100)

  output$rec_controls <- renderUI({
    r <- biv()
    layout_columns(col_widths = c(6, 6),
                   radioButtons("rec_series", "Series", choices = r$biv$pair, selected = r$biv$pair[1], inline = TRUE),
                   radioButtons("rec_type", "Display", c("Reconstruction" = "rec", "Significant waves" = "waves"), inline = TRUE))
  })
  output$rec_plot <- renderPlot({
    r <- biv(); s <- input$rec_series; if (is.null(s) || !(s %in% r$biv$pair)) s <- r$biv$pair[1]
    plot_reconstruct(r$biv, s, r$siglvl, waves = identical(input$rec_type, "waves"))
  }, res = 100)

  output$coh_table <- renderDT({
    r <- biv(); df <- coherence_table(r$biv)
    datatable(df, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "compact stripe") |>
      formatRound(names(df), 3)
  })
  output$phase_table <- renderDT({
    r <- biv()
    datatable(phase_summary(r$biv, r$siglvl), rownames = FALSE, options = list(pageLength = 15), class = "compact stripe") |>
      formatRound("Period", 2)
  })

  # ---- multiple outputs --------------------------------------------------
  mult <- reactive({ r <- res(); validate(need(r$mode == "mult", "Run a multiple-wavelet analysis to see this tab.")); r })

  output$modwt_controls <- renderUI({
    r <- mult()
    selectInput("modwt_var", "Variable", choices = r$mult$vars, selected = r$mult$vars[1], width = "300px")
  })
  output$modwt_plot <- renderPlot({
    r <- mult(); v <- input$modwt_var; if (is.null(v) || !(v %in% r$mult$vars)) v <- r$mult$vars[1]
    plot_modwt(modwt_decomposition(r$data[[v]], r$date, r$mult$wf, r$mult$J),
               sprintf("MODWT decomposition of %s (%s filter, J = %d)", v, r$mult$wf, r$mult$J))
  }, res = 100)
  output$wc_controls <- renderUI({
    r <- mult()
    radioButtons("wc_layout", NULL, inline = TRUE,
                 choices = c("Sorted pairs (W2CWM2C style)" = "pairs", "Correlation matrix per level" = "matrix"),
                 selected = if (length(r$mult$vars) > 7) "matrix" else "pairs")
  })
  output$wc_plot <- renderPlot({
    r <- mult(); m <- r$mult
    layout <- input$wc_layout %||% (if (length(m$vars) > 7) "matrix" else "pairs")
    if (layout == "matrix") plot_wc_matrix(m$wc, m$vars, m$J, unit_plural(r$unit))
    else plot_pairwise_wc(m$wc, m$J, unit_plural(r$unit))
  }, res = 100)
  output$wcc_plot <- renderPlot({
    r <- mult(); m <- r$mult
    plot_lag_heatmap(m$wcc, "wcc", m$J, unit_plural(r$unit),
                     sprintf("Wavelet cross-correlation: %s vs %s", m$wcc_pair[1], m$wcc_pair[2]),
                     lag_unit = unit_plural(r$unit))
  }, res = 100)
  output$wcc_lines <- renderPlot({
    r <- mult(); m <- r$mult
    plot_lag_lines(m$wcc, "wcc", m$J, unit_plural(r$unit),
                   sprintf("Wavelet cross-correlation by level: %s vs %s", m$wcc_pair[1], m$wcc_pair[2]),
                   lag_unit = unit_plural(r$unit))
  }, res = 100)
  output$wmc_plot <- renderPlot({
    r <- mult(); plot_wmc(r$mult$wmc, unit_plural(r$unit),
                          sprintf("Wavelet multiple correlation (%s)", paste(r$mult$vars, collapse = ", ")))
  }, res = 100)
  output$wmcc_plot <- renderPlot({
    r <- mult(); m <- r$mult
    plot_lag_heatmap(m$wmcc, "wmcc", m$J, unit_plural(r$unit), "Wavelet multiple cross-correlation",
                     side_labels = attr(m$wmcc, "max_var"), lag_unit = unit_plural(r$unit))
  }, res = 100)
  output$wmcc_lines <- renderPlot({
    r <- mult(); m <- r$mult
    plot_lag_lines(m$wmcc, "wmcc", m$J, unit_plural(r$unit), "Wavelet multiple cross-correlation by level",
                   lag_unit = unit_plural(r$unit))
  }, res = 100)
  output$wlmc_plot <- renderPlot({
    r <- mult(); validate(need(!is.null(r$mult$wlmc), "Local multiple correlation could not be computed."))
    plot_wlmc(r$mult$wlmc, r$mult$J, unit_plural(r$unit),
              sprintf("Wavelet local multiple correlation (window M = %d)", r$mult$M))
  }, res = 100)
  output$wlmc_lines <- renderPlot({
    r <- mult(); validate(need(!is.null(r$mult$wlmc), ""))
    plot_wlmc_lines(r$mult$wlmc, r$mult$J, unit_plural(r$unit))
  }, res = 100)

  output$wmc_table <- renderDT({
    r <- mult()
    datatable(r$mult$wmc, rownames = FALSE, options = list(dom = "t"), class = "compact stripe") |>
      formatRound(c("wmc", "lower", "upper"), 3)
  })
  output$wc_table <- renderDT({
    r <- mult(); df <- r$mult$wc[, c("pair", "level", "scale", "wavecor", "lower", "upper", "significant")]
    datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "compact stripe") |>
      formatRound(c("wavecor", "lower", "upper"), 3)
  })
  output$wmcc_table <- renderDT({
    r <- mult(); df <- r$mult$wmcc
    datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE), class = "compact stripe") |>
      formatRound(c("wmcc", "lower", "upper"), 3)
  })

  # ---- download ----------------------------------------------------------
  output$download <- downloadHandler(
    filename = function() paste0("wavelet_results_", format(Sys.time(), "%Y%m%d_%H%M"), ".zip"),
    content = function(file) {
      r <- res()
      td <- file.path(tempdir(), paste0("wavelet_", as.integer(Sys.time())))
      dir.create(td, showWarnings = FALSE)
      old <- setwd(td); on.exit(setwd(old), add = TRUE)
      write_results(r, td)
      zip::zip(file, files = list.files(td), root = td)
    }
  )

  output$about <- renderUI(about_ui())
}

# ---------------------------------------------------------------------------
# Results export
# ---------------------------------------------------------------------------
write_results <- function(r, dir) {
  save_png <- function(name, expr, w = 11, h = 6.5) {
    grDevices::png(file.path(dir, name), width = w, height = h, units = "in", res = 130)
    on.exit(grDevices::dev.off(), add = TRUE)
    force(expr)
  }
  save_gg <- function(name, p, w = 11, h = 6.5) ggplot2::ggsave(file.path(dir, name), p, width = w, height = h, dpi = 130)
  dat <- if (!is.null(r$date)) cbind(Date = format(r$date), r$data) else r$data
  utils::write.csv(dat, file.path(dir, "transformed_data.csv"), row.names = FALSE)
  utils::write.csv(r$stats, file.path(dir, "summary_statistics.csv"), row.names = FALSE)
  save_gg("transformed_data.png", plot_series(r$date, r$data, r$transform))
  u <- unit_plural(r$unit)
  if (r$mode == "biv") {
    b <- r$biv
    utils::write.csv(coherence_table(b), file.path(dir, "average_power_coherence.csv"), row.names = FALSE)
    utils::write.csv(phase_summary(b, r$siglvl), file.path(dir, "phase_summary.csv"), row.names = FALSE)
    cm <- b$wc$Coherence; dimnames(cm) <- list(round(b$wc$Period, 3), if (!is.null(r$date)) format(r$date) else seq_len(ncol(cm)))
    utils::write.csv(cm, file.path(dir, "coherence_matrix_period_by_time.csv"))
    am <- b$wc$Angle; dimnames(am) <- dimnames(cm)
    utils::write.csv(am, file.path(dir, "phase_difference_matrix_period_by_time.csv"))
    save_png("wavelet_power_x.png", plot_wt_power(b, b$pair[1], r$unit, r$siglvl))
    save_png("wavelet_power_y.png", plot_wt_power(b, b$pair[2], r$unit, r$siglvl))
    save_png("cross_wavelet_power.png", plot_cross_power(b, r$unit, r$siglvl))
    save_png("wavelet_coherence.png", plot_coherence(b, r$unit, r$siglvl))
    save_png("average_coherence.png", plot_avg(b, "wc", r$unit))
    save_png("phase_difference_image.png", plot_phasediff_image(b, r$unit, r$siglvl))
  } else {
    m <- r$mult
    utils::write.csv(m$wc, file.path(dir, "pairwise_wavelet_correlation.csv"), row.names = FALSE)
    utils::write.csv(m$wcc, file.path(dir, sprintf("wavelet_cross_correlation_%s_%s.csv", m$wcc_pair[1], m$wcc_pair[2])), row.names = FALSE)
    utils::write.csv(m$wmc, file.path(dir, "wavelet_multiple_correlation.csv"), row.names = FALSE)
    utils::write.csv(m$wmcc, file.path(dir, "wavelet_multiple_cross_correlation.csv"), row.names = FALSE)
    if (!is.null(m$wlmc)) utils::write.csv(m$wlmc, file.path(dir, "wavelet_local_multiple_correlation.csv"), row.names = FALSE)
    save_gg("pairwise_wavelet_correlation.png",
            if (length(m$vars) > 7) plot_wc_matrix(m$wc, m$vars, m$J, u) else plot_pairwise_wc(m$wc, m$J, u))
    save_gg("pairwise_wavelet_correlation_matrix.png", plot_wc_matrix(m$wc, m$vars, m$J, u), w = 12, h = 8)
    save_gg("wavelet_cross_correlation.png", plot_lag_heatmap(m$wcc, "wcc", m$J, u, sprintf("Wavelet cross-correlation: %s vs %s", m$wcc_pair[1], m$wcc_pair[2]), lag_unit = u))
    save_gg("wavelet_multiple_correlation.png", plot_wmc(m$wmc, u))
    save_gg("wavelet_multiple_cross_correlation.png", plot_lag_heatmap(m$wmcc, "wmcc", m$J, u, "Wavelet multiple cross-correlation", side_labels = attr(m$wmcc, "max_var"), lag_unit = u))
    if (!is.null(m$wlmc)) save_gg("wavelet_local_multiple_correlation.png", plot_wlmc(m$wlmc, m$J, u))
  }
  writeLines(c("Wavelet Approach - results export",
               paste0("Application by ", APP_AUTHOR$name, " (", APP_AUTHOR$affiliation, "), ", paste(APP_AUTHOR$emails, collapse = " / ")),
               paste("Generated:", format(Sys.time())),
               paste("Mode:", if (r$mode == "biv") "Bivariate (WaveletComp)" else "Multiple (MODWT)"),
               paste("Transformation:", r$transform), paste("Observations:", nrow(r$data))),
             file.path(dir, "README.txt"))
}

# ---------------------------------------------------------------------------
# About page
# ---------------------------------------------------------------------------
about_ui <- function() {
  tagList(
    div(class = "alert alert-light border mb-3",
        h4(class = "mb-1", "Author"),
        tags$p(class = "mb-1", HTML(sprintf("<strong>%s</strong>, %s", APP_AUTHOR$name, APP_AUTHOR$affiliation))),
        tags$p(class = "mb-1 small",
               icon("envelope"), " ",
               HTML(paste(sprintf('<a href="mailto:%s">%s</a>', APP_AUTHOR$emails, APP_AUTHOR$emails), collapse = " | "))),
        tags$p(class = "mb-0 small",
               icon("link"), " ",
               HTML(paste(sprintf('<a href="%s" target="_blank">%s</a>', APP_AUTHOR$links, names(APP_AUTHOR$links)), collapse = " | "))),
        tags$p(class = "mb-0 small text-muted mt-2",
               HTML(sprintf("Suggested citation: %s (%s). <em>Wavelet Approach: A Shiny application for bivariate and multiple wavelet analysis</em>. R Shiny application.",
                            APP_AUTHOR$name, format(Sys.Date(), "%Y"))))),
    h4("What this app does"),
    p("This app brings together two complementary wavelet toolkits for studying co-movement between financial and economic time series across time and frequency."),
    layout_columns(
      col_widths = c(6, 6),
      card(card_header("Bivariate wavelets (continuous Morlet transform)"),
           card_body(
             tags$ul(
               tags$li(tags$b("Wavelet power spectrum"), " of each series: where in time and at which periodicity the variance is concentrated."),
               tags$li(tags$b("Cross-wavelet power"), ": regions of common high power of the pair."),
               tags$li(tags$b("Wavelet coherence"), ": a localised correlation (0 to 1) in the time-frequency plane. Arrows give the phase difference: pointing right = in phase, left = anti-phase, up = x leads y by a quarter cycle, down = y leads x."),
               tags$li(tags$b("Average power & coherence"), " over time, by period, with Monte-Carlo significance dots."),
               tags$li(tags$b("Phase difference"), " averaged over a chosen period band, and the phase-difference image."),
               tags$li(tags$b("Reconstruction"), " of each series from its significant wavelet components.")
             ),
             p(class = "small text-muted", "Significance is assessed by simulation against a chosen null (white noise, shuffling, Fourier randomisation, AR(1) or ARIMA). The cone of influence marks edge-affected regions.")
           )),
      card(card_header("Multiple wavelets (MODWT)"),
           card_body(
             tags$ul(
               tags$li(tags$b("MODWT decomposition"), ": each series is split into detail series d1..dJ (period bands 2-4, 4-8, 8-16, ... observations) and a smooth sJ using the maximal overlap discrete wavelet transform. Boundary-affected coefficients are removed (brick wall)."),
               tags$li(tags$b("Wavelet correlation"), ": pairwise correlation at every scale with 95% confidence intervals."),
               tags$li(tags$b("Wavelet cross-correlation"), ": lead-lag structure of a chosen pair, scale by scale."),
               tags$li(tags$b("Wavelet multiple correlation (WMC)"), ": one coefficient per scale summarising the whole system; the label shows which variable is best explained by a linear combination of the others (Fernandez-Macho, 2012)."),
               tags$li(tags$b("Wavelet multiple cross-correlation (WMCC)"), ": the same, at leads and lags."),
               tags$li(tags$b("Wavelet local multiple correlation (WLMC)"), ": a moving-window, time-varying version of the WMC (Fernandez-Macho, 2018).")
             )
           ))
    ),
    h4(class = "mt-3", "How to use"),
    tags$ol(
      tags$li("Upload your own file (csv, txt or xlsx: one column per series, ideally with a date column) or pick an example dataset. ",
              "Choose which column holds the dates and, if auto-detection fails, the date format."),
      tags$li("Choose a transformation (e.g. log returns for prices) and, if needed, a volatility proxy."),
      tags$li("Select the approach and its settings, then click ", tags$b("Start analysis!"), "."),
      tags$li("Browse the tabs; use ", tags$b("Download results"), " for a zip with CSV tables and PNG figures.")
    ),
    h4(class = "mt-3", "Applications by Peterson Owusu Junior"),
    p(class = "small text-muted", "Studies in which the bivariate and multiple wavelet tools implemented here were applied to African and emerging-market data."),
    tags$ul(class = "small",
      lapply(OWUSU_REFERENCES, function(r) tags$li(r$text, if (!is.null(r$doi)) tags$a(href = paste0("https://doi.org/", r$doi), target = "_blank", paste0(" doi:", r$doi))))
    ),
    h4(class = "mt-3", "Software"),
    p(class = "small text-muted", "The app wraps these R packages. Peterson Owusu Junior is acknowledged in the W2CWM2C package documentation for reporting bugs that were fixed in the package."),
    tags$ul(class = "small",
      tags$li("Roesch, A., & Schmidbauer, H. (2025). WaveletComp: Computational wavelet analysis. R package version 1.2. ",
              tags$a(href = "https://doi.org/10.32614/CRAN.package.WaveletComp", target = "_blank", "doi:10.32614/CRAN.package.WaveletComp")),
      tags$li("Polanco-Martinez, J. M. (2021). W2CWM2C: A graphical tool for wavelet (cross) correlation and wavelet multiple (cross) correlation analysis. R package version 2.2. ",
              tags$a(href = "https://doi.org/10.32614/CRAN.package.W2CWM2C", target = "_blank", "doi:10.32614/CRAN.package.W2CWM2C"),
              " (acknowledgement: bug reports by Peterson Owusu Junior, among others)."),
      tags$li("Polanco-Martinez, J., & Fernandez-Macho, J. (2014). The package 'W2CWM2C': Description, features and applications. Computing in Science & Engineering, 16(6), 68-78. ",
              tags$a(href = "https://doi.org/10.1109/MCSE.2014.96", target = "_blank", "doi:10.1109/MCSE.2014.96")),
      tags$li("Fernandez-Macho, J. (2021). wavemulcor: Wavelet routines for global and local multiple regression and correlation. R package version 3.1.2. ",
              tags$a(href = "https://doi.org/10.32614/CRAN.package.wavemulcor", target = "_blank", "doi:10.32614/CRAN.package.wavemulcor")),
      tags$li("Whitcher, B. (2024). waveslim: Basic wavelet routines for one-, two-, and three-dimensional signal processing. R package version 1.8.5. ",
              tags$a(href = "https://doi.org/10.32614/CRAN.package.waveslim", target = "_blank", "doi:10.32614/CRAN.package.waveslim"))
    ),
    h4(class = "mt-3", "Methodological references"),
    tags$ul(class = "small",
      tags$li("Fernandez-Macho, J. (2012). Wavelet multiple correlation and cross-correlation: A multiscale analysis of Eurozone stock markets. Physica A, 391(4), 1097-1104. ",
              tags$a(href = "https://doi.org/10.1016/j.physa.2011.11.002", target = "_blank", "doi:10.1016/j.physa.2011.11.002")),
      tags$li("Fernandez-Macho, J. (2018). Time-localized wavelet multiple regression and correlation. Physica A, 492, 1226-1238."),
      tags$li("Gencay, R., Selcuk, F., & Whitcher, B. (2001). An introduction to wavelets and other filtering methods in finance and economics. Academic Press."),
      tags$li("Whitcher, B., Guttorp, P., & Percival, D. B. (2000). Wavelet analysis of covariance with application to atmospheric time series. Journal of Geophysical Research: Atmospheres, 105(D11), 14941-14962."),
      tags$li("Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet analysis. Bulletin of the American Meteorological Society, 79, 61-78."),
      tags$li("Ramsey, J. B. (2002). Wavelets in economics and finance: Past and future. Studies in Nonlinear Dynamics & Econometrics, 6(3).")
    )
  )
}

# Wavelet / time-frequency co-movement studies authored or co-authored by Peterson Owusu Junior
OWUSU_REFERENCES <- list(
  list(text = "Owusu Junior, P., Adam, A. M., & Tweneboah, G. (2017). Co-movement of real exchange rates in the West African Monetary Zone. Cogent Economics & Finance, 5(1), 1351807.",
       doi = "10.1080/23322039.2017.1351807"),
  list(text = "Owusu Junior, P., Boafo, B. K., Awuye, B. K., Bonsu, K., & Obeng-Tawiah, H. (2018). Co-movement of stock exchange indices and exchange rates in Ghana: A wavelet coherence analysis. Cogent Business & Management, 5(1), 1481559.",
       doi = "10.1080/23311975.2018.1481559"),
  list(text = "Owusu Junior, P., Tweneboah, G., & Adam, A. M. (2019). Interdependence of major exchange rates in Ghana: A wavelet coherence analysis. Journal of African Business, 20(3), 407-430.",
       doi = "10.1080/15228916.2019.1583973"),
  list(text = "Tweneboah, G., Owusu Junior, P., & Oseifuah, E. K. (2019). Integration of major African stock markets: Evidence from multi-scale wavelets correlation. Academy of Accounting and Financial Studies Journal, 23(6).",
       doi = NULL),
  list(text = "Owusu Junior, P., & Tweneboah, G. (2020). Are there asymmetric linkages between African stocks and exchange rates? Research in International Business and Finance, 54, 101245.",
       doi = "10.1016/j.ribaf.2020.101245"),
  list(text = "Owusu Junior, P., Alagidede, I., & Tweneboah, G. (2020). Shape-shift contagion in emerging markets equities: Evidence from frequency- and time-domain analysis. Economics and Business Letters, 9(3), 146-156.",
       doi = "10.17811/ebl.9.3.2020.146-156"),
  list(text = "Owusu Junior, P., Adam, A. M., Asafo-Adjei, E., Boateng, E., Hamidu, Z., & Awotwe, E. (2021). Time-frequency domain analysis of investor fear and expectations in stock markets of BRIC economies. Heliyon, 7(10), e08211.",
       doi = "10.1016/j.heliyon.2021.e08211"),
  list(text = "Asafo-Adjei, E., Boateng, E., Isshaq, Z., Idun, A. A.-A., Owusu Junior, P., & Adam, A. M. (2021). Financial sector and economic growth amid external uncertainty shocks: Insights into emerging economies. PLOS ONE, 16(11), e0259303.",
       doi = "10.1371/journal.pone.0259303"),
  list(text = "Agyei, S. K., Adam, A. M., Bossman, A., Asiamah, O., Owusu Junior, P., Asafo-Adjei, R., & Asafo-Adjei, E. (2022). Does volatility in cryptocurrencies drive the interconnectedness between the cryptocurrencies market? Insights from wavelets. Cogent Economics & Finance, 10(1), 2061682.",
       doi = "10.1080/23322039.2022.2061682"),
  list(text = "Agyei, S. K., Owusu Junior, P., Bossman, A., Asafo-Adjei, E., Asiamah, O., & Adam, A. M. (2022). Spillovers and contagion between BRIC and G7 markets: New evidence from time-frequency analysis. PLOS ONE, 17(7), e0271088.",
       doi = "10.1371/journal.pone.0271088"),
  list(text = "Barson, Z., Owusu Junior, P., & Adam, A. M. (2023). Comovement between commodity returns in Ghana: The role of exchange rates. Journal of Economic Structures, 12, 17.",
       doi = "10.1186/s40008-023-00312-z"),
  list(text = "Owusu Junior, P., Tetteh, J. E., Nkrumah-Boadu, B., & Adjei, A. N. K. (2024). Comovement of African stock markets: Any influence from the COVID-19 pandemic? Heliyon, 10(9), e29409.",
       doi = "10.1016/j.heliyon.2024.e29409"),
  list(text = "Yeboah, S. D., Owusu Junior, P., Idun, A. A.-A., Akorsu, P. K., & Fumey, M. P. (2025). Dynamic interdependence of major currencies and the US dollar: A wavelet coherence approach. Future Business Journal, 11(1).",
       doi = "10.1186/s43093-025-00636-1"),
  list(text = "Tetteh, J. E., & Owusu Junior, P. (2025). Co-movement between stock markets in advanced economies and Africa in times of uncertainty: A time-frequency domain approach. PLOS ONE, 20(11), e0334325.",
       doi = "10.1371/journal.pone.0334325")
)

shinyApp(ui, server)
