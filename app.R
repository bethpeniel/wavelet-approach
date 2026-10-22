# ===========================================================================
# Bi-Partial-MultipleWavelets_Explore - Bivariate, Partial & Multiple Wavelet Analysis
# Author : Peterson Owusu Junior (powusujunior@gmail.com)
#
# Engines:
#   * Bivariate  : continuous Morlet wavelet transform, cross-wavelet power,
#                  wavelet coherence and phase difference  (WaveletComp)
#   * Partial    : partial wavelet coherence with one control (biwavelet::pwtc,
#                  Mihanovic et al. 2009 / Ng & Chan 2012)
#   * Higher-order: multiple wavelet coherence and partial wavelet coherence
#                  with several controls (R/higher_order_coherence.R, after
#                  Aguiar-Conraria & Soares 2014, Ng & Chan 2012, Hu & Si 2016, 2021)
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
  library(biwavelet)
  library(shinycssloaders)
})

if (!exists("prepare_dataset")) {
  for (f in list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)) source(f, local = FALSE)
}

options(shiny.maxRequestSize = 30 * 1024^2, spinner.type = 6, spinner.color = "#2a78d6")

APP_AUTHOR <- list(
  name        = "Peterson Owusu Junior",
  affiliation = "University of Cape Coast, Ghana",
  email       = "peterson.owusu@ucc.edu.gh",          # institutional (used for the footer link)
  emails      = c("peterson.owusu@ucc.edu.gh", "powusujunior@gmail.com"),
  website     = "https://positivecarryacademy.com/",
  links       = c("Website"        = "https://positivecarryacademy.com/",
                  "Google Scholar" = "https://scholar.google.com/citations?user=r6ZueU4AAAAJ&hl=en",
                  "Scopus"         = "https://www.scopus.com/authid/detail.uri?authorId=58833238800",
                  "LinkedIn"       = "https://www.linkedin.com/in/peterson-owusu-junior-phd-80858334/",
                  "ResearchGate"   = "https://www.researchgate.net/profile/Peterson-Owusu-Junior",
                  "IDEAS/RePEc"    = "https://ideas.repec.org/e/pow25.html")
)
APP_AUTHOR$links <- APP_AUTHOR$links[nzchar(APP_AUTHOR$links)]   # drop links not yet filled in

# First-author wavelet papers cited in the footer (short form; full entries in OWUSU_REFERENCES)
FOOTER_PAPERS <- list(
  list(short = "Owusu Junior et al. (2018, Cogent Business & Management)", doi = "10.1080/23311975.2018.1481559"),
  list(short = "Owusu Junior, Tweneboah & Adam (2019, Journal of African Business)", doi = "10.1080/15228916.2019.1583973"),
  list(short = "Owusu Junior & Tweneboah (2020, Research in International Business and Finance)", doi = "10.1016/j.ribaf.2020.101245"),
  list(short = "Owusu Junior et al. (2021, Heliyon)", doi = "10.1016/j.heliyon.2021.e08211"),
  list(short = "Owusu Junior et al. (2024, Heliyon)", doi = "10.1016/j.heliyon.2024.e29409")
)

UNIT_CHOICES <- c("Daily" = "day", "Weekly" = "week", "Monthly" = "month",
                  "Quarterly" = "quarter", "Yearly" = "year", "Unspecified" = "observation")

BIV_TABS  <- c("tab_wt", "tab_xwt", "tab_coh", "tab_avg", "tab_phase", "tab_phimg", "tab_rec", "tab_bivtab")
MULT_TABS <- c("tab_modwt", "tab_wc", "tab_wcc", "tab_wmc", "tab_wmcc", "tab_wlmc", "tab_multtab")
PART_TABS <- c("tab_pw_wtc", "tab_pw_1", "tab_pw_2", "tab_pw_cmp", "tab_pw_tab")
HOC_TABS  <- c("tab_ho_mwc", "tab_ho_pwc", "tab_ho_cmp", "tab_ho_tab")

# One chart per card: a header that names it, optional controls, the plot, and
# an optional reading note. full_screen gives every chart an expand button.
plot_card <- function(title, id, height = "560px", controls = NULL, note = NULL) {
  card(full_screen = TRUE, class = "chart-card",
       card_header(title),
       card_body(controls, withSpinner(plotOutput(id, height = height)),
                 if (!is.null(note)) tags$p(class = "chart-note", note)))
}
table_card <- function(title, id) {
  card(class = "chart-card", card_header(title), card_body(withSpinner(DTOutput(id))))
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- page_sidebar(
  title = tags$span(tags$strong("Bi-Partial-MultipleWavelets_Explore"),
                    tags$small(class = "text-muted ms-2", "Bivariate, partial & multiple wavelet analysis")),
  # Neutral chrome, one accent. Colour is reserved for the data: a single blue
  # for actions and single-series marks, status colours only for status, and the
  # charts carry the sequential / diverging ramps defined in R/data_utils.R.
  theme = bs_theme(version = 5,
                   primary = "#2a78d6", secondary = "#52514e", success = "#0ca30c",
                   info = "#2a78d6", warning = "#fab219", danger = "#d03b3b",
                   bg = "#f9f9f7", fg = "#0b0b0b",
                   base_font = 'system-ui, -apple-system, "Segoe UI", Roboto, sans-serif',
                   "font-size-base" = "0.95rem", "border-radius" = "0.5rem"),
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
                       "Partial wavelet coherence, one control (biwavelet)" = "part",
                       "Multiple & partial wavelet coherence, several controls" = "hoc",
                       "Multiple wavelets (MODWT)" = "mult"), selected = "biv"),
        conditionalPanel(
          "input.mode == 'hoc'",
          selectInput("ho_y", "Dependent variable y", choices = NULL),
          selectizeInput("ho_x", "Predictors / controls X (2 or more)", choices = NULL, multiple = TRUE,
                         options = list(plugins = list("remove_button"))),
          helpText("Multiple wavelet coherence MWC(y | X) measures how much of y's wavelet variance all of X explain together. ",
                   "Partial wavelet coherence PWC(y, x | X\\x) is computed for every x in X, controlling for all the others."),
          radioButtons("ho_pwc_method", "Partial coherence formula", HO_PWC_METHODS, selected = "complex"),
          layout_columns(col_widths = c(6, 6),
                         selectInput("ho_dj", "Scale resolution (dj)", HO_DJ_CHOICES, selected = "1/12"),
                         numericInput("ho_nrands", "Monte-Carlo draws (0 = none)", value = 30, min = 0, max = 2000, step = 10)),
          layout_columns(col_widths = c(6, 6),
                         numericInput("ho_siglevel", "Significance level", value = 0.95, min = 0.8, max = 0.999, step = 0.01),
                         numericInput("ho_maxscale", "Max. period (blank = auto)", value = NA, min = 4, step = 1)),
          checkboxInput("ho_wtc", "Also compute plain wavelet coherence of each pair for comparison", value = TRUE),
          helpText("Each Monte-Carlo draw re-estimates every coherence map from AR(1) surrogates of all series (Hu & Si, 2016, 2021); ",
                   "30 draws on 1,500 daily observations with three predictors take about 2-4 minutes.")
        ),
        conditionalPanel(
          "input.mode == 'part'",
          selectInput("pw_y", "Dependent variable y", choices = NULL),
          layout_columns(col_widths = c(6, 6),
                         selectInput("pw_x1", "Variable of interest x1", choices = NULL),
                         selectInput("pw_x2", "Control variable x2", choices = NULL)),
          helpText("Computes PWTC(y, x1 | x2) and PWTC(y, x2 | x1): the coherence of y with one variable after removing the influence of the other."),
          layout_columns(col_widths = c(6, 6),
                         selectInput("pw_dj", "Scale resolution (dj)", PW_DJ_CHOICES, selected = "1/12"),
                         numericInput("pw_nrands", "Monte-Carlo draws (0 = none)", value = 50, min = 0, max = 2000, step = 10)),
          layout_columns(col_widths = c(6, 6),
                         numericInput("pw_siglevel", "Significance level", value = 0.95, min = 0.8, max = 0.999, step = 0.01),
                         numericInput("pw_maxscale", "Max. period (blank = auto)", value = NA, min = 4, step = 1)),
          checkboxInput("pw_wtc", "Also compute plain wavelet coherence for comparison", value = TRUE),
          helpText("Significance testing is Monte-Carlo based; 50-100 draws take roughly 1-3 minutes on 1,500 daily observations.")
        ),
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
      ),
      accordion_panel(
        "Display", icon = icon("palette"),
        selectInput("palette", "Colour scale for power / coherence maps", PALETTE_CHOICES, selected = "blue"),
        helpText("A single hue from light to dark reads magnitude correctly and survives colour-vision deficiency and greyscale printing. ",
                 "The rainbow option reproduces the classic look of the wavelet literature.")
      )
    ),
    # primary action stays visible however long the settings list gets
    div(class = "sidebar-actions",
        actionButton("run", "Start analysis", class = "btn-primary w-100", icon = icon("play")),
        downloadButton("download", "Download results", class = "btn-outline-secondary w-100 mt-2"),
        uiOutput("status"))
  ),

  uiOutput("data_info"),
  uiOutput("run_summary"),

  # Results: every tab visible in one strip (data tabs always; the tabs of the
  # two approaches that are not selected are hidden by the server).
  navset_underline(
    id = "tabs",
      nav_panel("Data table", value = "tab_data", table_card("Loaded data", "data_table")),
      nav_panel("Raw data", value = "tab_raw", plot_card("Raw series", "raw_plot", "620px")),
      nav_panel("Transformed data", value = "tab_trans", plot_card("Series as analysed", "trans_plot", "620px")),
      nav_panel("Summary statistics", value = "tab_stats", table_card("Summary statistics of the analysed series", "stats_table")),

    # ---- bivariate --------------------------------------------------------
      nav_panel("Wavelet power", value = "tab_wt",
                plot_card("Wavelet power spectrum", "wt_plot", controls = uiOutput("wt_series_ui"))),
      nav_panel("Cross-wavelet power", value = "tab_xwt", plot_card("Cross-wavelet power", "xwt_plot")),
      nav_panel("Wavelet coherence", value = "tab_coh",
                plot_card("Wavelet coherence", "coh_plot",
                          note = "Arrows give the phase difference: right = in phase, left = anti-phase, up = x leads y by a quarter cycle, down = y leads x. The shaded cone marks edge effects.")),
      nav_panel("Average power & coherence", value = "tab_avg",
                layout_columns(col_widths = c(6, 6),
                               plot_card("Average cross-wavelet power", "avg_wp_plot", "480px"),
                               plot_card("Average coherence", "avg_wc_plot", "480px"))),
      nav_panel("Phase difference", value = "tab_phase",
                plot_card("Phase difference in a period band", "phase_plot", "520px", controls = uiOutput("phase_controls"),
                          note = "Phase difference = phase(x) - phase(y). (0, π/2): in phase, x leads; (-π/2, 0): in phase, y leads; (π/2, π): anti-phase, y leads; (-π, -π/2): anti-phase, x leads.")),
      nav_panel("Phase-difference image", value = "tab_phimg", plot_card("Phase difference across time and period", "phimg_plot")),
      nav_panel("Reconstruction", value = "tab_rec",
                plot_card("Reconstruction from significant components", "rec_plot", "520px", controls = uiOutput("rec_controls"))),
      nav_panel("Tables", value = "tab_bivtab",
                table_card("Average cross-wavelet power and coherence by period", "coh_table"),
                table_card("Phase summary inside significant coherence regions", "phase_table")),

    # ---- partial ----------------------------------------------------------
      nav_panel("Wavelet coherence (pairs)", value = "tab_pw_wtc",
                plot_card("Wavelet coherence of each pair", "pw_wtc_plot", controls = uiOutput("pw_wtc_controls"))),
      nav_panel("Partial coherence y~x1 | x2", value = "tab_pw_1", plot_card("Partial wavelet coherence, y ~ x1 | x2", "pw_plot1")),
      nav_panel("Partial coherence y~x2 | x1", value = "tab_pw_2", plot_card("Partial wavelet coherence, y ~ x2 | x1", "pw_plot2")),
      nav_panel("WTC vs PWTC", value = "tab_pw_cmp",
                plot_card("Coherence before and after conditioning", "pw_cmp_plot", "900px", controls = uiOutput("pw_cmp_controls"),
                          note = "Top: ordinary wavelet coherence of the pair. Bottom: the same pair after partialling out the control variable. Regions that stay dark in the bottom panel are co-movements not explained by the control; regions that fade were driven by it.")),
      nav_panel("Partial tables", value = "tab_pw_tab",
                table_card("Average coherence by period (inside the cone of influence)", "pw_table")),

    # ---- higher-order: multiple & partial coherence -------------------------
      nav_panel("Multiple coherence (MWC)", value = "tab_ho_mwc",
                plot_card("Multiple wavelet coherence of y on all predictors", "ho_mwc_plot",
                          note = "Dark regions: at that time and period the predictors jointly explain y's wavelet variance. Adding predictors can only raise the MWC, so compare predictor sets rather than reading the level alone (Hu & Si, 2016). The thick contour marks the Monte-Carlo significance level.")),
      nav_panel("Partial coherence (several controls)", value = "tab_ho_pwc",
                plot_card("Partial wavelet coherence of y with one predictor, controlling for the others", "ho_pwc_plot",
                          controls = uiOutput("ho_pwc_controls"),
                          note = "With the complex formula the arrows show the partial phase difference (right = in phase, left = anti-phase). The classical formula has no phase information.")),
      nav_panel("WTC vs PWC vs MWC", value = "tab_ho_cmp",
                plot_card("Bivariate, partial and multiple coherence for one predictor", "ho_cmp_plot", "1200px",
                          controls = uiOutput("ho_cmp_controls"),
                          note = "Top: plain coherence of the pair. Middle: the same pair after partialling out the other predictors. Bottom: the multiple coherence of y on all predictors. Co-movement that survives in the middle panel is specific to this predictor; what fades was carried by the controls.")),
      nav_panel("Higher-order tables", value = "tab_ho_tab",
                table_card("Average coherence by period (inside the cone of influence)", "ho_table"),
                table_card("Mean partial phase inside significant regions (complex formula only)", "ho_phase_table")),

    # ---- multiple ---------------------------------------------------------
      nav_panel("MODWT decomposition", value = "tab_modwt",
                plot_card("MODWT decomposition", "modwt_plot", "760px", controls = uiOutput("modwt_controls"))),
      nav_panel("Wavelet correlation", value = "tab_wc",
                plot_card("Pairwise wavelet correlation", "wc_plot", "640px", controls = uiOutput("wc_controls"))),
      nav_panel("Wavelet cross-correlation", value = "tab_wcc",
                plot_card("Wavelet cross-correlation by lag and level", "wcc_plot", "520px"),
                plot_card("Wavelet cross-correlation by level, with confidence bands", "wcc_lines", "520px")),
      nav_panel("Multiple correlation", value = "tab_wmc", plot_card("Wavelet multiple correlation", "wmc_plot", "600px")),
      nav_panel("Multiple cross-correlation", value = "tab_wmcc",
                plot_card("Wavelet multiple cross-correlation by lag and level", "wmcc_plot", "520px"),
                plot_card("Wavelet multiple cross-correlation by level, with confidence bands", "wmcc_lines", "520px")),
      nav_panel("Local multiple correlation", value = "tab_wlmc",
                plot_card("Wavelet local multiple correlation over time", "wlmc_plot", "520px"),
                plot_card("Wavelet local multiple correlation by level", "wlmc_lines", "620px")),
      nav_panel("Tables", value = "tab_multtab",
                table_card("Wavelet multiple correlation", "wmc_table"),
                table_card("Pairwise wavelet correlation", "wc_table"),
                table_card("Wavelet multiple cross-correlation", "wmcc_table")),

    nav_spacer(),
    nav_panel("About", value = "tab_about", icon = icon("circle-info"), uiOutput("about"))
  ),

  tags$footer(class = "app-footer",
              HTML(sprintf('<strong>Bi-Partial-MultipleWavelets_Explore</strong> by <a href="mailto:%s">%s</a>, %s (<a href="%s" target="_blank">%s</a>). ',
                           APP_AUTHOR$email, APP_AUTHOR$name, APP_AUTHOR$affiliation,
                           APP_AUTHOR$website, sub("^https?://", "", sub("/$", "", APP_AUTHOR$website)))),
              HTML(paste0("A companion to the author's published wavelet research, including ",
                          paste(vapply(FOOTER_PAPERS, function(p) sprintf('<a href="https://doi.org/%s" target="_blank">%s</a>', p$doi, p$short), ""),
                                collapse = "; "),
                          "; the full list is on the About tab. Built with R Shiny, WaveletComp, biwavelet, waveslim and wavemulcor.")))
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
    unit <- switch(input$rep_dataset, resilience = "day", africa = "week", crypto = "day", reit = "day", ghana = "month", "observation")
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
    updateSelectInput(session, "pw_y",  choices = vars, selected = vars[1])
    updateSelectInput(session, "pw_x1", choices = vars, selected = vars[min(2, length(vars))])
    updateSelectInput(session, "pw_x2", choices = vars, selected = vars[min(3, length(vars))])
    updateSelectInput(session, "ho_y", choices = vars, selected = vars[1])
    updateSelectizeInput(session, "ho_x", choices = vars, selected = vars[seq.int(2, min(3, length(vars)))])
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
    for (t in PART_TABS) if (input$mode == "part") nav_show("tabs", t) else nav_hide("tabs", t)
    for (t in HOC_TABS)  if (input$mode == "hoc")  nav_show("tabs", t) else nav_hide("tabs", t)
    for (t in MULT_TABS) if (input$mode == "mult") nav_show("tabs", t) else nav_hide("tabs", t)
  }, ignoreNULL = FALSE)

  # After a run, jump to the first result tab of the chosen approach
  observeEvent(res(), {
    first <- switch(res()$mode, biv = "tab_coh", part = "tab_pw_1", hoc = "tab_ho_mwc", mult = "tab_wc")
    nav_select("tabs", first)
  })

  output$data_info <- renderUI({
    if (input$source == "own" && is.null(input$file))
      return(div(class = "empty-state",
                 h5("Get started"),
                 tags$ol(class = "mb-0",
                   tags$li(tags$b("Load data"), ": upload a csv, txt or xlsx file in the sidebar, or switch to an example dataset."),
                   tags$li(tags$b("Choose the approach"), ": bivariate coherence, partial coherence with one control, multiple and partial coherence with several controls, or multiple wavelets (MODWT)."),
                   tags$li(tags$b("Start analysis"), ": each result has its own tab above; Download results bundles tables and figures."))))
    r <- raw(); req(r)
    rng <- if (is.null(r$date)) "no dates" else paste(format(range(r$date)), collapse = " to ")
    div(class = "dataset-strip",
        tags$span(class = "ds-main", sprintf("%d observations, %d variables", nrow(r$data), ncol(r$data))),
        tags$span(class = "ds-sep"), tags$span(rng),
        tags$span(class = "ds-sep"), tags$span(class = "ds-vars", paste(names(r$data), collapse = ", ")),
        if (length(r$note)) tags$div(class = "ds-note", icon("circle-info"), " ", paste(r$note, collapse = " ")))
  })

  # Stat tiles summarising the last run (the numbers people look for first)
  output$run_summary <- renderUI({
    r <- res()
    mode_lab <- switch(r$mode, biv = "Bivariate coherence", part = "Partial coherence", hoc = "Multiple & partial coherence", mult = "Multiple wavelets")
    detail <- switch(r$mode,
                     biv  = paste(r$biv$pair, collapse = " vs "),
                     part = sprintf("%s | %s, %s", r$part$y, r$part$x1, r$part$x2),
                     hoc  = sprintf("%s | %s", r$hoc$y, paste(r$hoc$X, collapse = ", ")),
                     mult = sprintf("%d series, J = %d", length(r$mult$vars), r$mult$J))
    span_lab <- if (is.null(r$date)) "index 1..N" else paste(format(range(r$date), "%b %Y"), collapse = " to ")
    tile <- function(label, value, sub = NULL)
      div(class = "stat-tile", div(class = "stat-label", label), div(class = "stat-value", value),
          if (!is.null(sub)) div(class = "stat-sub", sub))
    div(class = "stat-row",
        tile("Approach", mode_lab, detail),
        tile("Observations", format(nrow(r$data), big.mark = ","), r$transform),
        tile("Sample", span_lab, paste0(ncol(r$data), " variables, ", unit_plural(r$unit))),
        tile("Run time", paste0(r$elapsed, " s"), if (length(r$transform_note)) "rows dropped, see sidebar" else "complete"))
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
                  palette = input$palette %||% "blue",
                  transform = transform_label(input$transformation, input$lag_n, input$volatility),
                  transform_note = prep$note,
                  stats = summary_stats(prep$data))
      if (length(prep$note)) showNotification(prep$note, type = "warning", duration = 12)
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
      } else if (input$mode == "part") {
        validate(need(ncol(prep$data) >= 3, "Partial wavelet coherence needs at least three variables."),
                 need(length(unique(c(input$pw_y, input$pw_x1, input$pw_x2))) == 3, "y, x1 and x2 must be three different variables."))
        nr <- max(0, as.integer(input$pw_nrands))
        incProgress(0.1, detail = sprintf("Partial wavelet coherence (%d Monte-Carlo draws)", nr))
        out$part <- run_partial(prep$date, prep$data, input$pw_y, input$pw_x1, input$pw_x2,
                                dj = eval(parse(text = input$pw_dj)), nrands = nr,
                                sig.level = input$pw_siglevel, max.scale = input$pw_maxscale,
                                compute_wtc = isTRUE(input$pw_wtc))
      } else if (input$mode == "hoc") {
        X <- intersect(input$ho_x, names(prep$data))
        X <- setdiff(X, input$ho_y)
        validate(need(ncol(prep$data) >= 3, "Multiple / partial wavelet coherence needs at least three variables."),
                 need(length(X) >= 2, "Select at least two predictor / control variables different from y."),
                 need(input$ho_pwc_method != "classical" || length(X) == 2,
                      "The classical partial-coherence formula handles one control: select exactly two variables in X, or switch to the complex formula."))
        nr <- max(0, as.integer(input$ho_nrands))
        incProgress(0.05, detail = sprintf("Multiple & partial wavelet coherence (%d Monte-Carlo draws)", nr))
        out$hoc <- run_higher_order(prep$date, prep$data, input$ho_y, X,
                                    dj = eval(parse(text = input$ho_dj)), nrands = nr,
                                    sig.level = input$ho_siglevel, max.scale = input$ho_maxscale,
                                    pwc_method = input$ho_pwc_method, compute_wtc = isTRUE(input$ho_wtc),
                                    progress = function(f) setProgress(0.15 + 0.75 * f,
                                                                       detail = sprintf("Monte-Carlo significance: %d%%", round(100 * f))))
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
    tagList(
      div(class = "status-line status-ok", icon("circle-check"),
          sprintf(" Done in %s s. %s, N = %d.", r$elapsed, r$transform, nrow(r$data))),
      if (length(r$transform_note)) div(class = "status-line status-warn", icon("triangle-exclamation"), " ", r$transform_note))
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
    plot_wt_power(r$biv, s, r$unit, r$siglvl, r$palette)
  }, res = 100)
  output$xwt_plot <- renderPlot({ r <- biv(); plot_cross_power(r$biv, r$unit, r$siglvl, r$palette) }, res = 100)
  output$coh_plot <- renderPlot({ r <- biv(); plot_coherence(r$biv, r$unit, r$siglvl, r$palette) }, res = 100)
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
  output$phimg_plot <- renderPlot({ r <- biv(); plot_phasediff_image(r$biv, r$unit, r$siglvl, r$palette) }, res = 100)

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

  # ---- partial outputs ---------------------------------------------------
  part <- reactive({ r <- res(); validate(need(r$mode == "part", "Run a partial wavelet coherence analysis to see this tab.")); r })

  pw_pair_choices <- function(p) stats::setNames(c(1, 2), c(sprintf("%s ~ %s | %s", p$y, p$x1, p$x2),
                                                             sprintf("%s ~ %s | %s", p$y, p$x2, p$x1)))
  output$pw_wtc_controls <- renderUI({
    r <- part(); p <- r$part
    validate(need(!is.null(p$wtc_y_x1), "Plain wavelet coherence was not computed (tick the option in the sidebar)."))
    radioButtons("pw_wtc_pair", NULL, inline = TRUE,
                 choices = stats::setNames(c(1, 2), c(sprintf("%s and %s", p$y, p$x1), sprintf("%s and %s", p$y, p$x2))))
  })
  output$pw_wtc_plot <- renderPlot({
    r <- part(); validate(need(!is.null(r$part$wtc_y_x1), ""))
    plot_wtc_pair(r$part, as.integer(input$pw_wtc_pair %||% 1), r$unit, r$palette)
  }, res = 100)
  output$pw_plot1 <- renderPlot({ r <- part(); plot_partial_pair(r$part, 1, r$unit, r$palette) }, res = 100)
  output$pw_plot2 <- renderPlot({ r <- part(); plot_partial_pair(r$part, 2, r$unit, r$palette) }, res = 100)
  output$pw_cmp_controls <- renderUI({
    r <- part(); validate(need(!is.null(r$part$wtc_y_x1), "Plain wavelet coherence was not computed (tick the option in the sidebar)."))
    radioButtons("pw_cmp_pair", NULL, choices = pw_pair_choices(r$part), inline = TRUE)
  })
  output$pw_cmp_plot <- renderPlot({
    r <- part(); validate(need(!is.null(r$part$wtc_y_x1), ""))
    plot_partial_compare(r$part, as.integer(input$pw_cmp_pair %||% 1), r$unit, r$palette)
  }, res = 100)
  output$pw_table <- renderDT({
    r <- part()
    datatable(partial_table(r$part), rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "compact stripe")
  })

  # ---- higher-order outputs ----------------------------------------------
  hoc <- reactive({ r <- res(); validate(need(r$mode == "hoc", "Run a multiple & partial wavelet coherence analysis to see this tab.")); r })

  ho_choices <- function(h) stats::setNames(seq_along(h$X), vapply(seq_along(h$X), function(k)
    sprintf("%s ~ %s | %s", h$y, h$X[k], paste(attr(h$pwc[[k]], "controls"), collapse = ", ")), ""))
  output$ho_pwc_controls <- renderUI({
    r <- hoc(); radioButtons("ho_pwc_which", NULL, choices = ho_choices(r$hoc), inline = TRUE)
  })
  output$ho_cmp_controls <- renderUI({
    r <- hoc(); radioButtons("ho_cmp_which", NULL, choices = ho_choices(r$hoc), inline = TRUE)
  })
  output$ho_mwc_plot <- renderPlot({ r <- hoc(); plot_ho_mwc(r$hoc, r$unit, r$palette) }, res = 100)
  output$ho_pwc_plot <- renderPlot({
    r <- hoc(); k <- as.integer(input$ho_pwc_which %||% 1); if (is.na(k) || k > length(r$hoc$X)) k <- 1
    plot_ho_pwc(r$hoc, k, r$unit, r$palette)
  }, res = 100)
  output$ho_cmp_plot <- renderPlot({
    r <- hoc(); k <- as.integer(input$ho_cmp_which %||% 1); if (is.na(k) || k > length(r$hoc$X)) k <- 1
    plot_ho_compare(r$hoc, k, r$unit, r$palette)
  }, res = 100)
  output$ho_table <- renderDT({
    r <- hoc()
    datatable(ho_table(r$hoc), rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "compact stripe")
  })
  output$ho_phase_table <- renderDT({
    r <- hoc(); df <- ho_phase_table(r$hoc)
    validate(need(!is.null(df), "Partial phase is only defined for the complex formula."))
    datatable(df, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE), class = "compact stripe")
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
    save_png("wavelet_power_x.png", plot_wt_power(b, b$pair[1], r$unit, r$siglvl, r$palette))
    save_png("wavelet_power_y.png", plot_wt_power(b, b$pair[2], r$unit, r$siglvl, r$palette))
    save_png("cross_wavelet_power.png", plot_cross_power(b, r$unit, r$siglvl, r$palette))
    save_png("wavelet_coherence.png", plot_coherence(b, r$unit, r$siglvl, r$palette))
    save_png("average_coherence.png", plot_avg(b, "wc", r$unit))
    save_png("phase_difference_image.png", plot_phasediff_image(b, r$unit, r$siglvl, r$palette))
  } else if (r$mode == "part") {
    p <- r$part
    utils::write.csv(partial_table(p), file.path(dir, "partial_coherence_by_period.csv"), row.names = FALSE)
    utils::write.csv(partial_matrix(p$pwtc_y_x1, r$date), file.path(dir, sprintf("pwtc_%s_%s_given_%s_period_by_time.csv", p$y, p$x1, p$x2)))
    utils::write.csv(partial_matrix(p$pwtc_y_x2, r$date), file.path(dir, sprintf("pwtc_%s_%s_given_%s_period_by_time.csv", p$y, p$x2, p$x1)))
    save_png(sprintf("partial_coherence_%s_%s_given_%s.png", p$y, p$x1, p$x2), plot_partial_pair(p, 1, r$unit, r$palette))
    save_png(sprintf("partial_coherence_%s_%s_given_%s.png", p$y, p$x2, p$x1), plot_partial_pair(p, 2, r$unit, r$palette))
    if (!is.null(p$wtc_y_x1)) {
      save_png(sprintf("wavelet_coherence_%s_%s.png", p$y, p$x1), plot_wtc_pair(p, 1, r$unit, r$palette))
      save_png(sprintf("wavelet_coherence_%s_%s.png", p$y, p$x2), plot_wtc_pair(p, 2, r$unit, r$palette))
      save_png(sprintf("wtc_vs_pwtc_%s_%s_given_%s.png", p$y, p$x1, p$x2), plot_partial_compare(p, 1, r$unit, r$palette), h = 10)
      save_png(sprintf("wtc_vs_pwtc_%s_%s_given_%s.png", p$y, p$x2, p$x1), plot_partial_compare(p, 2, r$unit, r$palette), h = 10)
    }
  } else if (r$mode == "hoc") {
    h <- r$hoc
    utils::write.csv(ho_table(h), file.path(dir, "higher_order_coherence_by_period.csv"), row.names = FALSE)
    pt <- ho_phase_table(h); if (!is.null(pt)) utils::write.csv(pt, file.path(dir, "partial_phase_by_period.csv"), row.names = FALSE)
    utils::write.csv(ho_matrix(h$mwc, r$date), file.path(dir, sprintf("mwc_%s_given_%s_period_by_time.csv", h$y, paste(h$X, collapse = "_"))))
    save_png(sprintf("mwc_%s_given_%s.png", h$y, paste(h$X, collapse = "_")), plot_ho_mwc(h, r$unit, r$palette))
    for (k in seq_along(h$X)) {
      ctrl <- paste(attr(h$pwc[[k]], "controls"), collapse = "_")
      utils::write.csv(ho_matrix(h$pwc[[k]], r$date), file.path(dir, sprintf("pwc_%s_%s_given_%s_period_by_time.csv", h$y, h$X[k], ctrl)))
      save_png(sprintf("pwc_%s_%s_given_%s.png", h$y, h$X[k], ctrl), plot_ho_pwc(h, k, r$unit, r$palette))
      save_png(sprintf("wtc_pwc_mwc_%s_%s.png", h$y, h$X[k]), plot_ho_compare(h, k, r$unit, r$palette), h = if (is.null(h$wtc)) 10 else 14)
    }
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
  writeLines(c("Bi-Partial-MultipleWavelets_Explore - results export",
               paste0("Application by ", APP_AUTHOR$name, " (", APP_AUTHOR$affiliation, "), ", paste(APP_AUTHOR$emails, collapse = " / ")),
               paste("Generated:", format(Sys.time())),
               paste("Mode:", switch(r$mode, biv = "Bivariate (WaveletComp)", part = "Partial wavelet coherence, one control (biwavelet)",
                                     hoc = "Multiple & partial wavelet coherence, several controls", "Multiple (MODWT)")),
               if (r$mode == "hoc") c("",
                 "Method: multiple wavelet coherence and complex partial wavelet coherence computed from the smoothed",
                 "cross-wavelet spectral matrix (Aguiar-Conraria & Soares 2014, Appendix A; Hu & Si 2016, 2021; Ng & Chan 2012).",
                 "Please acknowledge in publications: 'Partial and multiple wavelet coherence software was provided by",
                 "E. K. W. Ng and T. W. Kwok (http://www.cityu.edu.hk/lccic/wavelet), by W. Hu and B. Si (Supplement of",
                 "Hu & Si 2016; https://doi.org/10.6084/m9.figshare.13031123), and by L. Aguiar-Conraria and M. J. Soares",
                 "(ASToolbox, https://sites.google.com/site/aguiarconraria/joanasoares-wavelets); the R implementation",
                 "builds on biwavelet (Gouhier, Grinsted & Simko) and the wavelet-coherence code of A. Grinsted.'"),
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
               HTML(sprintf("Suggested citation: %s (%s). <em>Bi-Partial-MultipleWavelets_Explore: A Shiny application for bivariate, partial and multiple wavelet analysis</em>. R Shiny application.",
                            APP_AUTHOR$name, format(Sys.Date(), "%Y"))))),
    h4("What this app does"),
    p("This app brings together complementary wavelet toolkits for studying co-movement between financial and economic time series across time and frequency: bivariate coherence, partial coherence with one or several controls, multiple wavelet coherence, and MODWT-based multiple correlation."),
    layout_columns(
      col_widths = c(6, 6),
      card(card_header("Multiple & partial wavelet coherence with several controls (higher-order coherencies)"),
           card_body(
             p("Ordinary wavelet coherence is bivariate. Following the Fourier-domain concepts of multiple and partial coherency (Koopmans, 1995; Priestley, 1981), ",
               "Aguiar-Conraria & Soares (2014, Section 4 and Appendix A) and Ng & Chan (2012) carried them to the time-frequency plane, and Hu & Si (2016, 2021) ",
               "gave the general matrix forms with Monte-Carlo significance. All quantities come from the Hermitian matrix S of smoothed (cross-)wavelet spectra at each time and period:"),
             tags$ul(
               tags$li(tags$b("Multiple wavelet coherence MWC(y | X)"), ": the share of y's wavelet variance explained jointly by all predictors X, ",
                       "R² = S", tags$sub("yX"), " S", tags$sub("XX"), tags$sup("-1"), " S", tags$sub("Xy"), " / S", tags$sub("yy"),
                       " (akin to the squared multiple correlation). With two predictors it reduces to Ng & Chan's (2012) eq. 7. Adding predictors never lowers it, so use it to compare predictor sets."),
               tags$li(tags$b("Partial wavelet coherence PWC(y, x | Z)"), ": the coherence of y with one predictor x after removing every other selected variable Z, ",
                       "computed from the partial spectral matrix S", tags$sub("AA·Z"), " = S", tags$sub("AA"), " - S", tags$sub("AZ"), " S", tags$sub("ZZ"), tags$sup("-1"), " S", tags$sub("ZA"),
                       " (akin to partial correlation). The ", tags$b("complex formula"), " (Aguiar-Conraria & Soares 2014, eq. A2; Hu & Si 2021, eq. 1) keeps the phase, so arrows show the partial lead-lag; ",
                       "the ", tags$b("classical formula"), " (Mihanovic et al. 2009 as coded by Ng & Chan 2012, and in biwavelet) uses real coherences, handles one control and has no phase. Hu & Si (2021) show it understates the partial coherence unless the two complex terms are collinear."),
               tags$li(tags$b("Significance"), ": AR(1) surrogates of every series are drawn, each map is recomputed, and the chosen quantile of the values outside the cone of influence at each period is the threshold (Grinsted et al. 2004; Hu & Si 2016, 2021: about 300 draws are enough; 1,000 recommended when time allows)."),
               tags$li(tags$b("Tabs"), ": the MWC map; a PWC map for each predictor with the others as controls; a stacked WTC / PWC / MWC comparison; and tables of average coherence, share of significant area and mean partial phase by period.")
             ),
             p(class = "small text-muted",
               "Source acknowledgement: partial and multiple wavelet coherence software was provided by E. K. W. Ng and T. W. Kwok (City University of Hong Kong, ",
               tags$a(href = "https://www.cityu.edu.hk/lccic/wavelet/", target = "_blank", "cityu.edu.hk/lccic/wavelet"), "), by W. Hu and B. Si (Supplement of Hu & Si 2016 and ",
               tags$a(href = "https://doi.org/10.6084/m9.figshare.13031123", target = "_blank", "figshare"), "), and by L. Aguiar-Conraria and M. J. Soares (ASToolbox, ",
               tags$a(href = "https://sites.google.com/site/aguiarconraria/joanasoares-wavelets", target = "_blank", "ASToolbox"), "). ",
               "This R implementation re-derives their formulas on top of biwavelet's wavelet transform and smoothing, which in turn follow Torrence & Compo (1998) and A. Grinsted's wavelet-coherence package.")
           )),
      card(card_header("Partial wavelet coherence, one control (biwavelet)"),
           card_body(
             p("Ordinary wavelet coherence between y and x1 can be inflated by a third variable that drives both. ",
               "Partial wavelet coherence, PWTC(y, x1 | x2), removes the influence of x2 and shows the co-movement that remains, ",
               "in the same time-frequency plane (Mihanovic et al., 2009; Ng & Chan, 2012). This mode calls biwavelet::pwtc, i.e. the classical real-coherence formula."),
             tags$ul(
               tags$li(tags$b("Wavelet coherence (pairs)"), ": the plain coherence of y with x1 and with x2, for reference."),
               tags$li(tags$b("Partial coherence y~x1 | x2"), " and ", tags$b("y~x2 | x1"), ": the two conditional coherences, each with the significance contour, cone of influence and phase arrows (the arrows are the plain cross-wavelet phase)."),
               tags$li(tags$b("WTC vs PWTC"), ": the pair's coherence before and after conditioning, stacked for direct comparison."),
               tags$li(tags$b("Tables"), ": average coherence and partial coherence by period, and the share of time each is significant.")
             ),
             p(class = "small text-muted", "Significance comes from a Monte-Carlo test against AR(1) surrogates (biwavelet); set the number of draws to 0 to skip it when exploring. For several controls, or for partial phase information, use the higher-order mode.")
           ))
    ),
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
      tags$li("Gouhier, T. C., Grinsted, A., & Simko, V. (2024). biwavelet: Conduct univariate and bivariate wavelet analyses. R package version 0.20.22. ",
              tags$a(href = "https://cran.r-project.org/package=biwavelet", target = "_blank", "CRAN"), " | ",
              tags$a(href = "https://github.com/tgouhier/biwavelet", target = "_blank", "GitHub")),
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
      tags$li("Torrence, C., & Webster, P. J. (1999). Interdecadal changes in the ENSO-monsoon system. Journal of Climate, 12(8), 2679-2690. (wavelet coherence)"),
      tags$li("Grinsted, A., Moore, J. C., & Jevrejeva, S. (2004). Application of the cross wavelet transform and wavelet coherence to geophysical time series. Nonlinear Processes in Geophysics, 11, 561-566. ",
              tags$a(href = "https://doi.org/10.5194/npg-11-561-2004", target = "_blank", "doi:10.5194/npg-11-561-2004")),
      tags$li("Mihanovic, H., Orlic, M., & Pasaric, Z. (2009). Diurnal thermocline oscillations driven by tidal flow around an island in the Middle Adriatic. Journal of Marine Systems, 78, S157-S168. (partial wavelet coherence) ",
              tags$a(href = "https://doi.org/10.1016/j.jmarsys.2009.01.021", target = "_blank", "doi:10.1016/j.jmarsys.2009.01.021")),
      tags$li("Ng, E. K. W., & Chan, J. C. L. (2012). Geophysical applications of partial wavelet coherence and multiple wavelet coherence. Journal of Atmospheric and Oceanic Technology, 29(12), 1845-1853. ",
              tags$a(href = "https://doi.org/10.1175/JTECH-D-12-00056.1", target = "_blank", "doi:10.1175/JTECH-D-12-00056.1")),
      tags$li("Aguiar-Conraria, L., & Soares, M. J. (2014). The continuous wavelet transform: Moving beyond uni- and bivariate analysis. Journal of Economic Surveys, 28(2), 344-375. (Section 4 'Higher order coherencies: partial and multiple coherencies' and Appendix A) ",
              tags$a(href = "https://doi.org/10.1111/joes.12012", target = "_blank", "doi:10.1111/joes.12012")),
      tags$li("Aguiar-Conraria, L., & Soares, M. J. (2018). ASToolbox: Wavelet tools for economists (MATLAB toolbox, version April 2018, with the function MPAWCOG for multiple and partial wavelet coherency). ",
              tags$a(href = "https://sites.google.com/site/aguiarconraria/joanasoares-wavelets", target = "_blank", "ASToolbox")),
      tags$li("Aguiar-Conraria, L., Azevedo, N., & Soares, M. J. (2008). Using wavelets to decompose the time-frequency effects of monetary policy. Physica A, 387(12), 2863-2878. ",
              tags$a(href = "https://doi.org/10.1016/j.physa.2008.01.063", target = "_blank", "doi:10.1016/j.physa.2008.01.063")),
      tags$li("Hu, W., & Si, B. C. (2016). Technical note: Multiple wavelet coherence for untangling scale-specific and localized multivariate relationships in geosciences. Hydrology and Earth System Sciences, 20, 3183-3191. ",
              tags$a(href = "https://doi.org/10.5194/hess-20-3183-2016", target = "_blank", "doi:10.5194/hess-20-3183-2016")),
      tags$li("Hu, W., & Si, B. (2021). Technical note: Improved partial wavelet coherency for understanding scale-specific and localized bivariate relationships in geosciences. Hydrology and Earth System Sciences, 25, 321-331. ",
              tags$a(href = "https://doi.org/10.5194/hess-25-321-2021", target = "_blank", "doi:10.5194/hess-25-321-2021"), " | code: ",
              tags$a(href = "https://doi.org/10.6084/m9.figshare.13031123", target = "_blank", "figshare")),
      tags$li("Ng, E. K. W., & Kwok, T. W. (2012). Wavelet package v1.0: Partial wavelet coherence, multiple wavelet coherence and bias-rectified wavelet spectra (MATLAB). City University of Hong Kong, Low-Carbon and Climate Impact Research Centre. ",
              tags$a(href = "https://www.cityu.edu.hk/lccic/wavelet/", target = "_blank", "cityu.edu.hk/lccic/wavelet")),
      tags$li("Koopmans, L. H. (1995). The spectral analysis of time series (2nd ed.). Academic Press. (multivariate partial and multiple coherency)"),
      tags$li("Priestley, M. B. (1981). Spectral analysis and time series. Academic Press."),
      tags$li("Liu, Y., Liang, X. S., & Weisberg, R. H. (2007). Rectification of the bias in the wavelet power spectrum. Journal of Atmospheric and Oceanic Technology, 24(12), 2093-2102."),
      tags$li("Veleda, D., Montagne, R., & Araujo, M. (2012). Cross-wavelet bias corrected by normalizing scales. Journal of Atmospheric and Oceanic Technology, 29(9), 1401-1408."),
      tags$li("Dima, B., Dima, S. M., & Barna, F. (2015). A wavelet analysis of capital markets' integration in Latin America. Applied Economics, 47(10), 1019-1036. ",
              tags$a(href = "https://doi.org/10.1080/00036846.2014.987916", target = "_blank", "doi:10.1080/00036846.2014.987916")),
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
  list(text = "Adjei, A. N. K., Tweneboah, G., & Owusu Junior, P. (2022). Interdependence of economic policy uncertainty and business cycles in selected emerging market economies. Journal of Financial Economic Policy, 14(5), 601-632. (wavelet multiple correlation and cross-correlation)",
       doi = "10.1108/JFEP-07-2021-0193"),
  list(text = "Kyei, C. B., Asafo-Adjei, E., Owusu Junior, P., Adam, A. M., Idun, A. A.-A., & Boateng, J. K. G. A. (2023). Nexus between commodities and banking sector financial soundness: The role of general macroeconomic setting in Ghana. Cogent Economics & Finance, 11(1), 2173871. (bi-wavelet, wavelet multiple and partial wavelet coherence)",
       doi = "10.1080/23322039.2023.2173871"),
  list(text = "Barson, Z., Owusu Junior, P., & Adam, A. M. (2023). Comovement between commodity returns in Ghana: The role of exchange rates. Journal of Economic Structures, 12, 17.",
       doi = "10.1186/s40008-023-00312-z"),
  list(text = "Owusu Junior, P., Tetteh, J. E., Nkrumah-Boadu, B., & Adjei, A. N. K. (2024). Comovement of African stock markets: Any influence from the COVID-19 pandemic? Heliyon, 10(9), e29409.",
       doi = "10.1016/j.heliyon.2024.e29409"),
  list(text = "Yeboah, S. D., Owusu Junior, P., Idun, A. A.-A., Akorsu, P. K., & Fumey, M. P. (2025). Dynamic interdependence of major currencies and the US dollar: A wavelet coherence approach. Future Business Journal, 11(1).",
       doi = "10.1186/s43093-025-00636-1"),
  list(text = "Tetteh, J. E., & Owusu Junior, P. (2025). Co-movement between stock markets in advanced economies and Africa in times of uncertainty: A time-frequency domain approach. PLOS ONE, 20(11), e0334325.",
       doi = "10.1371/journal.pone.0334325"),
  list(text = "Woode, J. K., Idun, A. A.-A., Kawor, S., Owusu Junior, P., & Adam, A. M. (2026). Commodity imports and exports dependence, market sentiment, and pressure currencies in Sub-Saharan Africa: A time-frequency wavelet analysis. Strategic Business Research, 2, 100054. (bi-wavelet, partial wavelet and multiple wavelet analysis)",
       doi = "10.1016/j.sbr.2026.100054")
)

shinyApp(ui, server)
