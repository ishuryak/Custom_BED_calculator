library(shiny)
library(bslib)
library(ggplot2)
library(DT)

source("R/bed_models.R", local = TRUE)

num <- function(id, label, value, min, max, step) {
  numericInput(id, label, value = value, min = min, max = max, step = step)
}

# Allowed range and readable name of every numeric input, and whether it must
# be a whole number. Checked on the server, because numericInput min/max do
# not stop typed values.
LIMITS <- list(
  m = list("Number of fractions", 1, 100, TRUE),
  d = list("Dose per fraction", 0.1, 30, FALSE),
  T = list("Overall treatment time", 1, 150, FALSE),
  ref_m = list("Reference number of fractions", 2, 100, TRUE),
  ref_d = list("Reference dose per fraction", 0.1, 30, FALSE),
  ref_T = list("Reference overall time", 1, 150, FALSE),
  alpha = list("Tumor α", 0.01, 2, FALSE),
  ab = list("Tumor α/β", 0.5, 100, FALSE),
  lambda_di = list("DI repopulation rate", 0, 2, FALSE),
  lambda_dd = list("DD maximum repopulation rate", 0, 5, FALSE),
  Tk = list("DI repopulation onset", 0, 100, FALSE),
  ab_e = list("Mucosa α/β", 0.5, 100, FALSE),
  alpha_e = list("Mucosa α", 0.01, 2, FALSE),
  Tp_e = list("Mucosal doubling time", 0.1, 50, FALSE),
  Tk_e = list("Mucosal repopulation onset", 0, 100, FALSE),
  ab_l = list("Late-tissue α/β", 0.5, 20, FALSE),
  tol_early = list("Allowed early BED above reference", 0, 30, FALSE),
  tol_late = list("Allowed late BED above reference", 0, 30, FALSE),
  min_gain = list("Minimum tumor BED gain", 0, 30, FALSE),
  min_days = list("Minimum days saved", 1, 140, TRUE),
  d_max = list("Highest dose per fraction searched", 2, 30, FALSE))

PRESET_LABELS <- c(
  "Calculator defaults (version 1)" = "calculator",
  "Withers et al. 1988: 0.6 Gy/day (2-Gy fractions) after 4 weeks" = "withers",
  "HYHOPE protocol: α 0.35/Gy, doubling time 3 d, onset day 21" = "hyhope",
  "Custom (edited values)" = "custom")

ui <- page_sidebar(
  fillable = FALSE,
  title = "Biologically effective dose (BED) calculator for head and neck squamous cell carcinoma (HNSCC), version 2",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  sidebar = sidebar(
    width = 350,
    accordion(
      open = c("Regimen to evaluate", "Tumor (DI and DD models)"),
      accordion_panel("Regimen to evaluate",
        num("m", "Number of fractions", 20, 1, 100, 1),
        num("d", "Dose per fraction (Gy)", 2.65, 0.1, 30, 0.05),
        checkboxInput("auto_T", "Compute overall time (once daily, Mon-Fri, starting Monday)", TRUE),
        conditionalPanel("!input.auto_T",
          num("T", "Overall treatment time (days)", 26, 1, 150, 1)),
        uiOutput("T_shown")),
      accordion_panel("Reference regimen",
        num("ref_m", "Number of fractions", 35, 2, 100, 1),
        num("ref_d", "Dose per fraction (Gy)", 2, 0.1, 30, 0.05),
        checkboxInput("ref_auto_T", "Compute overall time (once daily, Mon-Fri, starting Monday)", TRUE),
        conditionalPanel("!input.ref_auto_T",
          num("ref_T", "Overall treatment time (days)", 47, 1, 150, 1)),
        uiOutput("ref_T_shown"),
        helpText("All differences and limits are relative to this regimen.",
                 "The DD model is calibrated so that its tumor BED equals the DI model's here.")),
      accordion_panel("Tumor (DI and DD models)",
        selectInput("preset", "Parameter preset", PRESET_LABELS),
        num("alpha", "α, tumor cell killing (1/Gy)", 0.2, 0.01, 2, 0.01),
        num("ab", "α/β, tumor (Gy)", 10, 0.5, 100, 1),
        num("lambda_di", "Repopulation rate, DI model (1/day)", 0.2, 0, 2, 0.01),
        num("lambda_dd", "Maximum repopulation rate, DD model (1/day)", 0.5, 0, 5, 0.01),
        num("Tk", "Onset of accelerated repopulation, DI model (days)", 28, 0, 100, 1),
        uiOutput("tumor_note")),
      accordion_panel("Early toxicity (acute mucosa)",
        num("ab_e", "α/β, mucosa (Gy)", 10, 0.5, 100, 1),
        num("alpha_e", "α, mucosa (1/Gy)", 0.35, 0.01, 2, 0.01),
        num("Tp_e", "Mucosal doubling time during repopulation (days)", 2.5, 0.1, 50, 0.1),
        num("Tk_e", "Onset of mucosal repopulation (days)", 7, 0, 100, 1),
        uiOutput("lambda_e_shown"),
        checkboxInput("chemo", "Concurrent chemotherapy planned", FALSE),
        helpText("Lowers the Fowler tolerance zone shown on the Results tab, and the",
                 "early limit when the Fowler option is used for the search. It does",
                 "not change any BED, and with the default early limit the reference",
                 "regimen is assumed to have the same chemotherapy.")),
      accordion_panel("Late toxicity",
        num("ab_l", "α/β, late-responding tissue (Gy)", 3, 0.5, 20, 0.5),
        helpText("No repopulation. Late BED = n d (1 + d/(α/β)), which does",
                 "not depend on α, so no α is needed here. It is evaluated at",
                 "the prescription dose per fraction."))
    )
  ),
  navset_card_tab(
    nav_panel("Results",
      uiOutput("dd_warning"),
      uiOutput("neg_warning"),
      helpText("Tumor BED is shown under two repopulation models: dose-independent (DI),",
               "with a fixed onset and rate, and dose-dependent (DD), where onset and rate",
               "depend on cell killing. EQD2 is the equivalent dose in 2-Gy fractions.",
               "Model outputs for hypothesis generation, not treatment recommendations."),
      tableOutput("results"),
      uiOutput("verdicts"),
      uiOutput("acute_zone"),
      h5("Published regimens for comparison"),
      tableOutput("benchmarks"),
      uiOutput("hypno_note")),
    nav_panel("Regimens meeting the BED criteria",
      uiOutput("dd_warning2"),
      div(class = "alert alert-secondary",
          "These are model outputs for hypothesis generation, not treatment recommendations.",
          "Differences of a few Gy are within the uncertainty of the model parameters.",
          "Any trial regimen needs clinical precedent and toxicity stopping rules."),
      layout_columns(
        col_widths = c(4, 4, 4),
        selectInput("criterion", "Tumor BED gain required under",
                    c("DI or DD model" = "either", "both models" = "both",
                      "DI model" = "di", "DD model" = "dd")),
        selectInput("rank_by", "Rank by",
                    c("larger of ΔDI and ΔDD" = "best",
                      "smaller of ΔDI and ΔDD" = "worst",
                      "ΔBED DI" = "di", "ΔBED DD" = "dd")),
        sliderInput("m_range", "Number of fractions searched", 5, 34, c(5, 34), step = 1),
        num("min_gain", "Minimum tumor BED gain (Gy)", 1, 0, 30, 0.5),
        num("min_days", "Minimum days saved", 7, 1, 140, 1),
        selectInput("d_step", "Dose-per-fraction step (Gy)", c(0.05, 0.1, 0.2, 0.25, 0.5),
                    selected = 0.05),
        num("d_max", "Highest dose per fraction searched (Gy)", 8, 2, 30, 0.5),
        radioButtons("early_mode", "Early toxicity limit",
                     c("Reference regimen, plus allowed excess" = "reference",
                       "Lower edge of the Fowler tolerance zone" = "fowler")),
        conditionalPanel("input.early_mode == 'reference'",
          num("tol_early", "Allowed early BED above reference (Gy)", 0, 0, 30, 0.5)),
        num("tol_late", "Allowed late BED above reference (Gy)", 0, 0, 30, 0.5)),
      uiOutput("search_note"),
      uiOutput("hypno_cap_note"),
      uiOutput("suggest_msg"),
      DTOutput("suggest_table", fill = FALSE),
      br(),
      downloadButton("download", "Download table (CSV)"),
      plotOutput("suggest_plot", height = "380px")),
    nav_panel("Model and sources", uiOutput("about"))
  )
)

server <- function(input, output, session) {

  checked <- function(id) {
    v <- input[[id]]
    lim <- LIMITS[[id]]
    validate(need(is.numeric(v) && length(v) == 1 && is.finite(v),
                  paste("Enter a number for:", lim[[1]])))
    validate(need(v >= lim[[2]] && v <= lim[[3]],
                  sprintf("%s must be between %g and %g.", lim[[1]], lim[[2]], lim[[3]])),
             need(!lim[[4]] || v == round(v),
                  paste(lim[[1]], "must be a whole number.")))
    v
  }

  observeEvent(input$preset, {
    if (input$preset == "custom") return()
    pr <- tumor_presets()[[input$preset]]
    updateNumericInput(session, "alpha", value = pr$alpha)
    updateNumericInput(session, "lambda_di", value = pr$lambda_di)
    updateNumericInput(session, "Tk", value = pr$Tk)
  }, ignoreInit = TRUE)

  # Show which preset the current tumor values match, or "custom".
  observeEvent(list(input$alpha, input$lambda_di, input$Tk), {
    now <- preset_matching(input$alpha, input$lambda_di, input$Tk)
    if (!identical(now, input$preset)) updateSelectInput(session, "preset", selected = now)
  }, ignoreInit = TRUE)

  params <- reactive({
    list(alpha = checked("alpha"), ab = checked("ab"),
         lambda_di = checked("lambda_di"), lambda_dd = checked("lambda_dd"),
         Tk = checked("Tk"), ab_e = checked("ab_e"), alpha_e = checked("alpha_e"),
         lambda_e = log(2) / checked("Tp_e"), Tk_e = checked("Tk_e"),
         ab_l = checked("ab_l"))
  })

  # Overall time: from the calendar, or typed but not shorter than once-daily
  # Mon-Fri treatment allows.
  time_for <- function(m, auto, id) {
    if (isTRUE(auto)) return(overall_time(m))
    T <- checked(id)
    validate(need(T >= overall_time(m), sprintf(
      "%d fractions given once daily (Mon-Fri) need at least %d days. %g days was entered. This calculator covers once-daily Mon-Fri treatment.",
      m, overall_time(m), T)))
    T
  }
  ref <- reactive({
    m <- checked("ref_m")
    list(m = m, d = checked("ref_d"), T = time_for(m, input$ref_auto_T, "ref_T"))
  })
  regimen <- reactive({
    m <- checked("m")
    list(m = m, d = checked("d"), T = time_for(m, input$auto_T, "T"))
  })
  calibration <- reactive({
    p <- params()
    dd_calibration(ref(), p$alpha, p$ab, p$lambda_di, p$lambda_dd, p$Tk)
  })

  # Keep the fraction search range below the reference number of fractions.
  observeEvent(input$ref_m, {
    m <- input$ref_m
    if (is.numeric(m) && length(m) == 1 && is.finite(m) && m >= 2 && m == round(m)) {
      rng <- search_range_for_ref(m)
      updateSliderInput(session, "m_range", min = 1, max = rng[2], value = rng)
    }
  }, ignoreInit = TRUE)

  output$T_shown <- renderUI({
    r <- regimen()
    helpText(sprintf("Used: %d fractions, %g days, total dose %.2f Gy.", r$m, r$T, r$m * r$d))
  })
  output$ref_T_shown <- renderUI({
    r <- ref()
    helpText(sprintf("Used: %d fractions, %g days, total dose %.2f Gy.", r$m, r$T, r$m * r$d))
  })
  output$lambda_e_shown <- renderUI(
    helpText(sprintf("Mucosal repopulation rate = ln2 / doubling time = %.3f per day.",
                     params()$lambda_e)))
  output$tumor_note <- renderUI({
    p <- params()
    helpText(sprintf(paste("DI repopulation costs %.2f Gy of tumor BED per day after day %g.",
                           "The DI parameters also set the DD calibration."),
                     p$lambda_di / p$alpha, p$Tk))
  })

  dd_message <- reactive({
    cal <- calibration()
    if (cal$status != "infeasible") return(NULL)
    onset <- if (is.finite(cal$onset_ref)) sprintf(" (it would have to start on day %.1f)", cal$onset_ref) else ""
    div(class = "alert alert-warning",
        paste0("DD model unavailable: with these parameters, no repopulation onset after the ",
               "start of treatment makes the DD tumor BED equal the DI tumor BED for the ",
               "reference regimen", onset, ". Increase the DD maximum repopulation rate or ",
               "change the DI parameters. DD values are hidden and DD criteria are switched off."))
  })
  output$dd_warning <- renderUI(dd_message())
  output$dd_warning2 <- renderUI(dd_message())

  both <- reactive({
    r <- regimen(); rf <- ref()
    rbind(compute_all(r$m, r$d, r$T, params(), rf),
          compute_all(rf$m, rf$d, rf$T, params(), rf))
  })

  fmt <- function(x, plus = FALSE) {
    ifelse(is.na(x), "not available", sprintf(if (plus) "%+.2f" else "%.2f", x))
  }

  output$results <- renderTable({
    b <- both()
    rows <- c("Overall time (days)" = "T", "Total dose (Gy)" = "total_dose",
              "Tumor BED, DI model (Gy)" = "bed_di", "Tumor BED, DD model (Gy)" = "bed_dd",
              "Early toxicity BED (Gy)" = "bed_early", "Early toxicity EQD2 (Gy)" = "eqd2_early",
              "Late toxicity BED (Gy)" = "bed_late", "Late toxicity EQD2 (Gy)" = "eqd2_late")
    v1 <- unlist(b[1, rows]); v2 <- unlist(b[2, rows])
    out <- data.frame(Quantity = names(rows), `This regimen` = fmt(v1), Reference = fmt(v2),
                      Difference = fmt(v1 - v2, plus = TRUE), check.names = FALSE)
    out[1, 2:4] <- c(sprintf("%g", v1[1]), sprintf("%g", v2[1]), sprintf("%+g", v1[1] - v2[1]))
    out
  }, striped = TRUE, bordered = TRUE, align = "lrrr")

  output$verdicts <- renderUI({
    b <- both()
    word <- function(x, y) {
      if (is.na(x)) "not available" else if (x > y + 1e-9) "higher" else if (x < y - 1e-9) "lower" else "the same"
    }
    tags$p(class = "text-muted",
      sprintf(paste("Compared with the reference (model values): tumor BED is %s under DI",
                    "and %s under DD. Early toxicity BED is %s, and late toxicity BED is %s."),
              word(b$bed_di[1], b$bed_di[2]), word(b$bed_dd[1], b$bed_dd[2]),
              word(b$bed_early[1], b$bed_early[2]), word(b$bed_late[1], b$bed_late[2])))
  })

  # Unit label for early BED: Gy10 only when the mucosal alpha/beta is 10 Gy.
  early_unit <- reactive({
    ab_e <- params()$ab_e
    if (isTRUE(all.equal(ab_e, 10))) "Gy₁₀" else sprintf("Gy (α/β = %g Gy)", ab_e)
  })

  output$neg_warning <- renderUI({
    b <- both()
    v <- unlist(b[, c("bed_di", "bed_dd", "bed_early")])
    if (any(v < 0, na.rm = TRUE))
      div(class = "alert alert-warning",
          "Some BEDs are negative: with these parameters the repopulation loss exceeds the cell-killing term.",
          "Check the inputs; comparisons between negative BEDs are not meaningful.")
  })

  fowler_defaults <- reactive({
    p <- params()
    isTRUE(all.equal(c(p$ab_e, p$alpha_e, p$lambda_e, p$Tk_e), c(10, 0.35, log(2) / 2.5, 7)))
  })

  output$acute_zone <- renderUI({
    z <- fowler_zone(isTRUE(input$chemo))
    b <- both()
    u <- early_unit()
    txt <- sprintf(paste("Acute tolerance zone of Fowler et al. (2003): %g–%g Gy₁₀",
                         "%s. This regimen: %.1f %s. Reference: %.1f %s."),
                   z[1], z[2],
                   if (isTRUE(input$chemo)) "with concurrent chemotherapy (derived from Fowler's possible decrement of 3–5 Gy₁₀)" else "for radiotherapy alone (Fowler 2010 later gave a 59–63 Gy₁₀ gray zone)",
                   b$bed_early[1], u, b$bed_early[2], u)
    if (!fowler_defaults())
      txt <- paste(txt, "The zone applies only with Fowler's acute parameters (the defaults), which have been changed.")
    tags$p(class = "text-muted", txt)
  })

  output$benchmarks <- renderTable({
    bm <- benchmark_regimens(); rf <- ref()
    a <- compute_all(bm$m, bm$d, bm$T, params(), rf)
    r0 <- compute_all(rf$m, rf$d, rf$T, params(), rf)
    data.frame(Regimen = bm$name, Days = sprintf("%g", bm$T),
               `Δ BED DI` = fmt(a$bed_di - r0$bed_di, TRUE),
               `Δ BED DD` = fmt(a$bed_dd - r0$bed_dd, TRUE),
               `Early BED` = fmt(a$bed_early),
               `Δ early BED` = fmt(a$bed_early - r0$bed_early, TRUE),
               `Δ late BED` = fmt(a$bed_late - r0$bed_late, TRUE), check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)

  output$hypno_note <- renderUI({
    bm <- benchmark_regimens(); rf <- ref()
    a <- compute_all(bm$m, bm$d, bm$T, params(), rf)
    dd_txt <- if (is.na(a$bed_dd[1])) "not available" else sprintf("%+.1f Gy", a$bed_dd[1] - a$bed_dd[2])
    helpText(sprintf(paste("HYPNO randomized 792 patients. 55 Gy in 20 fractions was non-inferior to",
                           "66 Gy in 33 fractions for locoregional control and grade 3+ late effects.",
                           "With the current parameters, test arm minus control arm in tumor BED is",
                           "%+.1f Gy under DI and %s under DD. Overall times assume a Monday start."),
                     a$bed_di[1] - a$bed_di[2], dd_txt))
  })

  output$hypno_cap_note <- renderUI({
    sr <- search(); rf <- ref()
    if (isTRUE(sr$blocked)) return(NULL)
    bm <- benchmark_regimens()[1, ]
    e <- compute_all(bm$m, bm$d, bm$T, params(), rf)$bed_early
    r0 <- compute_all(rf$m, rf$d, rf$T, params(), rf)$bed_early
    cap <- if (is.null(sr$cap)) r0 + checked("tol_early") else sr$cap
    u <- early_unit()
    helpText(sprintf(paste("For comparison, HYPNO's 55 Gy in 20 fractions has early BED %.1f %s,",
                           "which is %s the current early limit of %.1f %s."),
                     e, u, if (e > cap + 1e-9) "above" else "within", cap, u))
  })

  search <- reactive({
    rf <- ref()
    validate(need(input$m_range[1] < rf$m,
                  "The fraction range must include numbers below the reference."))
    m_hi <- min(input$m_range[2], rf$m - 1)
    if (input$early_mode == "fowler" && !fowler_defaults()) {
      s <- empty_suggestions()
      attr(s, "reason") <- "fowler_params_changed"
      return(list(s = s, cap = NULL, blocked = TRUE))
    }
    cap <- if (input$early_mode == "fowler") fowler_zone(isTRUE(input$chemo))[1] else NULL
    s <- suggest_regimens(params(), rf, m_range = input$m_range[1]:m_hi,
                          d_step = as.numeric(input$d_step), d_max = checked("d_max"),
                          criterion = input$criterion, rank_by = input$rank_by,
                          tol_early = if (is.null(cap)) checked("tol_early") else 0,
                          tol_late = checked("tol_late"), early_cap = cap,
                          min_gain = checked("min_gain"), min_days_saved = checked("min_days"))
    list(s = s, cap = cap, blocked = FALSE)
  })
  suggestions <- reactive(search()$s)

  output$search_note <- renderUI({
    p <- params(); sr <- search()
    if (isTRUE(sr$blocked)) return(NULL)
    early_txt <- if (is.null(sr$cap)) "early BED at most the reference's plus the allowed excess" else
      sprintf("early BED at most %g Gy₁₀ (Fowler zone, lower edge%s)", sr$cap,
              if (isTRUE(input$chemo)) ", with concurrent chemotherapy" else "")
    helpText(sprintf(paste("Searched: once-daily Mon-Fri regimens with a dose per fraction above the",
                           "reference's (multiples of %s Gy, up to %g Gy), %s, and late BED at most the",
                           "reference's plus the allowed excess. For each number of fractions, the",
                           "highest-ranked dose per fraction is listed. Δ = regimen minus reference.",
                           "The list depends strongly on the tumor repopulation settings (now %.2f Gy/day under DI)."),
                     input$d_step, checked("d_max"), early_txt, p$lambda_di / p$alpha))
  })

  pretty_suggestions <- reactive({
    s <- suggestions()
    dd_col <- if (nrow(s) > 0 && all(is.na(s$d_dd))) rep("not available", nrow(s)) else s$d_dd
    data.frame(Fractions = s$m, `Gy per fraction` = s$d, Days = s$T, `Weeks saved` = s$weeks_saved,
               `Total Gy` = s$total_dose,
               `Δ BED DI` = s$d_di, `Δ BED DD` = dd_col,
               `Δ early BED` = s$d_early, `Δ late BED` = s$d_late,
               `Early EQD2` = s$eqd2_early, `Late EQD2` = s$eqd2_late,
               `Next dose step would exceed` = s$next_step, check.names = FALSE)
  })

  output$suggest_msg <- renderUI({
    s <- suggestions()
    if (nrow(s) > 0) return(NULL)
    reason <- attr(s, "reason")
    if (identical(reason, "dd_unavailable"))
      return(div(class = "alert alert-info",
                 "The DD model is unavailable with these parameters, so DD-based criteria and rankings",
                 "(including 'smaller of ΔDI and ΔDD') cannot be used."))
    if (identical(reason, "empty_dose_range"))
      return(div(class = "alert alert-info",
                 "No dose per fraction lies between the reference dose per fraction and the highest dose searched."))
    if (identical(reason, "fowler_params_changed"))
      return(div(class = "alert alert-info",
                 "The Fowler tolerance zone applies only with Fowler's acute parameters",
                 "(the defaults: α/β 10 Gy, α 0.35/Gy, doubling time 2.5 d, onset day 7).",
                 "Restore them or use the reference-regimen limit."))
    if (identical(reason, "no_shorter_regimen"))
      return(div(class = "alert alert-info",
                 "No number of fractions in the searched range saves the minimum number of days."))
    div(class = "alert alert-info",
        "No regimen in the searched range meets all limits with these settings.",
        if (input$criterion != "either") "Options: require the gain under 'DI or DD model'," else "Options:",
        "lower the minimum gain or days saved, use the Fowler early limit,",
        "or allow some early or late BED above the reference.")
  })

  output$suggest_table <- renderDT({
    s <- pretty_suggestions()
    two_dp <- c("Gy per fraction", "Total Gy", "Δ BED DI", "Δ BED DD", "Δ early BED",
                "Δ late BED", "Early EQD2", "Late EQD2")
    two_dp <- two_dp[vapply(s[two_dp], is.numeric, logical(1))]
    datatable(s, rownames = FALSE, fillContainer = FALSE,
              options = list(pageLength = 15, dom = "tip")) |>
      formatRound(two_dp, 2) |>
      formatRound("Weeks saved", 1)
  })

  output$download <- downloadHandler(
    filename = function() sprintf("bed_v2_regimens_%s.csv", format(Sys.Date())),
    content = function(file) {
      s <- pretty_suggestions()
      req(nrow(s) > 0)
      num_cols <- vapply(s, is.numeric, logical(1))
      s[num_cols] <- lapply(s[num_cols], round, 3)
      names(s) <- sub("Δ", "Delta", names(s))
      write.csv(s, file, row.names = FALSE)
    })

  output$suggest_plot <- renderPlot({
    s <- suggestions()
    validate(need(nrow(s) > 0, ""))
    long <- rbind(data.frame(m = s$m, delta = s$d_di, Model = "DI"),
                  data.frame(m = s$m, delta = s$d_dd, Model = "DD"))
    long <- long[!is.na(long$delta), ]
    ggplot(long, aes(m, delta, colour = Model, shape = Model, linetype = Model)) +
      geom_hline(yintercept = 0, colour = "grey40") +
      geom_line() + geom_point(size = 3) +
      scale_colour_manual(values = c(DI = "#1f77b4", DD = "#d62728")) +
      scale_shape_manual(values = c(DI = 16, DD = 17)) +
      scale_linetype_manual(values = c(DI = "solid", DD = "dashed")) +
      labs(x = "Number of fractions (listed dose per fraction for each)",
           y = "Tumor BED minus reference (Gy)",
           title = "Tumor BED difference from the reference for the listed regimens") +
      theme_minimal(base_size = 15) +
      theme(axis.title = element_text(face = "bold"), legend.position = "top")
  })

  output$about <- renderUI(withMathJax(HTML(
    markdown::mark(file = "about.md", options = "+latex_math"))))
}

shinyApp(ui, server)
