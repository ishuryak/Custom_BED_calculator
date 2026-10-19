library(shiny)
library(bslib)
library(ggplot2)
library(DT)

source("R/bed_models.R", local = TRUE)
source("R/export.R", local = TRUE)

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
  ref_m = list("Comparison number of fractions", 2, 100, TRUE),
  ref_d = list("Comparison dose per fraction", 0.1, 30, FALSE),
  ref_T = list("Comparison overall time", 1, 150, FALSE),
  anchor_m = list("Anchor number of fractions", 2, 100, TRUE),
  anchor_d = list("Anchor dose per fraction", 0.1, 30, FALSE),
  alpha = list("Tumor α", 0.01, 2, FALSE),
  ab = list("Tumor α/β", 0.5, 100, FALSE),
  lambda_di = list("DI repopulation rate", 0, 2, FALSE),
  lambda_dd = list("DD maximum repopulation rate", 0, 5, FALSE),
  Tk = list("DI repopulation onset", 0, 100, FALSE),
  g = list("Background tumor growth rate", 0, 0.1, FALSE),
  ab_e = list("Mucosa α/β", 0.5, 100, FALSE),
  alpha_e = list("Mucosa α", 0.01, 2, FALSE),
  Tp_e = list("Mucosal doubling time", 0.1, 50, FALSE),
  Tk_e = list("Mucosal repopulation onset", 0, 100, FALSE),
  ab_l = list("Late-tissue α/β", 0.5, 20, FALSE),
  tol_early = list("Allowed early BED above the comparison", 0, 30, FALSE),
  tol_late = list("Allowed late BED above the comparison", 0, 30, FALSE),
  min_gain = list("Minimum tumor BED change", -20, 30, FALSE),
  min_days = list("Minimum days saved", 1, 140, TRUE),
  d_max = list("Highest dose per fraction searched", 2, 30, FALSE))

PRESET_LABELS <- c(
  "Calculator defaults (version 1)" = "calculator",
  "Withers et al. 1988: 0.6 Gy/day in 2-Gy fractions after 4 weeks (follows α/β)" = "withers",
  "Hartley et al. 2011: α 0.3/Gy, doubling time 3 d, onset day 22 (elapsed)" = "hartley",
  "HYHOPE protocol: α 0.35/Gy, doubling time 3 d, onset day 21" = "hyhope",
  "Scenario D18 (2026-10-18 proposal, not fitted)" = "d18",
  "Custom (edited values)" = "custom")

# Short names for the sensitivity table.
PRESET_SHORT <- c(calculator = "Version 1 defaults", withers = "Withers 1988 (rule)",
                  hartley = "Hartley 2011", hyhope = "HYHOPE protocol", d18 = "Scenario D18")

COMPARATOR_LABELS <- c(
  "Enter below (default 35 × 2 Gy)" = "custom",
  "HYPNO test arm: 55 Gy in 20 fractions (26 days)" = "hypno_test",
  "HYPNO control arm: 66 Gy in 33 fractions, 6 fractions/week (38 days)" = "hypno_control")

ui <- page_sidebar(
  fillable = FALSE,
  title = "Biologically effective dose (BED) calculator for head and neck squamous cell carcinoma (HNSCC), version 2.2",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  sidebar = sidebar(
    width = 360,
    accordion(
      open = c("Regimen to evaluate", "Tumor (DI and DD models)"),
      accordion_panel("Regimen to evaluate",
        num("m", "Number of fractions", 20, 1, 100, 1),
        num("d", "Dose per fraction (Gy)", 2.65, 0.1, 30, 0.05),
        checkboxInput("auto_T", "Compute overall time (once daily, Mon-Fri, starting Monday)", TRUE),
        conditionalPanel("!input.auto_T",
          num("T", "Overall treatment time (days)", 26, 1, 150, 1)),
        uiOutput("T_shown")),
      accordion_panel("Comparison regimen",
        selectInput("comparator", "Compare with", COMPARATOR_LABELS),
        conditionalPanel("input.comparator == 'custom'",
          num("ref_m", "Number of fractions", 35, 2, 100, 1),
          num("ref_d", "Dose per fraction (Gy)", 2, 0.1, 30, 0.05),
          checkboxInput("ref_auto_T", "Compute overall time (once daily, Mon-Fri, starting Monday)", TRUE),
          conditionalPanel("!input.ref_auto_T",
            num("ref_T", "Overall treatment time (days)", 47, 1, 150, 1))),
        uiOutput("ref_T_shown"),
        helpText("All differences and toxicity limits are relative to this regimen.",
                 "It does not change the DD model, which is calibrated at the anchor below.")),
      accordion_panel("Tumor (DI and DD models)",
        selectInput("preset", "Parameter preset", PRESET_LABELS),
        num("alpha", "α, tumor cell killing (1/Gy)", 0.2, 0.01, 2, 0.01),
        num("ab", "α/β, tumor (Gy)", 10, 0.5, 100, 1),
        num("lambda_di", "Repopulation rate, DI model (1/day)", 0.2, 0, 2, 0.01),
        num("lambda_dd", "Maximum repopulation rate, DD model (1/day)", 0.5, 0, 5, 0.01),
        num("Tk", "Onset of accelerated repopulation, DI model (days)", 28, 0, 100, 1),
        num("g", "Background tumor growth rate g (1/day, 0 = off)", 0, 0, 0.1, 0.001),
        uiOutput("tumor_note")),
      accordion_panel("DD calibration anchor",
        num("anchor_m", "Anchor number of fractions", 35, 2, 100, 1),
        num("anchor_d", "Anchor dose per fraction (Gy)", 2, 0.1, 30, 0.05),
        checkboxInput("coupled", "Legacy: calibrate at the comparison regimen instead (versions 1 and 2.0)", FALSE),
        uiOutput("anchor_note")),
      accordion_panel("Early toxicity (acute mucosa)",
        num("ab_e", "α/β, mucosa (Gy)", 10, 0.5, 100, 1),
        num("alpha_e", "α, mucosa (1/Gy)", 0.35, 0.01, 2, 0.01),
        num("Tp_e", "Mucosal doubling time during repopulation (days)", 2.5, 0.1, 50, 0.1),
        num("Tk_e", "Onset of mucosal repopulation (days, elapsed)", 7, 0, 100, 1),
        uiOutput("lambda_e_shown"),
        checkboxInput("chemo", "Concurrent chemotherapy planned", FALSE),
        helpText("Affects acute mucositis tolerance only: it lowers the Fowler zone shown",
                 "on the Results tab, and the early limit when the Fowler option is used",
                 "in the search. No BED (tumor, early or late) includes an effect of",
                 "chemotherapy. A shorter course usually means fewer chemotherapy cycles",
                 "(HYPNO: 4 vs 5 weekly courses), which this calculator does not model.")),
      accordion_panel("Late toxicity",
        num("ab_l", "α/β, late-responding tissue (Gy)", 3, 0.5, 20, 0.5),
        helpText("No repopulation. Late BED = n d (1 + d/(α/β)), which does",
                 "not depend on α. It is evaluated at the prescription dose per",
                 "fraction, not at the dose to any organ at risk."))
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
    nav_panel("Candidate regimens (model screen)",
      uiOutput("dd_warning2"),
      div(class = "alert alert-secondary",
          "These are model outputs for hypothesis generation, not treatment recommendations.",
          "Differences of a few Gy are within the uncertainty of the model parameters.",
          "A listed regimen is not shown to be safe or effective. Any trial regimen needs",
          "clinical precedent, organ-at-risk planning and toxicity stopping rules."),
      layout_columns(
        col_widths = c(4, 4, 4),
        selectInput("criterion", "Tumor BED criterion applied under",
                    c("DI or DD model" = "either", "both models" = "both",
                      "DI model" = "di", "DD model" = "dd")),
        selectInput("rank_by", "Rank by",
                    c("larger of ΔDI and ΔDD (optimistic)" = "best",
                      "smaller of ΔDI and ΔDD (conservative)" = "worst",
                      "ΔBED DI" = "di", "ΔBED DD" = "dd")),
        sliderInput("m_range", "Number of fractions searched", 5, 34, c(5, 34), step = 1),
        num("min_gain", "Minimum tumor BED change (Gy; negative = allowed loss)", 1, -20, 30, 0.5),
        num("min_days", "Minimum days saved", 7, 1, 140, 1),
        selectInput("d_step", "Dose-per-fraction step (Gy)", c(0.05, 0.1, 0.2, 0.25, 0.5),
                    selected = 0.05),
        num("d_max", "Highest dose per fraction searched (Gy)", 5, 2, 30, 0.5),
        radioButtons("early_mode", "Early toxicity limit",
                     c("Comparison regimen, plus allowed excess" = "reference",
                       "Lower edge of the Fowler tolerance zone" = "fowler",
                       "Acute BED of HYPNO 55 Gy/20 fractions" = "hypno",
                       "Fenwick et al. 2008 boundary (cell-kill BED, up to 50 days)" = "fenwick")),
        conditionalPanel("input.early_mode == 'reference'",
          num("tol_early", "Allowed early BED above the comparison (Gy)", 0, 0, 30, 0.5)),
        radioButtons("late_mode", "Late toxicity limit",
                     c("Comparison regimen, plus allowed excess" = "reference",
                       "Late BED of HYPNO 55 Gy/20 fractions" = "hypno")),
        conditionalPanel("input.late_mode == 'reference'",
          num("tol_late", "Allowed late BED above the comparison (Gy)", 0, 0, 30, 0.5))),
      uiOutput("search_note"),
      uiOutput("hypno_cap_note"),
      uiOutput("suggest_msg"),
      DTOutput("suggest_table", fill = FALSE),
      br(),
      downloadButton("download", "Download results and scenario (ZIP)"),
      plotOutput("suggest_plot", height = "380px")),
    nav_panel("Sensitivity",
      helpText("Tumor BED of the regimen being evaluated minus the comparison regimen, under",
               "each tumor preset and three tumor α/β values. These are selected scenarios,",
               "not a probability distribution: a count of scenarios is not a chance of benefit.",
               "Presets that set only DI values keep the DD maximum rate entered in the sidebar."),
      uiOutput("sens_summary"),
      tableOutput("sens_table"),
      h5("Toxicity BED differences (do not depend on the tumor scenario)"),
      tableOutput("sens_tox")),
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
  # Inputs added in version 2.1 fall back to their defaults when absent.
  checked_or <- function(id, default) if (is.null(input[[id]])) default else checked(id)
  choice_or <- function(id, default) if (is.null(input[[id]])) default else input[[id]]

  # Presets. The Withers rule follows the tumor alpha/beta.
  last_ab <- reactiveVal(10)
  observeEvent(input$preset, {
    if (input$preset == "custom") return()
    ab <- input$ab
    pr <- resolve_preset(input$preset, if (is.numeric(ab) && length(ab) == 1 && is.finite(ab)) ab else 10)
    for (k in intersect(names(pr), c("alpha", "ab", "lambda_di", "lambda_dd", "Tk")))
      updateNumericInput(session, k, value = pr[[k]])
  }, ignoreInit = TRUE)

  # Show which preset the current tumor values match, or "custom".
  observeEvent(list(input$alpha, input$lambda_di, input$Tk, input$ab, input$lambda_dd), {
    ab <- input$ab
    # While a field is being retyped it is empty (NA); keep the preset as it is.
    vals <- c(input$alpha, input$lambda_di, input$Tk, ab)
    if (length(vals) < 4 || any(!is.finite(vals))) return()
    if (identical(input$preset, "withers") &&
        !isTRUE(all.equal(ab, last_ab()))) {
      last_ab(ab)
      updateNumericInput(session, "lambda_di", value = withers_lambda(input$alpha, ab))
      return()
    }
    if (is.numeric(ab) && is.finite(ab)) last_ab(ab)
    now <- preset_matching(input$alpha, input$lambda_di, input$Tk, ab = ab, lambda_dd = input$lambda_dd)
    if (!identical(now, input$preset)) updateSelectInput(session, "preset", selected = now)
  }, ignoreInit = TRUE)

  params <- reactive({
    list(alpha = checked("alpha"), ab = checked("ab"),
         lambda_di = checked("lambda_di"), lambda_dd = checked("lambda_dd"),
         Tk = checked("Tk"), g = checked_or("g", 0),
         ab_e = checked("ab_e"), alpha_e = checked("alpha_e"),
         lambda_e = log(2) / checked("Tp_e"), Tk_e = checked("Tk_e"),
         clock_e = "elapsed", ab_l = checked("ab_l"))
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
    cmp <- choice_or("comparator", "custom")
    if (cmp != "custom") {
      r <- comparator_presets()[[cmp]]
      return(list(m = r$m, d = r$d, T = r$T, label = r$label))
    }
    m <- checked("ref_m")
    list(m = m, d = checked("ref_d"), T = time_for(m, input$ref_auto_T, "ref_T"))
  })
  regimen <- reactive({
    m <- checked("m")
    list(m = m, d = checked("d"), T = time_for(m, input$auto_T, "T"))
  })
  coupled <- reactive(isTRUE(input$coupled))
  anchor <- reactive({
    if (coupled()) { r <- ref(); return(list(m = r$m, d = r$d, T = r$T)) }
    m <- checked_or("anchor_m", 35)
    list(m = m, d = checked_or("anchor_d", 2), T = overall_time(m))
  })
  calibration <- reactive({
    p <- params(); a <- anchor()
    dd_calibration(a, p$alpha, p$ab, p$lambda_di, p$lambda_dd, p$Tk)
  })

  # Keep the fraction search range below the comparison number of fractions.
  observeEvent(list(input$ref_m, input$comparator), {
    m <- tryCatch(ref()$m, error = function(e) NA)
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
    helpText(sprintf(paste("Mucosal repopulation rate = ln2 / doubling time = %.3f per day.",
                           "Overall time is converted to Fowler's elapsed days (calendar days minus 1)."),
                     params()$lambda_e)))
  output$tumor_note <- renderUI({
    p <- params()
    helpText(sprintf(paste("DI repopulation costs %.2f Gy of tumor BED per day after day %g",
                           "(%.2f Gy per day in 2-Gy-equivalent terms).%s"),
                     p$lambda_di / p$alpha, p$Tk, p$lambda_di / p$alpha / (1 + 2 / p$ab),
                     if (p$g > 0) sprintf(" Background growth costs a further %.3f Gy per day over the whole course.", p$g / p$alpha) else ""))
  })
  output$anchor_note <- renderUI({
    a <- anchor(); cal <- calibration()
    where <- sprintf("%d × %g Gy in %g days", a$m, a$d, a$T)
    txt <- switch(cal$status,
      ok = sprintf("DD and DI tumor BEDs are equal at the anchor (%s). DD repopulation there starts on day %.1f.", where, cal$onset_ref),
      no_repopulation = sprintf("No DI repopulation, so DD has none either (anchor %s).", where),
      infeasible = sprintf("No DD calibration exists at the anchor (%s).", where))
    helpText(txt, if (coupled()) "Legacy mode: changing the comparison regimen changes the DD model." else
      "Changing the comparison regimen does not change the DD model.")
  })

  dd_message <- reactive({
    cal <- calibration()
    if (cal$status != "infeasible") return(NULL)
    onset <- if (is.finite(cal$onset_ref)) sprintf(" (it would have to start on day %.1f)", cal$onset_ref) else ""
    div(class = "alert alert-warning",
        paste0("DD model unavailable: with these parameters, no repopulation onset after the ",
               "start of treatment makes the DD tumor BED equal the DI tumor BED for the ",
               "anchor regimen", onset, ". Increase the DD maximum repopulation rate or ",
               "change the DI parameters. DD values are hidden and DD criteria are switched off."))
  })
  output$dd_warning <- renderUI(dd_message())
  output$dd_warning2 <- renderUI(dd_message())

  both <- reactive({
    r <- regimen(); rf <- ref(); a <- anchor()
    rbind(compute_all(r$m, r$d, r$T, params(), a),
          compute_all(rf$m, rf$d, rf$T, params(), a))
  })

  fmt <- function(x, plus = FALSE) {
    ifelse(is.na(x), "not available", sprintf(if (plus) "%+.2f" else "%.2f", x))
  }

  output$results <- renderTable({
    b <- both(); r <- regimen(); rf <- ref()
    rows <- c("Overall time (days)" = "T", "Total dose (Gy)" = "total_dose",
              "Tumor BED, DI model (Gy)" = "bed_di", "Tumor BED, DD model (Gy)" = "bed_dd",
              "Early toxicity BED (Gy)" = "bed_early", "Early toxicity EQD2 (Gy)" = "eqd2_early",
              "Late toxicity BED (Gy)" = "bed_late", "Late toxicity EQD2 (Gy)" = "eqd2_late")
    v1 <- unlist(b[1, rows]); v2 <- unlist(b[2, rows])
    out <- data.frame(Quantity = names(rows), `This regimen` = fmt(v1), Comparison = fmt(v2),
                      Difference = fmt(v1 - v2, plus = TRUE), check.names = FALSE)
    out[1, 2:4] <- c(sprintf("%g", v1[1]), sprintf("%g", v2[1]), sprintf("%+g", v1[1] - v2[1]))
    ab_l <- params()$ab_l
    for (ab in setdiff(c(2, 4), ab_l)) {
      x <- bed_late(r$m, r$d, ab); y <- bed_late(rf$m, rf$d, ab)
      out <- rbind(out, data.frame(Quantity = sprintf("Late toxicity BED at α/β %g Gy (sensitivity)", ab),
                                   `This regimen` = fmt(x), Comparison = fmt(y),
                                   Difference = fmt(x - y, TRUE), check.names = FALSE))
    }
    out
  }, striped = TRUE, bordered = TRUE, align = "lrrr")

  output$verdicts <- renderUI({
    b <- both()
    word <- function(x, y) {
      if (is.na(x)) "not available" else if (x > y + 1e-9) "higher" else if (x < y - 1e-9) "lower" else "the same"
    }
    tags$p(class = "text-muted",
      sprintf(paste("Compared with the comparison regimen (model values): tumor BED is %s under DI",
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
          "Check the inputs; comparisons between negative BEDs are not meaningful, and the search",
          "does not use a model whose comparison BED is negative.")
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
                         "%s. This regimen: %.1f %s. Comparison: %.1f %s."),
                   z[1], z[2],
                   if (isTRUE(input$chemo)) "with concurrent chemotherapy (derived from Fowler's possible decrement of 3–5 Gy₁₀)" else "for radiotherapy alone (Fowler 2010 later gave a 59–63 Gy₁₀ gray zone)",
                   b$bed_early[1], u, b$bed_early[2], u)
    if (!fowler_defaults())
      txt <- paste(txt, "The zone applies only with Fowler's acute parameters (the defaults), which have been changed.")
    r <- regimen(); rf <- ref()
    fw <- function(m, d, T) {
      lim <- fenwick_limit(T)
      if (is.na(lim)) return("outside its 50-day range")
      ck <- fenwick_cell_kill(m, d)
      sprintf("%.1f vs %.1f Gy₁₀, %s", ck, lim, if (ck <= lim + 1e-9) "below" else "above")
    }
    txt2 <- sprintf(paste("Fenwick et al. (2008) boundary on cell-kill BED (α/β 10 Gy, no repopulation",
                          "term, chemotherapy not included): this regimen %s; comparison %s."),
                    fw(r$m, r$d, r$T), fw(rf$m, rf$d, rf$T))
    tagList(tags$p(class = "text-muted", txt), tags$p(class = "text-muted", txt2))
  })

  bench <- reactive({
    bm <- benchmark_regimens(); rf <- ref(); a <- anchor(); p <- params()
    x <- compute_all(bm$m, bm$d, bm$T, p, a)
    r0 <- compute_all(rf$m, rf$d, rf$T, p, a)
    list(bm = bm, x = x, r0 = r0)
  })

  output$benchmarks <- renderTable({
    z <- bench(); x <- z$x; r0 <- z$r0
    data.frame(Regimen = z$bm$name, Days = sprintf("%g", z$bm$T),
               `Δ BED DI` = fmt(x$bed_di - r0$bed_di, TRUE),
               `Δ BED DD` = fmt(x$bed_dd - r0$bed_dd, TRUE),
               `Early BED` = fmt(x$bed_early),
               `Δ early BED` = fmt(x$bed_early - r0$bed_early, TRUE),
               `Late BED` = fmt(x$bed_late),
               `Δ late BED` = fmt(x$bed_late - r0$bed_late, TRUE), check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)

  output$hypno_note <- renderUI({
    x <- bench()$x
    dd_txt <- if (is.na(x$bed_dd[1])) "not available" else sprintf("%+.1f Gy", x$bed_dd[1] - x$bed_dd[2])
    tagList(
      helpText(sprintf(paste("Δ values are relative to the comparison regimen chosen in the sidebar.",
                             "HYPNO compared its test arm with its own control arm (66 Gy in 33 fractions",
                             "at 6 fractions/week), not with 70 Gy in 35 fractions. With the current",
                             "parameters, test arm minus control arm in tumor BED is %+.1f Gy under DI",
                             "and %s under DD."),
                       x$bed_di[1] - x$bed_di[2], dd_txt)),
      helpText(paste("HYPNO (792 patients, 75% with weekly cisplatin 35 mg/m², HPV not tested):",
                     "3-year locoregional control 50.7% vs 51.2%; grade 3+ late effects 18.8% vs 20.2%;",
                     "both passed a 10-percentage-point non-inferiority margin. Grade 3+ acute mucositis",
                     "50.9% (190/373) vs 54.7% (208/380). The margin does not support escalation beyond",
                     "55 Gy in 20 fractions.")),
      helpText(paste("HYHOPE (phase 1, 24 patients, favorable-risk HPV-positive oropharynx, reduced",
                     "elective volumes, 3 × weekly cisplatin 40 mg/m², designed against 60 Gy in 30",
                     "fractions): grade 3 oral mucositis 5/24, any acute grade 3 7/24, 2-year",
                     "locoregional recurrence 13.8% (95% CI 3.2–32%). Not a comparison for unselected HNSCC.")))
  })

  hypno_caps <- reactive({
    a <- compute_all(20, 2.75, 26, params(), anchor())
    list(early = a$bed_early, late = a$bed_late)
  })

  output$hypno_cap_note <- renderUI({
    sr <- search(); rf <- ref()
    if (isTRUE(sr$blocked) || input$early_mode %in% c("hypno", "fenwick")) return(NULL)
    e <- hypno_caps()$early
    r0 <- compute_all(rf$m, rf$d, rf$T, params(), anchor())$bed_early
    cap <- if (is.null(sr$cap)) r0 + checked("tol_early") else sr$cap
    u <- early_unit()
    helpText(sprintf(paste("For comparison, HYPNO's 55 Gy in 20 fractions has early BED %.1f %s,",
                           "which is %s the current early limit of %.1f %s."),
                     e, u, if (e > cap + 1e-9) "above" else "within", cap, u))
  })

  search_args <- reactive({
    rf <- ref()
    validate(need(input$m_range[1] < rf$m,
                  "The fraction range must include numbers below the comparison regimen."))
    m_hi <- min(input$m_range[2], rf$m - 1)
    early_mode <- input$early_mode
    late_mode <- choice_or("late_mode", "reference")
    cap <- switch(early_mode,
      fowler = fowler_zone(isTRUE(input$chemo))[1],
      hypno = hypno_caps()$early,
      NULL)
    late_cap <- if (late_mode == "hypno") hypno_caps()$late else NULL
    fenwick <- identical(early_mode, "fenwick")
    list(m_range = input$m_range[1]:m_hi, d_step = as.numeric(input$d_step),
         d_max = checked("d_max"), criterion = input$criterion, rank_by = input$rank_by,
         tol_early = if (is.null(cap) && !fenwick) checked("tol_early") else 0,
         tol_late = if (is.null(late_cap)) checked("tol_late") else 0,
         early_cap = cap, late_cap = late_cap,
         min_gain = checked("min_gain"), min_days_saved = checked("min_days"),
         fenwick = fenwick)
  })

  search <- reactive({
    if (input$early_mode == "fowler" && !fowler_defaults()) {
      s <- empty_suggestions()
      attr(s, "reason") <- "fowler_params_changed"
      return(list(s = s, cap = NULL, blocked = TRUE))
    }
    args <- search_args()
    s <- do.call(suggest_regimens, c(list(p = params(), ref = ref(), anchor = anchor()), args))
    list(s = s, cap = args$early_cap, late_cap = args$late_cap, blocked = FALSE)
  })
  suggestions <- reactive(search()$s)

  output$search_note <- renderUI({
    p <- params(); sr <- search()
    if (isTRUE(sr$blocked)) return(NULL)
    early_txt <- if (identical(input$early_mode, "fenwick")) "cell-kill BED (α/β 10 Gy) at or below the Fenwick et al. 2008 boundary for its overall time" else
      if (is.null(sr$cap)) "early BED at most the comparison's plus the allowed excess" else
      if (input$early_mode == "hypno") sprintf("early BED at most %.1f Gy₁₀ (HYPNO 55 Gy/20 fractions)", sr$cap) else
      sprintf("early BED at most %g Gy₁₀ (Fowler zone, lower edge%s)", sr$cap,
              if (isTRUE(input$chemo)) ", with concurrent chemotherapy" else "")
    late_txt <- if (is.null(sr$late_cap)) "late BED at most the comparison's plus the allowed excess" else
      sprintf("late BED at most %.1f Gy (HYPNO 55 Gy/20 fractions)", sr$late_cap)
    helpText(sprintf(paste("Searched: once-daily Mon-Fri regimens with a dose per fraction above the",
                           "comparison's (multiples of %s Gy, up to %g Gy), %s, and %s. For each number",
                           "of fractions, the highest-ranked dose per fraction is listed; the ZIP download",
                           "also holds every qualifying regimen. Δ = regimen minus comparison. Rows above",
                           "%g Gy per fraction are beyond the fraction sizes of the benchmark trials.",
                           "The list depends strongly on the tumor repopulation settings (now %.2f Gy/day under DI)."),
                     input$d_step, checked("d_max"), early_txt, late_txt, PRECEDENT_MAX_D,
                     p$lambda_di / p$alpha))
  })

  pretty_suggestions <- reactive({
    s <- suggestions()
    dd_col <- if (nrow(s) > 0 && all(is.na(s$d_dd))) rep("not available", nrow(s)) else s$d_dd
    data.frame(Fractions = s$m, `Gy per fraction` = s$d, Days = s$T, `Weeks saved` = s$weeks_saved,
               `Total Gy` = s$total_dose,
               `Δ BED DI` = s$d_di, `Δ BED DD` = dd_col, `Credited by` = s$credited_by,
               `Δ early BED` = s$d_early, `Δ late BED` = s$d_late,
               `Early EQD2` = s$eqd2_early, `Late EQD2` = s$eqd2_late,
               `Beyond precedent` = ifelse(s$beyond_precedent, "yes", ""),
               `Next dose step would exceed` = s$next_step, check.names = FALSE)
  })

  output$suggest_msg <- renderUI({
    s <- suggestions()
    if (nrow(s) > 0) return(NULL)
    reason <- attr(s, "reason")
    msg <- function(...) div(class = "alert alert-info", ...)
    if (identical(reason, "dd_unavailable"))
      return(msg("The DD model is unavailable with these parameters, so DD-based criteria and rankings",
                 "(including 'smaller of ΔDI and ΔDD') cannot be used."))
    if (identical(reason, "negative_reference_di") || identical(reason, "negative_reference_dd") ||
        identical(reason, "negative_reference_both"))
      return(msg("The comparison regimen's tumor BED is negative under the",
                 switch(reason, negative_reference_di = "DI model", negative_reference_dd = "DD model",
                        negative_reference_both = "DI and DD models"),
                 "with these parameters, so a negative model cannot credit any regimen.",
                 "Check the repopulation settings."))
    if (identical(reason, "negative_reference_early") || identical(reason, "negative_reference_late"))
      return(msg("The comparison regimen's toxicity BED is negative, so a limit relative to it is not meaningful.",
                 "Check the toxicity parameters or choose an absolute limit."))
    if (identical(reason, "empty_dose_range"))
      return(msg("No dose per fraction lies between the comparison dose per fraction and the highest dose searched."))
    if (identical(reason, "fowler_params_changed"))
      return(msg("The Fowler tolerance zone applies only with Fowler's acute parameters",
                 "(the defaults: α/β 10 Gy, α 0.35/Gy, doubling time 2.5 d, onset day 7 elapsed).",
                 "Restore them or use another early limit."))
    if (identical(reason, "no_shorter_regimen"))
      return(msg("No number of fractions in the searched range saves the minimum number of days."))
    msg("No regimen in the searched range meets all limits with these settings.",
        "This is a result, not an error: with these parameters, shortening the course costs",
        "tumor BED faster than the limits allow.",
        if (input$criterion != "either") "Options: require the change under 'DI or DD model'," else "Options:",
        "allow a tumor BED loss (a negative minimum change), lower the days saved,",
        "or choose another early or late limit.")
  })

  output$suggest_table <- renderDT({
    s <- pretty_suggestions()
    two_dp <- c("Gy per fraction", "Total Gy", "Δ BED DI", "Δ BED DD", "Δ early BED",
                "Δ late BED", "Early EQD2", "Late EQD2")
    two_dp <- two_dp[vapply(s[two_dp], is.numeric, logical(1))]
    datatable(s, rownames = FALSE, fillContainer = FALSE,
              options = list(pageLength = 15, dom = "tip", scrollX = TRUE)) |>
      formatRound(two_dp, 2) |>
      formatRound("Weeks saved", 1)
  })

  scenario <- reactive({
    sr <- search()
    validate(need(!isTRUE(sr$blocked), "Nothing to export with the current settings."))
    build_scenario(p = params(), ref = ref(), anchor = anchor(), search = search_args(),
                   context = list(comparator = choice_or("comparator", "custom"),
                                  preset = input$preset, early_mode = input$early_mode,
                                  late_mode = choice_or("late_mode", "reference"),
                                  chemo = isTRUE(input$chemo), coupled = coupled(),
                                  regimen_evaluated = regimen()))
  })

  output$download <- downloadHandler(
    filename = function() sprintf("bed_v2.1_%s_%s.zip", scenario_id(scenario()), format(Sys.Date())),
    content = function(file) {
      s <- suggestions()
      write_export_zip(file, scenario(), s, attr(s, "feasible"))
    },
    contentType = "application/zip")

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
           y = "Tumor BED minus comparison (Gy)",
           title = "Tumor BED difference from the comparison for the listed regimens") +
      theme_minimal(base_size = 15) +
      theme(axis.title = element_text(face = "bold"), legend.position = "top")
  })

  sens <- reactive(sensitivity_table(regimen(), ref(), anchor(), params()))

  output$sens_summary <- renderUI({
    st <- sens(); mg <- checked("min_gain")
    n <- nrow(st)
    di_n <- sum(st$d_di >= mg - 1e-9)
    dd_ok <- !is.na(st$d_dd)
    dd_n <- sum(st$d_dd[dd_ok] >= mg - 1e-9)
    tags$p(sprintf(paste("Tumor BED change at least %g Gy (the minimum set on the candidate tab):",
                         "%d of %d scenarios under DI, %d of %d calibratable scenarios under DD."),
                   mg, di_n, n, dd_n, sum(dd_ok)))
  })

  output$sens_table <- renderTable({
    st <- sens()
    data.frame(Scenario = unname(PRESET_SHORT[st$preset]), `Tumor α/β (Gy)` = sprintf("%g", st$ab),
               `α (1/Gy)` = sprintf("%.3g", st$alpha), `DI loss (Gy/day)` = sprintf("%.2f", st$lambda_di / st$alpha),
               `DI onset (day, inclusive)` = sprintf("%g", st$Tk), `DD max rate (1/day)` = sprintf("%.3g", st$lambda_dd),
               `Δ BED DI` = fmt(st$d_di, TRUE), `Δ BED DD` = fmt(st$d_dd, TRUE), check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)

  output$sens_tox <- renderTable({
    r <- regimen(); rf <- ref(); p <- params(); b <- both()
    data.frame(Quantity = c("Early BED (current acute parameters)", "Late BED, α/β 2 Gy",
                            "Late BED, α/β 3 Gy", "Late BED, α/β 4 Gy"),
               `Δ (regimen − comparison)` = fmt(c(b$bed_early[1] - b$bed_early[2],
                 bed_late(r$m, r$d, 2) - bed_late(rf$m, rf$d, 2),
                 bed_late(r$m, r$d, 3) - bed_late(rf$m, rf$d, 3),
                 bed_late(r$m, r$d, 4) - bed_late(rf$m, rf$d, 4)), TRUE), check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)

  output$about <- renderUI(withMathJax(HTML(
    markdown::mark(file = "about.md", options = "+latex_math"))))
}

shinyApp(ui, server)
