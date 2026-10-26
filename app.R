library(shiny)
library(bslib)
library(ggplot2)
library(DT)

source("R/bed_models.R", local = TRUE)
source("R/schedule.R", local = TRUE)
source("R/uncertainty.R", local = TRUE)
source("R/oar.R", local = TRUE)
source("R/export.R", local = TRUE)
source("R/mod_oar.R", local = TRUE)
source("R/mod_uncertainty.R", local = TRUE)
# Identity of the code as loaded, so an export made after the files on disk
# have changed says so (code$changed_since_load).
CODE_AT_LOAD <- code_identity()

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
  d_max = list("Highest dose per fraction searched", 2, 30, FALSE),
  agree_margin = list("Own threshold for the both-model screen", -20, 30, FALSE),
  alpha_dd = list("DD α (fitted)", 0.001, 2, FALSE),
  c_over_alpha_dd = list("DD onset constant C/α (fitted)", 0.01, 1e5, FALSE),
  g_dd = list("DD background growth (fitted)", 0, 0.1, FALSE),
  cmp1_m = list("Regimen A fractions", 1, 100, TRUE), cmp1_d = list("Regimen A dose per fraction", 0.1, 30, FALSE),
  cmp2_m = list("Regimen B fractions", 1, 100, TRUE), cmp2_d = list("Regimen B dose per fraction", 0.1, 30, FALSE),
  cmp3_m = list("Regimen C fractions", 1, 100, TRUE), cmp3_d = list("Regimen C dose per fraction", 0.1, 30, FALSE),
  agree_days = list("Minimum days saved (agreement)", 1, 140, TRUE))

PRESET_LABELS <- c(
  "Version 1 values (defaults of the published version 1 website)" = "calculator",
  "Withers et al. 1988 rule: 0.6 Gy/day in 2-Gy fractions after 4 weeks (α 0.2 chosen, follows α/β)" = "withers",
  "Hartley et al. 2011: α 0.3/Gy, doubling time 3 d, onset day 22 (elapsed)" = "hartley",
  "HYHOPE protocol: α 0.35/Gy, doubling time 3 d, onset day 21" = "hyhope",
  "Scenario D18 (default): α 0.35/Gy, doubling time 3 d, onset day 28" = "d18",
  "Verified-trial refit 2026-10-25: posterior medians (DI; DD calibrated)" = "refit",
  "Custom (edited values)" = "custom")

# Point values of the verified-trial refit (fit/export_draws.py): medians of the
# identified DI quantities, and the calibrated-DD rate median (DD is calibrated to
# DI at the anchor, as in the sidebar's default DD mode). Offered as a preset;
# not part of the sensitivity set.
REFIT <- jsonlite::fromJSON(file.path(DRAWS_DIR, "refit_point_values.json"))$values
refit_preset <- function() list(alpha = REFIT$di$alpha, ab = REFIT$di$ab, lambda_di = REFIT$di$lambda,
                                Tk = REFIT$di$Tk, g = REFIT$di$g, lambda_dd = REFIT$dd$lambda)
refit_dd_values <- function() list(alpha_dd = REFIT$dd$alpha, lambda_dd = REFIT$dd$lambda,
                                   c_over_alpha_dd = REFIT$dd$c_over_alpha, g_dd = REFIT$dd$g)

WEEKDAYS_IN <- c("Monday" = 1, "Tuesday" = 2, "Wednesday" = 3, "Thursday" = 4, "Friday" = 5, "Saturday" = 6)

# Short names for the sensitivity table.
PRESET_SHORT <- c(calculator = "Version 1 defaults", withers = "Withers 1988 (rule)",
                  hartley = "Hartley 2011", hyhope = "HYHOPE protocol", d18 = "Scenario D18")

COMPARATOR_LABELS <- c(
  "Enter below (default 35 × 2 Gy)" = "custom",
  "HYPNO test arm: 55 Gy in 20 fractions (26 days)" = "hypno_test",
  "HYPNO control arm: 66 Gy in 33 fractions, 6 fractions/week (38 days)" = "hypno_control")

DEF <- default_params()
# A difference smaller than this is reported as "about the same" in words
# (plan workstream C, item 7); tables keep full values.
VERDICT_TOL <- 0.5
# Benchmark fraction-size flag, derived from the constant (not typed).
PRECEDENT_HEADER <- sprintf("Above benchmark fraction sizes (> %g Gy)", PRECEDENT_MAX_D)

# Download file name: version, analysis, scenario ID and date.
export_filename <- function(analysis, id, date = Sys.Date()) {
  sprintf("bed_v%s_%s_%s_%s.zip", APP_VERSION, analysis, id, format(date))
}

# BED unit for a tissue alpha/beta: Gy with a subscript for the conventional
# 10 Gy and 3 Gy, otherwise the alpha/beta spelled out.
bed_unit <- function(ab) {
  if (isTRUE(all.equal(ab, 10))) "Gy\u2081\u2080" else if (isTRUE(all.equal(ab, 3))) "Gy\u2083" else
    sprintf("Gy (\u03b1/\u03b2 = %g Gy)", ab)
}

# Explanation of an empty search, shared by the candidate and agreement tabs.
reason_text <- function(reason, criterion = "either") {
  if (is.null(reason)) return(NULL)
  txt <- switch(reason,
    dd_unavailable = c("The DD model is unavailable with these parameters, so DD-based criteria and rankings",
                       "(including 'smaller of \u0394DI and \u0394DD' and the two-model agreement) cannot be used."),
    negative_reference_di = , negative_reference_dd = , negative_reference_both =
      c("The comparison regimen's tumor BED is not positive under the",
        switch(reason, negative_reference_di = "DI model", negative_reference_dd = "DD model",
               negative_reference_both = "DI and DD models"),
        "with these parameters, so that model cannot credit any regimen. Check the repopulation settings."),
    negative_reference_early = , negative_reference_late =
      c("The comparison regimen's toxicity BED is not positive, so a limit relative to it is not meaningful.",
        "Check the toxicity parameters or choose an absolute limit."),
    empty_dose_range = "No dose per fraction lies between the comparison dose per fraction and the highest dose searched.",
    fowler_params_changed = c("The Fowler tolerance zone applies only with Fowler's acute parameters",
                              "(the defaults: \u03b1/\u03b2 10 Gy, \u03b1 0.35/Gy, doubling time 2.5 d, onset day 7 elapsed).",
                              "Restore them or use another early limit."),
    no_shorter_regimen = "No number of fractions in the searched range saves the minimum number of days.",
    no_rankable_row = c("Some regimens meet every limit, but none can be ranked as chosen: the ranking needs a",
                        "tumor model whose difference is not evaluable for them (a tumor BED is not positive, or DD",
                        "is unavailable). They are listed in the ZIP download (all_qualifying.csv, rank_eligible",
                        "FALSE). Rank by the larger change or by the other model, or check the repopulation settings."),
    paste("No search was possible:", reason))
  paste(txt, collapse = " ")
}

# Short list of what the model leaves out, shown beside the main results.
# Defined at top level because the UI uses it.
caveat_strip <- function(id) {
  div(id = id, class = "bed-caveat alert alert-light", style = "padding: 0.5rem 0.75rem;",
      tags$strong("Model limits: "),
      "Chemotherapy effects, HPV status, stage, tumor volume and hypoxia are not modeled; the tumor parameters come",
      "mostly from pre-HPV-era trials without chemotherapy. Overall duration is included, but individual gap timing",
      "is not modeled separately. Organ-at-risk doses count only where entered on the Organs at risk tab.",
      "Parameter uncertainty is shown on the Uncertainty tab. Default tumor parameters are exploratory.",
      "See Model and sources.")
}

# Five-line orientation card (plan workstream C, item 12).
how_to_read <- function() {
  div(class = "alert alert-info", style = "padding: 0.6rem 0.9rem;",
      tags$strong("How to read this. "),
      tags$ol(style = "margin-bottom:0;",
        tags$li("Every difference (Δ) is the evaluated regimen minus the comparison regimen named in the banner above."),
        tags$li("Tumor BED is shown under two repopulation models, DI and DD; they can disagree, and neither is averaged."),
        tags$li("Toxicity limits are model quantities (acute and late BED, organ EQD2), with their values shown next to them."),
        tags$li("'Within' and 'exceeds' refer to the stated limit under the stated parameters, never to safety."),
        tags$li("Outputs generate hypotheses for a trial; they are not treatment recommendations.")))
}

ui <- page_sidebar(
  fillable = FALSE,
  title = paste0("Biologically effective dose (BED) calculator for head and neck squamous cell carcinoma (HNSCC), version ",
                 sub("\\.0$", "", APP_VERSION)),
  theme = bs_theme(version = 5, bootswatch = "flatly") |>
    bs_add_rules(paste(
      ".bslib-page-title, .navbar-brand { white-space: normal; overflow-wrap: anywhere; font-size: clamp(1rem, 2.5vw, 1.25rem); }",
      "#about .MathJax_Display, #about mjx-container[display], #about .MathJax_SVG_Display { overflow-x: auto; overflow-y: hidden; max-width: 100%; }",
      "#about table { display: block; overflow-x: auto; max-width: 100%; }",
      ".bed-banner { position: sticky; top: 0; z-index: 10; }",
      ".bed-override { border-left: 4px solid #E69F00; padding-left: 0.5rem; }")),
  sidebar = sidebar(
    width = 380, open = list(desktop = "open", mobile = "closed"),
    accordion(
      open = c("Regimen to evaluate", "Comparison and tumor threshold"),
      accordion_panel("Regimen to evaluate",
        num("m", "Number of fractions", 18, 1, 100, 1),
        num("d", "Dose per fraction (Gy)", 2.9, 0.1, 30, 0.05),
        helpText("Default: 18 × 2.9 Gy (52.2 Gy), an example schedule chosen by the developer. It is not a published or",
                 "trial-tested regimen and not a recommendation. HYPNO's tested 55 Gy in 20 fractions is a benchmark on the Results tab."),
        radioButtons("fx_week", "Fractions per week", c("5 (Mon-Fri)" = 5, "6 (Mon-Sat)" = 6), inline = TRUE),
        selectInput("start_day", "First treatment day", WEEKDAYS_IN, selected = 1),
        textInput("skip_days", "Days without treatment (holidays, gaps): calendar day numbers, day 1 = first treatment", ""),
        helpText("Skipped days apply to the evaluated regimen only. Both searches build candidate schedules on the",
                 "regular calendar (fractions per week and first day), without skipped days."),
        checkboxInput("auto_T", "Compute overall time from this calendar", TRUE),
        conditionalPanel("!input.auto_T",
          num("T", "Overall treatment time (days)", 24, 1, 150, 1)),
        uiOutput("T_shown")),
      accordion_panel("Comparison and tumor threshold",
        selectInput("comparator", "Compare with (used by every tab)", COMPARATOR_LABELS),
        conditionalPanel("input.comparator == 'custom'",
          num("ref_m", "Number of fractions", 35, 2, 100, 1),
          num("ref_d", "Dose per fraction (Gy)", 2, 0.1, 30, 0.05),
          radioButtons("ref_fx_week", "Fractions per week", c("5" = 5, "6" = 6), inline = TRUE),
          checkboxInput("ref_auto_T", "Compute overall time (starting Monday)", TRUE),
          conditionalPanel("!input.ref_auto_T",
            num("ref_T", "Overall treatment time (days)", 47, 1, 150, 1))),
        uiOutput("ref_T_shown"),
        num("min_gain", "Tumor threshold: minimum tumor BED change from the comparison (Gy; negative = allowed loss)", 1, -20, 30, 0.5),
        helpText("One comparison and one threshold for every tab. The both-model screen can override them explicitly;",
                 "an override is highlighted. The HYPNO, Fowler and Fenwick toxicity limits are fixed values.",
                 "The comparison does not change the DD model, which is calibrated at the anchor (or fitted).")),
      accordion_panel("Tumor (DI and DD models)",
        selectInput("preset", "Parameter preset", PRESET_LABELS, selected = DEFAULT_PRESET),
        uiOutput("preset_note"),
        num("alpha", "α, tumor cell killing (1/Gy)", DEF$alpha, 0.01, 2, 0.01),
        num("ab", "α/β, tumor (Gy)", DEF$ab, 0.5, 100, 1),
        num("lambda_di", "Repopulation rate, DI model (1/day)", DEF$lambda_di, 0, 2, 0.001),
        num("lambda_dd", "Maximum repopulation rate, DD model (1/day)", DEF$lambda_dd, 0, 5, 0.01),
        num("Tk", "Onset of accelerated repopulation, DI model (days)", DEF$Tk, 0, 100, 1),
        num("g", "Background tumor growth rate g (1/day, 0 = off)", DEF$g, 0, 0.1, 0.001),
        radioButtons("dd_mode", "DD model",
                     c("Calibrated to DI at the anchor (versions 1-2; default)" = "calibrated",
                       "Refit point values, entered directly (no calibration)" = "fitted")),
        conditionalPanel("input.dd_mode == 'fitted'",
          num("alpha_dd", "DD α (1/Gy)", REFIT$dd$alpha, 0.001, 2, 0.001),
          num("c_over_alpha_dd", "DD onset constant C/α (Gy)", REFIT$dd$c_over_alpha, 0.01, 1e5, 1),
          num("g_dd", "DD background growth (1/day)", REFIT$dd$g, 0, 0.1, 0.001),
          helpText("This mode uses its own α, onset constant and growth, and the maximum rate above; no calibration.",
                   "The values shown are the refit's: its DD rate median with the onset constant that makes DD equal DI",
                   "at 35 × 2 Gy in 47 days under the refit's DI medians. Changing the DI values does not change them.")),
        uiOutput("tumor_note")),
      accordion_panel("DD calibration anchor",
        num("anchor_m", "Anchor number of fractions", 35, 2, 100, 1),
        num("anchor_d", "Anchor dose per fraction (Gy)", 2, 0.1, 30, 0.05),
        uiOutput("anchor_note")),
      accordion_panel("Early toxicity (acute mucosa)",
        num("ab_e", "α/β, mucosa (Gy)", 10, 0.5, 100, 1),
        num("alpha_e", "α, mucosa (1/Gy)", 0.35, 0.01, 2, 0.01),
        num("Tp_e", "Mucosal doubling time during repopulation (days)", 2.5, 0.1, 50, 0.1),
        num("Tk_e", "Onset of mucosal repopulation (days, elapsed)", 7, 0, 100, 1),
        uiOutput("lambda_e_shown"),
        checkboxInput("chemo", "Concurrent chemotherapy: lower the Fowler zone only (no BED changes)", FALSE),
        helpText("This box has one effect: it lowers the Fowler acute zone by 3-5 Gy₁₀ where that zone is shown or used",
                 "as the early limit. It changes nothing under the default (HYPNO) early limit. No BED includes a",
                 "chemotherapy effect."),
        uiOutput("cisplatin_context")),
      accordion_panel("Late toxicity",
        num("ab_l", "α/β, late-responding tissue (Gy)", 3, 0.5, 20, 0.5),
        helpText("A coarse screen at the prescription dose: late BED = n d (1 + d/(α/β)), no repopulation.",
                 "Organ doses are assessed on the Organs at risk tab."))
    )
  ),
  div(class = "bed-banner", uiOutput("banner")),
  navset_card_tab(
    nav_panel("Report card",
      how_to_read(),
      uiOutput("dd_warning0"),
      div(style = "overflow-x:auto;", tableOutput("report")),
      uiOutput("report_notes"),
      caveat_strip("caveat-report")),
    nav_panel("Results",
      uiOutput("dd_warning"),
      uiOutput("neg_warning"),
      helpText("Tumor BED is shown under two repopulation models: dose-independent (DI),",
               "with a fixed onset and rate, and dose-dependent (DD), where onset and rate",
               "depend on cell killing. EQD2 is the equivalent dose in 2-Gy fractions.",
               "Model outputs for hypothesis generation, not treatment recommendations."),
      div(class = "table-scroll", style = "overflow-x: auto;", tableOutput("results")),
      caveat_strip("caveat-results"),
      uiOutput("verdicts"),
      uiOutput("acute_zone"),
      h5("Published regimens, recomputed with the current parameters"),
      div(class = "table-scroll", style = "overflow-x: auto;", tableOutput("benchmarks")),
      uiOutput("hypno_note")),
    nav_panel("Compare regimens",
      helpText("Up to three further regimens side by side with the evaluated regimen and the published benchmarks,",
               "all on 5 fractions/week from a Monday unless they are benchmarks with their own calendar."),
      layout_columns(col_widths = c(4, 4, 4),
        div(num("cmp1_m", "Regimen A: fractions", 17, 1, 100, 1), num("cmp1_d", "Gy per fraction", 3.05, 0.1, 30, 0.05)),
        div(num("cmp2_m", "Regimen B: fractions", 15, 1, 100, 1), num("cmp2_d", "Gy per fraction", 3.1, 0.1, 30, 0.05)),
        div(num("cmp3_m", "Regimen C: fractions", 30, 1, 100, 1), num("cmp3_d", "Gy per fraction", 2, 0.1, 30, 0.05))),
      div(style = "overflow-x:auto;", tableOutput("compare_table"))),
    nav_panel("Search",
      uiOutput("dd_warning2"),
      div(class = "alert alert-secondary",
          "These are model outputs for hypothesis generation, not treatment recommendations.",
          "A BED difference is not a difference in tumor control. A listed regimen is not shown to be safe or effective.",
          "Any trial regimen needs clinical precedent, organ-at-risk planning and toxicity stopping rules."),
      caveat_strip("caveat-candidates"),
      layout_columns(
        col_widths = c(4, 4, 4),
        selectInput("criterion", "Tumor threshold applied under",
                    c("both models (the both-model screen)" = "both", "DI or DD model" = "either",
                      "DI model" = "di", "DD model" = "dd"), selected = "both"),
        selectInput("rank_by", "Rank by",
                    c("larger of ΔDI and ΔDD (optimistic)" = "best",
                      "smaller of ΔDI and ΔDD (conservative)" = "worst",
                      "ΔBED DI" = "di", "ΔBED DD" = "dd"), selected = "worst"),
        sliderInput("m_range", "Number of fractions searched", 5, 34, c(5, 34), step = 1),
        num("min_days", "Minimum days saved", 7, 1, 140, 1),
        selectInput("d_step", "Dose-per-fraction step (Gy)", c(0.05, 0.1, 0.2, 0.25, 0.5), selected = 0.05),
        num("d_max", "Highest dose per fraction searched (Gy)", 5, 2, 30, 0.5),
        div(radioButtons("early_mode", "Early toxicity limit",
                         c("Comparison regimen, plus allowed excess" = "reference",
                           "Lower edge of the Fowler tolerance zone" = "fowler",
                           "Acute BED of HYPNO 55 Gy/20 fractions" = "hypno",
                           "Fenwick et al. 2008 boundary (cell-kill BED, no repopulation; up to 50 days)" = "fenwick"),
                         selected = DEFAULT_EARLY_MODE),
            uiOutput("early_limit_value")),
        conditionalPanel("input.early_mode == 'reference'",
          num("tol_early", "Allowed early BED above the comparison (Gy)", 0, 0, 30, 0.5)),
        div(radioButtons("late_mode", "Late toxicity limit (coarse screen at the prescription dose)",
                         c("Comparison regimen, plus allowed excess" = "reference",
                           "Late BED of HYPNO 55 Gy/20 fractions" = "hypno")),
            uiOutput("late_limit_value")),
        conditionalPanel("input.late_mode == 'reference'",
          num("tol_late", "Allowed late BED above the comparison (Gy)", 0, 0, 30, 0.5))),
      uiOutput("search_context"),
      uiOutput("search_note"),
      uiOutput("hypno_cap_note"),
      uiOutput("suggest_msg"),
      DTOutput("suggest_table", fill = FALSE),
      br(),
      downloadButton("download", "Download results and scenario (ZIP)"),
      plotOutput("suggest_plot", height = "380px"),
      helpText("Each point is the listed regimen for that number of fractions, so the dose per fraction",
               "changes from point to point. The plot is not a dose-response curve, and a model's point is",
               "missing where its difference is not evaluable.")),
    nav_panel("Both-model screen",
      helpText("Regimens whose tumor BED change from the comparison is at least the threshold under DI and under DD,",
               "with the toxicity limits chosen on the Search tab. It searches every number of fractions below the",
               "comparison. The ranges and counts use 15 selected parameter settings: five tumor presets × three tumor",
               "α/β values. Only settings in which both model differences are evaluable enter the joint count and ranges.",
               "These are sensitivity checks, not confidence intervals or probabilities (see the Uncertainty tab)."),
      layout_columns(
        col_widths = c(4, 4, 4),
        selectInput("agree_cmp", "Comparison",
                    c("The global comparison (sidebar)" = "sidebar",
                      "Override: HYPNO test arm, 55 Gy in 20 fractions" = "hypno_test")),
        div(radioButtons("agree_threshold_mode", "Threshold",
                         c("The global tumor threshold (sidebar)" = "global", "Override for this tab" = "own")),
            conditionalPanel("input.agree_threshold_mode == 'own'",
              num("agree_margin", "Own threshold (Gy; negative = allowed loss)", -3, -20, 30, 0.5))),
        num("agree_days", "Minimum days saved", 1, 1, 140, 1)),
      uiOutput("agree_override"),
      caveat_strip("caveat-agreement"),
      uiOutput("agree_context"),
      uiOutput("agree_note"),
      uiOutput("agree_msg"),
      div(class = "table-scroll", style = "overflow-x: auto;", tableOutput("agree_table")),
      downloadButton("download_agree", "Download the both-model screen (ZIP)")),
    nav_panel("Sensitivity",
      helpText("Tumor BED of the regimen being evaluated minus the comparison regimen, under",
               "each tumor preset and three tumor α/β values. These are selected parameter settings,",
               "not a probability distribution: a count of settings is not a chance of benefit.",
               "Presets that set only DI values keep the DD maximum rate entered in the sidebar, and every",
               "setting keeps the sidebar's background growth g and DD anchor. A difference is shown only",
               "when both BEDs are positive (and the DD model is available)."),
      uiOutput("sens_context"),
      uiOutput("sens_summary"),
      div(class = "table-scroll", style = "overflow-x: auto;", tableOutput("sens_table")),
      h5("Toxicity BED differences (do not depend on the tumor setting)"),
      tableOutput("sens_tox")),
    nav_panel("Uncertainty", unc_ui("unc")),
    nav_panel("Organs at risk", oar_ui("oar")),
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
    if (input$preset == "refit") {
      pr <- refit_preset()
      for (k in names(pr)) updateNumericInput(session, k, value = pr[[k]])
      return()
    }
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
    rp <- refit_preset()
    if (isTRUE(all.equal(c(input$alpha, input$lambda_di, input$Tk, ab), c(rp$alpha, rp$lambda_di, rp$Tk, rp$ab), tolerance = 1e-6)))
      now <- "refit"
    if (!identical(now, input$preset)) updateSelectInput(session, "preset", selected = now)
  }, ignoreInit = TRUE)
  observeEvent(input$dd_mode, {
    if (identical(input$dd_mode, "fitted")) {
      v <- refit_dd_values()
      for (k in names(v)) updateNumericInput(session, k, value = v[[k]])
    }
  }, ignoreInit = TRUE)

  params <- reactive({
    list(alpha = checked("alpha"), ab = checked("ab"),
         lambda_di = checked("lambda_di"), lambda_dd = checked("lambda_dd"),
         Tk = checked("Tk"), g = checked_or("g", 0),
         ab_e = checked("ab_e"), alpha_e = checked("alpha_e"),
         lambda_e = log(2) / checked("Tp_e"), Tk_e = checked("Tk_e"),
         clock_e = "elapsed", ab_l = checked("ab_l"),
         dd_mode = choice_or("dd_mode", "calibrated"),
         alpha_dd = if (identical(input$dd_mode, "fitted")) checked("alpha_dd") else NULL,
         c_over_alpha_dd = if (identical(input$dd_mode, "fitted")) checked("c_over_alpha_dd") else NULL,
         g_dd = if (identical(input$dd_mode, "fitted")) checked("g_dd") else NULL)
  })

  # Overall time: from the calendar, or typed but not shorter than once-daily
  # Mon-Fri treatment allows.
  calendar <- reactive({
    fpw <- as.integer(choice_or("fx_week", 5))
    sd <- as.integer(choice_or("start_day", 1))
    validate(need(sd <= fpw, "With 5 fractions per week the first treatment day must be Monday to Friday."))
    skip <- tryCatch(parse_skip_days(input$skip_days), error = function(e) validate(need(FALSE, conditionMessage(e))))
    list(fx_per_week = fpw, start_weekday = sd, skip_days = skip)
  })
  time_for <- function(m, auto, id, cal = default_calendar()) {
    shortest <- schedule_time(m, cal$fx_per_week, cal$start_weekday, cal$skip_days)
    if (isTRUE(auto)) return(shortest)
    T <- checked(id)
    validate(need(T >= shortest, sprintf(
      "%d fractions on this calendar (%s) need at least %g days. %g days was entered. Twice-daily treatment is not covered.",
      m, calendar_text(cal), shortest, T)))
    T
  }
  ref <- reactive({
    cmp <- choice_or("comparator", "custom")
    if (cmp != "custom") {
      r <- comparator_presets()[[cmp]]
      return(list(m = r$m, d = r$d, T = r$T, label = r$label))
    }
    m <- checked("ref_m")
    cal <- list(fx_per_week = as.integer(choice_or("ref_fx_week", 5)), start_weekday = 1, skip_days = integer())
    list(m = m, d = checked("ref_d"), T = time_for(m, input$ref_auto_T, "ref_T", cal))
  })
  regimen <- reactive({
    m <- checked("m")
    list(m = m, d = checked("d"), T = time_for(m, input$auto_T, "T", calendar()))
  })
  fixed_anchor <- reactive({
    m <- checked_or("anchor_m", 35)
    list(m = m, d = checked_or("anchor_d", 2), T = overall_time(m))
  })
  # The DD anchor: the same for every tab, whatever the comparison regimen.
  anchor <- reactive(fixed_anchor())
  calibration <- reactive({
    dd_setup(params(), anchor())
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
    helpText(sprintf("Used: %d fractions, %g days (%s), total dose %.2f Gy.", r$m, r$T, calendar_text(calendar()), r$m * r$d))
  })
  output$ref_T_shown <- renderUI({
    r <- ref()
    helpText(sprintf("Used: %d fractions, %g days, total dose %.2f Gy.", r$m, r$T, r$m * r$d))
  })
  output$lambda_e_shown <- renderUI(
    helpText(sprintf(paste("Mucosal repopulation rate = ln2 / doubling time = %.3f per day.",
                           "Overall time is converted to Fowler's elapsed days (calendar days minus 1)."),
                     params()$lambda_e)))
  output$preset_note <- renderUI({
    keeps_dd <- "It keeps the DD maximum rate entered below."
    txt <- switch(input$preset,
      d18 = paste("The default. Scenario D18 was chosen by the developer on 2026-10-18 and is not fitted to",
                  "outcome data. Its DI values are within published ranges: repopulation at 0.23 per day",
                  "(doubling time 3.0 days) from day 28, a loss of 0.66 Gy₁₀ of tumor BED per day, at the low",
                  "end of published estimates. Its DD maximum rate (0.40 per day) is a constructed value, not",
                  "a published DD estimate. See \"Model and sources\"."),
      calculator = paste("The tumor values of the published version 1 website: α 0.2/Gy, DI repopulation 0.2 per",
                         "day from day 28 (a loss of 1.0 Gy₁₀ per day) and DD maximum rate 0.5 per day."),
      withers = paste("Withers et al. 1988 supply the rule: repopulation compensated by about 0.6 Gy per day, in",
                      "2-Gy fractions, after about 4 weeks (onset day 28 here). α 0.2/Gy is this calculator's",
                      "choice to complete the preset (the version 1 value). The DI loss per day in BED does not",
                      "depend on it. The DI rate follows the tumor α/β.", keeps_dd),
      hartley = paste("Hartley et al. 2011: α 0.3/Gy, α/β 10 Gy, doubling time 3 days, onset on elapsed day 22",
                      "(day 23 in this calculator's count).", keeps_dd),
      hyhope = paste("HYHOPE protocol: α 0.35/Gy, α/β 10 Gy, doubling time 3 days, onset day 21.", keeps_dd),
      refit = sprintf(paste("Medians of the Bayesian refit to randomized fractionation trials (mainly without chemotherapy), re-extracted from",
                            "the trial publications (2026-10-25; one intercept per trial, overdispersion): α %.3f/Gy,",
                            "repopulation loss %.2f Gy₁₀/day after day %.0f, plus background growth %.2f Gy₁₀/day over the",
                            "whole course, and DD maximum rate %.2f per day (DD calibrated to DI at the anchor). α/β is not",
                            "informed by these data and stays 10 Gy. A small aggregate α may reflect variation between",
                            "patients, mixed endpoints and the model's form. Not the default."),
                      REFIT$di$alpha, REFIT$di$loss_per_day_after_onset, REFIT$di$Tk, REFIT$di$growth_loss_per_day,
                      REFIT$dd$lambda),
      custom = "The tumor values below were edited and match no preset.",
      "")
    helpText(txt, paste("A preset fills only the values it sets (named above). The menu shows Custom when",
                        "one of those values is edited. Values a preset does not set (the DD maximum rate for",
                        "the Withers, Hartley and HYHOPE presets, and g for all presets except the refit) keep what is entered."))
  })
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
    intro <- sprintf(paste("The DD model has one constant that is not entered: it is set so that DD and DI give",
                           "the same tumor BED for this anchor course (%s). The two models therefore agree",
                           "for the anchor course and differ only for other schedules."), where)
    txt <- switch(cal$status,
      ok = sprintf("With the current tumor values, DD accelerated repopulation in the anchor course starts on day %.1f (DI: day %g).",
                   cal$onset_ref, params()$Tk),
      no_repopulation = "The DI repopulation rate is 0, so DD has no accelerated repopulation either.",
      infeasible = "With the current tumor values no DD calibration exists at this anchor, so DD results are not shown.",
      fitted = sprintf("The DD model now uses fitted values instead, so the anchor is not used for DD. Its onset in the anchor course would be day %.1f.",
                       cal$onset_ref))
    helpText(intro, txt,
             "The anchor is separate from the comparison regimen: changing the comparison changes the differences shown, not the DD model.",
             "Version 1 used a fixed 35 × 2 Gy, 47-day anchor (with an α/β error, see Model and sources). Version 2.0",
             "calibrated at the comparison regimen. To approximate it, enter the comparison's fractions and dose here.",
             "The anchor's time always comes from the once-daily Mon-Fri calendar, so a comparison with another",
             "overall time (such as HYPNO's 33 fractions in 38 days at 6 per week, 45 days here) cannot be reproduced exactly.")
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
  # Tumor BED difference under the shared validity rule (tumor_endpoint):
  # shown only when both BEDs are positive.
  # y is one reference BED for every x, or one per x; any other length is an error (review R04:
  # a scalar y indexed by row was NA after row 1, so benchmark rows 2-3 read "not available").
  tumor_delta_txt <- function(x, y) {
    if (length(y) == 1) y <- rep(y, length(x))
    if (length(y) != length(x)) stop("tumor_delta_txt: reference must have length 1 or length(x)")
    vapply(seq_along(x), function(i) {
      if (is.na(x[i]) || is.na(y[i])) return("not available")
      e <- tumor_endpoint(x[i], y[i])
      if (isTRUE(e$valid)) sprintf("%+.2f", e$delta) else "not evaluable (BED \u2264 0)"
    }, character(1))
  }

  output$results <- renderTable({
    b <- both(); r <- regimen(); rf <- ref()
    pp <- params()
    rows <- c("T", "total_dose", "bed_di", "bed_dd", "bed_early", "eqd2_early", "bed_late", "eqd2_late")
    names(rows) <- c("Overall time (days)", "Total dose (Gy)",
                     sprintf("Tumor BED, DI model (%s)", bed_unit(pp$ab)),
                     sprintf("Tumor BED, DD model (%s)", bed_unit(pp$ab)),
                     sprintf("Early toxicity BED (%s)", bed_unit(pp$ab_e)), "Early toxicity EQD2 (Gy)",
                     sprintf("Late toxicity BED (%s)", bed_unit(pp$ab_l)), "Late toxicity EQD2 (Gy)")
    v1 <- unlist(b[1, rows]); v2 <- unlist(b[2, rows])
    diff <- ifelse(rows %in% c("bed_di", "bed_dd"), tumor_delta_txt(v1, v2), fmt(v1 - v2, plus = TRUE))
    out <- data.frame(Quantity = names(rows), `This regimen` = fmt(v1), Comparison = fmt(v2),
                      Difference = diff, check.names = FALSE)
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
      if (is.na(x) || is.na(y)) "not available" else if (x <= 0 || y <= 0) "not evaluable (a BED is not positive)" else
        if (x > y + VERDICT_TOL) "higher" else if (x < y - VERDICT_TOL) "lower" else
          sprintf("about the same (within %g Gy)", VERDICT_TOL)
    }
    tags$p(class = "text-muted",
      sprintf(paste("Compared with the comparison regimen (model values): tumor BED is %s under DI",
                    "and %s under DD. Early toxicity BED is %s, and late toxicity BED is %s."),
              word(b$bed_di[1], b$bed_di[2]), word(b$bed_dd[1], b$bed_dd[2]),
              word(b$bed_early[1], b$bed_early[2]), word(b$bed_late[1], b$bed_late[2])))
  })

  # Unit label for early BED: Gy10 only when the mucosal alpha/beta is 10 Gy.
  early_unit <- reactive(bed_unit(params()$ab_e))

  output$neg_warning <- renderUI({
    b <- both()
    v <- unlist(b[, c("bed_di", "bed_dd", "bed_early")])
    if (any(v < 0, na.rm = TRUE))
      div(class = "alert alert-warning",
          "Some BEDs are negative: with these parameters the repopulation loss exceeds the cell-killing term.",
          "Check the inputs. Comparisons involving a BED that is not positive are not meaningful, and the",
          "search, sensitivity and agreement tabs do not use them.")
  })

  fowler_defaults <- reactive({
    p <- params()
    isTRUE(all.equal(c(p$ab_e, p$alpha_e, p$lambda_e, p$Tk_e), c(10, 0.35, log(2) / 2.5, 7)))
  })

  output$acute_zone <- renderUI({
    # Both zones are shown (audit 2026-10-26): the default acute limit is framed for chemoradiation,
    # so the radiotherapy-alone zone alone would understate the concern. The checkbox sets which
    # zone the search's Fowler limit uses.
    z0 <- fowler_zone(FALSE); z1 <- fowler_zone(TRUE)
    b <- both()
    u <- early_unit()
    txt <- sprintf(paste("Acute tolerance zones of Fowler et al. (2003): %g–%g Gy₁₀ for radiotherapy alone (Fowler 2010",
                         "later gave a 59–63 Gy₁₀ gray zone), and %g–%g Gy₁₀ with concurrent chemotherapy (derived from",
                         "Fowler's possible decrement of 3–5 Gy₁₀). The search's Fowler limit uses the %s zone.",
                         "This regimen: %.1f %s. Comparison: %.1f %s."),
                   z0[1], z0[2], z1[1], z1[2], if (isTRUE(input$chemo)) "chemotherapy" else "radiotherapy-alone",
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
    z <- bench(); x <- z$x; r0 <- z$r0; p <- params(); vs <- regimen_text(ref())
    out <- data.frame(Regimen = z$bm$name, Days = sprintf("%g", z$bm$T),
               a = tumor_delta_txt(x$bed_di, r0$bed_di),
               b = tumor_delta_txt(x$bed_dd, r0$bed_dd),
               `Early BED` = fmt(x$bed_early),
               c = fmt(x$bed_early - r0$bed_early, TRUE),
               `Late BED` = fmt(x$bed_late),
               e = fmt(x$bed_late - r0$bed_late, TRUE), check.names = FALSE)
    names(out)[c(3, 4, 6, 8)] <- c(sprintf("Δ BED DI vs %s (%s)", vs, bed_unit(p$ab)), sprintf("Δ BED DD vs %s (%s)", vs, bed_unit(p$ab)),
                                   sprintf("Δ early BED vs %s (%s)", vs, bed_unit(p$ab_e)), sprintf("Δ late BED vs %s (%s)", vs, bed_unit(p$ab_l)))
    out
  }, striped = TRUE, bordered = TRUE)

  output$hypno_note <- renderUI({
    x <- bench()$x
    one <- function(a, b) {
      t <- tumor_delta_txt(a, b)
      if (startsWith(t, "not")) t else sprintf("%+.1f Gy", a - b)
    }
    dd_txt <- one(x$bed_dd[1], x$bed_dd[2])
    di_txt <- one(x$bed_di[1], x$bed_di[2])
    tagList(
      helpText(sprintf(paste("Δ values are relative to the comparison regimen chosen in the sidebar.",
                             "HYPNO compared its test arm with its own control arm (66 Gy in 33 fractions",
                             "at 6 fractions/week), not with 70 Gy in 35 fractions. With the current",
                             "parameters, test arm minus control arm in tumor BED is %s under DI",
                             "and %s under DD."),
                       di_txt, dd_txt)),
      helpText(paste("HYPNO (792 patients, 75% with weekly cisplatin 35 mg/m², HPV not tested):",
                     "3-year locoregional control 50.7% vs 51.2%; grade 3+ late effects 18.8% vs 20.2%;",
                     "both passed a 10-percentage-point non-inferiority margin. Grade 3+ acute mucositis",
                     "50.9% (190/373) vs 54.7% (208/380). HYPNO tested one schedule; its result does not",
                     "by itself establish any other schedule.")),
      helpText(paste("HYHOPE (phase 1, 24 patients, favorable-risk HPV-positive oropharynx, reduced",
                     "elective volumes, 3 × weekly cisplatin 40 mg/m², designed against 60 Gy in 30",
                     "fractions): grade 3 oral mucositis 5/24, any acute grade 3 7/24, 2-year",
                     "locoregional recurrence 13.8% (95% CI 3.2–32%). Not a comparison for unselected HNSCC.")))
  })

  # Toxicity BEDs do not depend on the DD anchor, so the default anchor is used.
  hypno_caps <- reactive({
    h <- comparator_presets()$hypno_test
    a <- compute_all(h$m, h$d, h$T, params(), default_anchor())
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

  # Settings shared by the candidate and agreement tabs: dose grid and toxicity
  # limits. A limit relative to the comparison is relative to each analysis's
  # own comparison (suggest_regimens computes it).
  tox_args <- reactive({
    early_mode <- input$early_mode
    late_mode <- choice_or("late_mode", "reference")
    cap <- switch(early_mode,
      fowler = fowler_zone(isTRUE(input$chemo))[1],
      hypno = hypno_caps()$early,
      NULL)
    late_cap <- if (late_mode == "hypno") hypno_caps()$late else NULL
    fenwick <- identical(early_mode, "fenwick")
    list(d_step = as.numeric(input$d_step), d_max = checked("d_max"),
         tol_early = if (is.null(cap) && !fenwick) checked("tol_early") else 0,
         tol_late = if (is.null(late_cap)) checked("tol_late") else 0,
         early_cap = cap, late_cap = late_cap, fenwick = fenwick)
  })

  # Candidate-tab settings, which only the candidate search uses.
  search_args <- reactive({
    rf <- ref()
    validate(need(input$m_range[1] < rf$m,
                  "The fraction range must include numbers below the comparison regimen."))
    m_hi <- min(input$m_range[2], rf$m - 1)
    cal <- calendar()
    c(list(m_range = input$m_range[1]:m_hi, criterion = input$criterion, rank_by = input$rank_by,
           min_gain = checked("min_gain"), min_days_saved = checked("min_days")),
      tox_args(),
      list(calendar = list(fx_per_week = cal$fx_per_week, start_weekday = cal$start_weekday)),
      if (!is.null(oar_search())) list(oar = attr(oar_search(), "settings")))
  })

  search <- reactive({
    if (input$early_mode == "fowler" && !fowler_defaults()) {
      s <- empty_suggestions()
      attr(s, "reason") <- "fowler_params_changed"
      return(list(s = s, cap = NULL, blocked = TRUE))
    }
    args <- search_args()
    args$oar <- NULL
    args$oar_check <- oar_search()
    s <- do.call(suggest_regimens, c(list(p = params(), ref = ref(), anchor = anchor()), args))
    list(s = s, cap = args$early_cap, late_cap = args$late_cap, blocked = FALSE)
  })
  suggestions <- reactive(search()$s)

  # One line naming what an analysis compares against and with which limits.
  limits_text <- function(rf) {
    p <- params(); ta <- tox_args()
    r0 <- compute_all(rf$m, rf$d, rf$T, p, anchor())
    u_e <- early_unit(); u_l <- bed_unit(p$ab_l)
    early <- if (isTRUE(ta$fenwick)) "Fenwick 2008 boundary (cell-kill BED, depends on overall time)" else
      if (is.null(ta$early_cap)) sprintf("%.2f %s (comparison %.2f + %g)", r0$bed_early + ta$tol_early, u_e, r0$bed_early, ta$tol_early) else
      sprintf("%.2f %s (%s)", ta$early_cap, u_e, if (input$early_mode == "hypno") "HYPNO 55 Gy/20, fixed" else "Fowler zone lower edge, fixed")
    late <- if (is.null(ta$late_cap)) sprintf("%.2f %s (comparison %.2f + %g)", r0$bed_late + ta$tol_late, u_l, r0$bed_late, ta$tol_late) else
      sprintf("%.2f %s (HYPNO 55 Gy/20, fixed)", ta$late_cap, u_l)
    sprintf("Early toxicity limit: %s. Late toxicity limit: %s. Dose grid: multiples of %s Gy up to %g Gy.",
            early, late, input$d_step, ta$d_max)
  }
  context_strip <- function(...) div(class = "alert alert-light", style = "padding: 0.5rem 0.75rem;", ...)
  regimen_txt <- function(r) sprintf("%d × %g Gy in %g days", r$m, r$d, r$T)
  days_txt <- function(n) sprintf("%g %s", n, if (n == 1) "day" else "days")

  output$search_context <- renderUI({
    rf <- ref(); an <- anchor(); p <- params()
    context_strip(sprintf(paste("Comparison: %s (sidebar). DD anchor: %s. Tumor α/β %g Gy, mucosa %g Gy,",
                                "late tissue %g Gy. Tumor criterion: change of at least %g Gy, and at least %s saved."),
                          regimen_txt(rf), regimen_txt(an), p$ab, p$ab_e, p$ab_l, checked("min_gain"), days_txt(checked("min_days"))),
                  tags$br(), if (isTRUE(search()$blocked)) "Early toxicity limit: not available (see below)." else limits_text(rf))
  })

  output$sens_context <- renderUI({
    r <- regimen(); rf <- ref(); an <- anchor()
    context_strip(sprintf(paste("Regimen evaluated: %s (the sidebar regimen, not a candidate row). Comparison: %s.",
                                "DD anchor: %s. Summary threshold: a change of at least %g Gy (the global tumor",
                                "threshold)."), regimen_txt(r), regimen_txt(rf), regimen_txt(an), checked("min_gain")))
  })

  output$search_note <- renderUI({
    p <- params(); sr <- search()
    if (isTRUE(sr$blocked)) return(NULL)
    early_txt <- if (identical(input$early_mode, "fenwick")) "cell-kill BED (α/β 10 Gy) at or below the Fenwick et al. 2008 boundary for its overall time" else
      if (is.null(sr$cap)) "early BED at most the comparison's plus the allowed excess" else
      if (input$early_mode == "hypno") sprintf("early BED at most %.1f %s (HYPNO 55 Gy/20 fractions, same acute model)", sr$cap, early_unit()) else
      sprintf("early BED at most %g Gy₁₀ (Fowler zone, lower edge%s)", sr$cap,
              if (isTRUE(input$chemo)) ", with concurrent chemotherapy" else "")
    late_txt <- if (is.null(sr$late_cap)) "late BED at most the comparison's plus the allowed excess" else
      sprintf("late BED at most %.1f %s (HYPNO 55 Gy/20 fractions)", sr$late_cap, bed_unit(p$ab_l))
    helpText(sprintf(paste("Searched: once-daily Mon-Fri regimens with a dose per fraction above the",
                           "comparison's (multiples of %s Gy, up to %g Gy), %s, and %s. For each number",
                           "of fractions, the highest-ranked dose per fraction is listed; the ZIP download",
                           "also holds every qualifying regimen. Δ = regimen minus comparison. Rows above",
                           "%g Gy per fraction are above the largest fraction size among the benchmark regimens shown on the Results tab.",
                           "The list depends strongly on the tumor repopulation settings (now %.2f Gy/day under DI)."),
                     input$d_step, checked("d_max"), early_txt, late_txt, PRECEDENT_MAX_D,
                     p$lambda_di / p$alpha))
  })

  pretty_suggestions <- reactive({
    s <- suggestions()
    dd_col <- if (nrow(s) > 0 && all(is.na(s$d_dd))) rep("not available", nrow(s)) else s$d_dd
    vs <- regimen_text(ref()); p <- params()
    out <- data.frame(Fractions = s$m, `Gy per fraction` = s$d, Days = s$T, `Days saved` = ref()$T - s$T,
               `Fractions saved` = s$fractions_saved, `Total Gy` = s$total_dose,
               `Δ BED DI` = s$d_di, `Δ BED DD` = dd_col, `Credited by` = s$credited_by,
               `Δ early BED` = s$d_early, `Δ late BED` = s$d_late,
               `Early EQD2` = s$eqd2_early, `Late EQD2` = s$eqd2_late,
               `Early limit minus regimen` = s$early_slack, `Late limit minus regimen` = s$late_slack,
               flag = ifelse(s$beyond_precedent, "yes", ""),
               `Outside the tumor data` = ifelse(s$d > tumor_data_support()$d_max + 1e-9, "yes", ""),
               `Next dose step would exceed` = s$next_step, check.names = FALSE)
    names(out)[names(out) == "flag"] <- PRECEDENT_HEADER
    for (k in c("DI", "DD", "early", "late")) {
      u <- switch(k, DI = , DD = bed_unit(p$ab), early = bed_unit(p$ab_e), late = bed_unit(p$ab_l))
      names(out)[names(out) == sprintf("Δ %s BED", k) | names(out) == sprintf("Δ BED %s", k)] <-
        sprintf("Δ %s vs %s (%s)", if (k %in% c("DI", "DD")) paste("BED", k) else paste(k, "BED"), vs, u)
    }
    out
  })

  output$suggest_msg <- renderUI({
    s <- suggestions()
    if (nrow(s) > 0) return(NULL)
    reason <- attr(s, "reason")
    msg <- function(...) div(class = "alert alert-info", ...)
    if (!is.null(reason)) return(msg(reason_text(reason, input$criterion)))
    if (identical(input$criterion, "both"))
      return(msg(paste(
        "No regimen in the searched range meets the selected minimum tumor BED change under both models",
        "and all other limits. An empty result is a valid outcome of this screen. The Both-model screen uses",
        "the global comparison and threshold by default, and you can override them there. An allowed BED loss",
        "changes the hypothesis being screened. To explore an allowed loss here, explicitly enter a negative minimum change.",
        "\u201cDI or DD model\u201d requires the threshold under at least one evaluable model. It can include",
        "regimens that fail the threshold under the other model. Changing these settings changes the",
        "hypothesis being screened.")))
    msg("No regimen in the searched range meets all limits with these settings.",
        "This is a result, not an error: with these parameters, shortening the course costs",
        "tumor BED faster than the limits allow.",
        if (input$criterion != "either") "Options: require the change under 'DI or DD model'," else "Options:",
        "allow a tumor BED loss (a negative minimum change), lower the days saved,",
        "or choose another early or late limit.")
  })

  output$suggest_table <- renderDT({
    s <- pretty_suggestions()
    two_dp <- c("Gy per fraction", "Total Gy", grep("^Δ ", names(s), value = TRUE), "Early EQD2", "Late EQD2",
                "Early limit minus regimen", "Late limit minus regimen")
    two_dp <- two_dp[vapply(s[two_dp], is.numeric, logical(1))]
    datatable(s, rownames = FALSE, fillContainer = FALSE,
              options = list(pageLength = 15, dom = "tip", scrollX = TRUE)) |>
      formatRound(two_dp, 2)
  })

  # Interface choices saved with an export. They do not change the result, so
  # an invalid Results-tab regimen is recorded as absent rather than blocking.
  export_context <- function() {
    list(comparator = choice_or("comparator", "custom"), preset = input$preset,
         early_mode = input$early_mode, late_mode = choice_or("late_mode", "reference"),
         chemo = isTRUE(input$chemo), calendar = calendar(), dd_mode = choice_or("dd_mode", "calibrated"),
         regimen_evaluated = tryCatch(regimen(), error = function(e) NULL))
  }

  scenario <- reactive({
    sr <- search()
    validate(need(!isTRUE(sr$blocked), "Nothing to export with the current settings."))
    build_scenario(p = params(), ref = ref(), anchor = anchor(), search = search_args(),
                   context = export_context(), analysis = "candidate", code = CODE_AT_LOAD)
  })

  output$download <- downloadHandler(
    filename = function() export_filename("candidate", scenario_id(scenario())),
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
      scale_colour_manual(values = c(DI = "#0072B2", DD = "#D55E00")) +   # Okabe-Ito, colorblind-safe
      scale_shape_manual(values = c(DI = 16, DD = 17)) +
      scale_linetype_manual(values = c(DI = "solid", DD = "dashed")) +
      labs(x = "Number of fractions (listed dose per fraction for each)",
           y = sprintf("Tumor BED minus %s (%s)", regimen_text(ref()), bed_unit(params()$ab)),
           title = "Tumor BED difference from the comparison for the listed regimens") +
      theme_minimal(base_size = 15) +
      theme(axis.title = element_text(face = "bold"), legend.position = "top")
  })

  agree_ref <- reactive({
    if (!identical(input$agree_cmp, "hypno_test")) { r <- ref(); return(list(m = r$m, d = r$d, T = r$T)) }
    r <- comparator_presets()$hypno_test
    list(m = r$m, d = r$d, T = r$T)
  })
  agree_anchor <- reactive(fixed_anchor())

  # Arguments of the agreement search: the shared toxicity settings and the
  # agreement tab's own fraction range and days saved.
  agree_args <- reactive({
    validate(need(!(input$early_mode == "fowler" && !fowler_defaults()),
                  "The Fowler limit needs Fowler's acute parameters. Restore them or choose another early limit."))
    rf <- agree_ref()
    validate(need(rf$m > 1, "The comparison needs at least 2 fractions."))
    c(list(m_range = 1:(rf$m - 1), min_days_saved = checked("agree_days")), tox_args(),
      list(oar_check = oar_search()))
  })
  # The both-model screen's threshold: the global one unless explicitly overridden.
  agree_threshold <- reactive(if (identical(input$agree_threshold_mode, "own")) checked("agree_margin") else checked("min_gain"))
  output$agree_override <- renderUI({
    o <- c(if (identical(input$agree_cmp, "hypno_test")) sprintf("comparison %s instead of the global %s",
                                                                regimen_text(agree_ref()), regimen_text(ref())),
           if (identical(input$agree_threshold_mode, "own")) sprintf("threshold %g Gy instead of the global %g Gy",
                                                                    agree_threshold(), checked("min_gain")))
    if (length(o)) div(class = "alert alert-warning bed-override", tags$strong("Override on this tab: "),
                       paste(o, collapse = "; "), ".")
  })

  agreement <- reactive({
    do.call(agreement_regimens, c(list(p = params(), ref = agree_ref(), anchor = agree_anchor(),
                                       margin = agree_threshold()), agree_args()))
  })

  output$agree_note <- renderUI({
    rf <- agree_ref(); p <- params(); an <- agree_anchor()
    helpText(sprintf(paste("Comparison: %d × %g Gy in %g days. DD anchor: %d × %g Gy in %g days (%s).",
                           "Toxicity limits: as chosen on the candidate tab (a limit relative to the comparison",
                           "is relative to this regimen). DI accelerated repopulation starts on day %g, so a course",
                           "ending well before that day gains little from further shortening under DI.%s Under DD,",
                           "onset and rate depend on cell killing, so shortening can save part of the repopulation or",
                           "cost more, depending on the regimen."),
                     rf$m, rf$d, rf$T, an$m, an$d, an$T,
                     "fixed", p$Tk,
                     if (p$g > 0) sprintf(" Background growth g still gives %.3f %s per day saved, before and after that day.",
                                          p$g / p$alpha, bed_unit(p$ab)) else ""))
  })

  output$agree_context <- renderUI({
    rf <- agree_ref()
    context_strip(sprintf("Comparison: %s (%s). Allowed change: at least %g Gy under both models, and at least %s saved.",
                          regimen_txt(rf),
                          if (identical(input$agree_cmp, "hypno_test")) "override: HYPNO test arm" else "the global comparison",
                          agree_threshold(), days_txt(checked("agree_days"))),
                  tags$br(), limits_text(rf))
  })

  output$agree_msg <- renderUI({
    a <- agreement()
    if (nrow(a) > 0) return(NULL)
    reason <- attr(a, "reason")
    div(class = "alert alert-info",
        if (!is.null(reason)) reason_text(reason, "both") else
          paste("No regimen meets the allowed change under both models within the toxicity limits.",
                "Try a more negative allowed change, fewer days saved, or another toxicity limit."))
  })

  output$agree_table <- renderTable({
    a <- agreement()
    validate(need(nrow(a) > 0, ""))
    p <- params()
    rng <- function(lo, hi) ifelse(is.na(lo), "not evaluable", sprintf("%+.2f to %+.2f", lo, hi))
    counts <- ifelse(a$n_invalid > 0,
                     sprintf("%d of %d evaluable settings (%d of %d not evaluable)", a$n_pass_both, a$n_valid_both,
                             a$n_invalid, a$n_scenarios),
                     sprintf("%d of %d settings", a$n_pass_both, a$n_valid_both))
    data.frame(Fractions = sprintf("%d", a$m), `Gy per fraction` = sprintf("%.2f", a$d),
               Days = sprintf("%g", a$T), `Total Gy` = sprintf("%.2f", a$total_dose),
               `Δ BED DI vs the screen's comparison` = fmt(a$d_di, TRUE), `Δ BED DD vs the screen's comparison` = fmt(a$d_dd, TRUE),
               `Δ DI over evaluable settings` = rng(a$d_di_min, a$d_di_max),
               `Δ DD over evaluable settings` = rng(a$d_dd_min, a$d_dd_max),
               `Selected parameter settings with both ≥ allowed change` = counts,
               `DD onset (day)` = ifelse(is.finite(a$dd_onset), sprintf("%.1f", a$dd_onset), "none"),
               `Ends before DI onset` = ifelse(a$before_di_onset, "yes", "no"),
               `Early BED` = sprintf("%.2f %s", a$bed_early, bed_unit(p$ab_e)),
               `Late BED` = sprintf("%.2f %s", a$bed_late, bed_unit(p$ab_l)),
               flag = ifelse(a$beyond_precedent, "yes", ""),
               check.names = FALSE) |> (\(x) { names(x)[names(x) == "flag"] <- PRECEDENT_HEADER; x })()
  }, striped = TRUE, bordered = TRUE)

  scenario_agree <- reactive({
    args <- c(agree_args(), list(criterion = "both", rank_by = "worst", min_gain = agree_threshold()))
    args$oar_check <- NULL
    if (!is.null(oar_search())) args$oar <- attr(oar_search(), "settings")
    build_scenario(p = params(), ref = agree_ref(), anchor = agree_anchor(), search = args,
                   context = c(export_context(), list(agree_cmp = input$agree_cmp)),
                   analysis = "agreement",
                   agreement = list(margin = agree_threshold(), ab_values = SCENARIO_AB_VALUES),
                   code = CODE_AT_LOAD)
  })

  output$download_agree <- downloadHandler(
    filename = function() export_filename("agreement", scenario_id(scenario_agree())),
    content = function(file) {
      a <- agreement()
      write_export_zip(file, scenario_agree(), a, attr(a, "feasible"))
    },
    contentType = "application/zip")

  sens <- reactive(sensitivity_table(regimen(), ref(), anchor(), params()))

  output$sens_summary <- renderUI({
    st <- sens(); mg <- checked("min_gain")
    n <- nrow(st)
    di_n <- sum(st$valid_di & st$d_di >= mg - 1e-9, na.rm = TRUE)
    dd_n <- sum(st$valid_dd & st$d_dd >= mg - 1e-9, na.rm = TRUE)
    bad <- n - sum(st$valid_di & st$valid_dd)
    tags$p(sprintf(paste("Tumor BED change at least %g Gy (the global tumor threshold):",
                         "%d of %d evaluable parameter settings under DI, %d of %d evaluable parameter settings under DD.%s",
                         "These selected parameter settings are not independent confirmations or a probability sample.",
                         "Within each setting, DI and DD share \u03b1 and DD is calibrated to DI at the anchor. \u03b1 and the",
                         "calibrated constant can differ between settings. Uncertainty intervals are not calculated."),
                   mg, di_n, sum(st$valid_di), dd_n, sum(st$valid_dd),
                   if (bad > 0) sprintf(paste(" %d of %d settings are not evaluable under at least one model",
                                              "(a comparison or regimen BED is not positive, or DD is unavailable)."),
                                        bad, n) else ""))
  })

  show_delta <- function(d, reason) {
    ifelse(reason == "ok", sprintf("%+.2f", d),
      ifelse(reason == "negative_reference", "not evaluable (comparison BED ≤ 0)",
        ifelse(reason == "negative_candidate", "not evaluable (regimen BED ≤ 0)", "not available (no DD calibration)")))
  }

  output$sens_table <- renderTable({
    st <- sens()
    data.frame(Preset = unname(PRESET_SHORT[st$preset]), `Tumor α/β (Gy)` = sprintf("%g", st$ab),
               `α (1/Gy)` = sprintf("%.3g", st$alpha), `DI loss (Gy/day)` = sprintf("%.2f", st$lambda_di / st$alpha),
               `DI onset (day, inclusive)` = sprintf("%g", st$Tk),
               `DD max rate (1/day)` = sprintf("%.3g%s", st$lambda_dd, ifelse(st$lambda_dd_inherited, " (sidebar)", "")),
               `g (1/day)` = sprintf("%g", st$g),
               `Δ BED DI` = show_delta(st$d_di, st$reason_di), `Δ BED DD` = show_delta(st$d_dd, st$reason_dd),
               check.names = FALSE)
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

  # --- version 3: shared settings, modules, banner, report card, comparison ----

  # Everything the modules need, from one place.
  settings <- reactive({
    ta <- tox_args()
    list(p = params(), regimen = regimen(), ref = ref(), threshold = checked("min_gain"),
         early_mode = input$early_mode, late_mode = choice_or("late_mode", "reference"),
         tol_early = ta$tol_early, tol_late = ta$tol_late, calendar = calendar(),
         d_step = as.numeric(input$d_step), d_max = checked("d_max"),
         m_range = input$m_range[1]:max(input$m_range[1], min(input$m_range[2], ref()$m - 1)),
         min_days = checked("min_days"))
  })
  oar <- oar_server("oar", settings)
  # The organ policy of every search (entered as % of prescription, and only
  # when the box on the Organs at risk tab is ticked). NULL means not applied.
  # When it is applied, an invalid organ input or an entered organ without a
  # verified constraint blocks the searches with the reason; it never turns
  # into "no check".
  oar_policy <- reactive({
    applied <- isTRUE(input[["oar-use_in_search"]]) && identical(input[["oar-mode"]], "percent")
    if (!applied) return(NULL)
    o <- tryCatch(oar(), error = function(e) e)
    if (inherits(o, "error"))
      validate(need(FALSE, paste("Organ limits are applied (Organs at risk tab), but an organ input is invalid:",
                                 conditionMessage(o), "Correct it, or untick the box to search without organ limits.")))
    if (length(o$percents) == 0) return(NULL)
    lib <- load_constraints()
    miss <- names(o$percents)[vapply(names(o$percents), function(k) is.null(pick_constraint(lib, k, o$library)), logical(1))]
    if (length(miss))
      validate(need(FALSE, sprintf(paste("Organ limits are applied, but %s %s no verified constraint in the %s library.",
                                         "Clear %s, or choose another constraint source."),
                                   paste(OAR_ORGANS$label[match(miss, OAR_ORGANS$key)], collapse = ", "),
                                   if (length(miss) > 1) "have" else "has", o$library,
                                   if (length(miss) > 1) "those doses" else "its dose")))
    list(library = o$library, percents = o$percents, ab = o$ab, coverage = oar_coverage(lib, o$percents, o$library))
  })
  unc_server("unc", settings, oar_policy)
  oar_search <- reactive({
    pol <- oar_policy()
    if (is.null(pol)) return(NULL)
    build_oar_check(load_constraints(), pol$percents, pol$library, pol$ab)
  })

  output$banner <- renderUI({
    r <- tryCatch(regimen(), error = function(e) NULL); rf <- tryCatch(ref(), error = function(e) NULL)
    if (is.null(r) || is.null(rf)) return(div(class = "alert alert-danger", "Check the regimen inputs."))
    p <- params(); an <- anchor()
    dd <- if (identical(p$dd_mode, "fitted")) "DD values entered directly (not calibrated)" else sprintf("DD calibrated at %s", regimen_text(an))
    pre <- names(PRESET_LABELS)[PRESET_LABELS == input$preset]
    sup <- tumor_data_support()
    extrap <- if (r$d > sup$d_max + 1e-9 || r$d < sup$d_min - 1e-9)
      tags$div(tags$strong("Extrapolation: "), sprintf("%g Gy per fraction is outside the tumor data (%g-%g Gy per fraction).",
                                                      r$d, sup$d_min, sup$d_max))
    div(class = "alert alert-primary", style = "padding:0.45rem 0.8rem; margin-bottom:0.5rem;",
        tags$strong("Evaluating "), regimen_text(r), tags$strong(" against "), regimen_text(rf),
        sprintf(" · tumor threshold %+g Gy · %s · tumor preset: %s", checked("min_gain"), dd,
                if (length(pre)) sub(" \\(.*$", "", sub(":.*$", "", pre)) else input$preset), extrap)
  })

  limit_values <- reactive({
    ta <- tox_args(); rf <- ref(); p <- params()
    r0 <- compute_all(rf$m, rf$d, rf$T, p, anchor())
    early <- if (isTRUE(ta$fenwick)) NA_real_ else if (is.null(ta$early_cap)) r0$bed_early + ta$tol_early else ta$early_cap
    late <- if (is.null(ta$late_cap)) r0$bed_late + ta$tol_late else ta$late_cap
    list(early = early, late = late, fenwick = isTRUE(ta$fenwick))
  })
  output$early_limit_value <- renderUI({
    lv <- limit_values()
    helpText(if (lv$fenwick) "Value: depends on each regimen's overall time (Fenwick boundary)." else
      sprintf("Value now: %.2f %s.", lv$early, bed_unit(params()$ab_e)))
  })
  output$late_limit_value <- renderUI(helpText(sprintf("Value now: %.2f %s.", limit_values()$late, bed_unit(params()$ab_l))))

  output$cisplatin_context <- renderUI({
    r <- regimen(); rf <- ref()
    wk <- function(T) ceiling(T / 7)
    helpText(sprintf(paste("Context only: with weekly cisplatin the evaluated course spans %d treatment weeks",
                           "(%d for the comparison). HYPNO gave 4 vs 5 weekly courses. Cumulative cisplatin is not modeled."),
                     wk(r$T), wk(rf$T)))
  })

  output$dd_warning0 <- renderUI(dd_message())

  # Report card: one row per quantity, with value, difference, limit and status.
  output$report <- renderTable({
    b <- both(); r <- regimen(); rf <- ref(); p <- params(); lv <- limit_values()
    vs <- regimen_text(rf)
    unc <- tryCatch({
      di <- load_draws("di"); dd <- load_draws("dd")
      summarize_draws(data.frame(m = r$m, d = r$d, T = r$T),
                      tumor_deltas(data.frame(m = r$m, d = r$d, T = r$T), rf, di, dd))
    }, error = function(e) NULL)
    rng <- function(med, lo, hi) if (is.null(unc)) "not available" else sprintf("%s (%s to %s)", fmt_half(med), fmt_half(lo), fmt_half(hi))
    st <- function(x, lim) if (is.na(lim)) "see Results tab" else if (x <= lim + 1e-9) sprintf("within (limit %.1f)", lim) else sprintf("exceeds (limit %.1f)", lim)
    dd_same <- !identical(p$dd_mode, "fitted") && isTRUE(all.equal(c(rf$m, rf$d, rf$T), c(anchor()$m, anchor()$d, anchor()$T)))
    early_st <- if (lv$fenwick) {
      lim <- fenwick_limit(r$T)
      if (is.na(lim)) "outside the Fenwick range" else st(fenwick_cell_kill(r$m, r$d), lim)
    } else st(b$bed_early[1], lv$early)
    rows <- data.frame(
      Quantity = c(sprintf("Tumor BED, DI (%s)", bed_unit(p$ab)), sprintf("Tumor BED, DD (%s)", bed_unit(p$ab)),
                   sprintf("Acute mucosal BED (%s)", bed_unit(p$ab_e)), sprintf("Late BED at the prescription dose (%s)", bed_unit(p$ab_l))),
      `This regimen` = fmt(c(b$bed_di[1], b$bed_dd[1], b$bed_early[1], b$bed_late[1])),
      Comparison = c(fmt(b$bed_di[2]), paste0(fmt(b$bed_dd[2]), if (dd_same) " (equal to DI by calibration)" else ""),
                     fmt(b$bed_early[2]), fmt(b$bed_late[2])),
      d = c(tumor_delta_txt(b$bed_di[1], b$bed_di[2]), tumor_delta_txt(b$bed_dd[1], b$bed_dd[2]),
            fmt(b$bed_early[1] - b$bed_early[2], TRUE), fmt(b$bed_late[1] - b$bed_late[2], TRUE)),
      u = c(if (is.null(unc)) "not available" else rng(unc$di_median, unc$di_q10, unc$di_q90),
            if (is.null(unc)) "not available" else paste(rng(unc$dd_median, unc$dd_q10, unc$dd_q90), "(calibrated DD)"),
            "see Uncertainty tab", "see Uncertainty tab"),
      `Limit and status` = c(sprintf("threshold %+g: %s", checked("min_gain"),
                                     ifelse(is.finite(b$bed_di[1] - b$bed_di[2]) && b$bed_di[1] - b$bed_di[2] >= checked("min_gain") - 1e-9, "met", "not met")),
                             sprintf("threshold %+g: %s", checked("min_gain"),
                                     ifelse(is.finite(b$bed_dd[1] - b$bed_dd[2]) && b$bed_dd[1] - b$bed_dd[2] >= checked("min_gain") - 1e-9, "met", "not met")),
                             early_st, st(b$bed_late[1], lv$late)),
      check.names = FALSE)
    names(rows)[4:5] <- c(sprintf("Δ vs %s", vs), "Verified-trial refit 2026-10-25: Δ median (80% range), Gy at each draw's α/β")
    o <- tryCatch(oar()$table, error = function(e) e)
    if (inherits(o, "error")) {
      rows <- rbind(rows, stats::setNames(data.frame("Organs at risk", "", "", "", "",
                                                     paste("not assessed: an organ input is invalid.", conditionMessage(o))), names(rows)))
    } else if (!is.null(o)) {
      o <- o[o$entered, ]
      if (nrow(o)) {
        mean_m <- o$metric == "Dmean"
        rows <- rbind(rows, data.frame(
          Quantity = ifelse(mean_m, sprintf("%s, mean physical dose (Gy)", o$organ),
                            sprintf("%s, maximum EQD2 (Gy, α/β %g)", o$organ, o$ab)),
          `This regimen` = ifelse(is.na(o$planned_physical), "", sprintf("%.1f", ifelse(mean_m, o$planned_physical, o$planned_eqd2))),
          Comparison = "", d = "", u = "",
          `Limit and status` = paste0(
            ifelse(is.na(o$limit_physical), o$status,
                   sprintf("%s (limit %s%.1f)", o$status, o$operator, ifelse(mean_m, o$limit_physical, o$limit_eqd2))),
            # audit 2026-10-26: the status alone can hide an alpha/beta-dependent result or the
            # cross-fractionation limitation of a physical mean comparison
            ifelse(mean_m, "; physical source limit, the effect of the changed fractionation is not validated",
                   ifelse(nzchar(o$flip), paste0("; over the α/β range: ", o$flip), ""))),
          check.names = FALSE) |>
          stats::setNames(names(rows)))
      }
    }
    rows
  }, striped = TRUE, bordered = TRUE)

  output$report_notes <- renderUI({
    r <- regimen(); sup <- tumor_data_support()
    tagList(
      if (r$d > sup$d_max + 1e-9 || r$d < sup$d_min - 1e-9)
        div(class = "alert alert-warning", sprintf(paste("%g Gy per fraction is outside the tumor data the models were fitted to (%g-%g Gy",
                                                         "per fraction): tumor BEDs at this fraction size are extrapolations."), r$d, sup$d_min, sup$d_max)),
      if (isTRUE(tryCatch(any(oar()$table$entered), error = function(e) FALSE)))
        helpText(oar_coverage_text(oar()$coverage)),
      helpText("Organ rows appear when organ doses are entered on the Organs at risk tab. The refit column applies the",
               sprintf(paste("refit's parameter draws (not the sidebar values) to the same two regimens: DI posterior draws and",
                             "calibrated-DD draws, fitted to %d arms of %d randomized trials (%s). Its BED differences use each",
                             "draw's α/β, so they are not in Gy₁₀."), sup$n_arms, sup$n_studies, sup$file)))
  })

  # Compare regimens side by side (quantities as rows).
  output$compare_table <- renderTable({
    p <- params(); rf <- ref(); r <- regimen(); a <- anchor()
    regs <- list(list(name = paste("Evaluated:", regimen_text(r)), m = r$m, d = r$d, T = r$T))
    for (k in 1:3) {
      m <- checked(sprintf("cmp%d_m", k)); d <- checked(sprintf("cmp%d_d", k))
      regs[[length(regs) + 1]] <- list(name = sprintf("%s: %s", LETTERS[k], regimen_text(list(m = m, d = d, T = overall_time(m)))),
                                       m = m, d = d, T = overall_time(m))
    }
    bm <- benchmark_regimens()
    for (i in seq_len(nrow(bm))) regs[[length(regs) + 1]] <- list(name = bm$name[i], m = bm$m[i], d = bm$d[i], T = bm$T[i])
    r0 <- compute_all(rf$m, rf$d, rf$T, p, a)
    cols <- lapply(regs, function(x) {
      v <- compute_all(x$m, x$d, x$T, p, a)
      c(sprintf("%d × %g Gy", x$m, x$d), sprintf("%g", x$m * x$d), sprintf("%g", x$T),
        tumor_delta_txt(v$bed_di, r0$bed_di), tumor_delta_txt(v$bed_dd, r0$bed_dd),
        sprintf("%.1f", v$bed_early), sprintf("%.1f", v$bed_late),
        if (x$d > tumor_data_support()$d_max + 1e-9) "yes" else "no")
    })
    out <- data.frame(Quantity = c("Schedule", "Total dose (Gy)", "Overall time (days)",
                                   sprintf("Δ tumor BED DI vs %s (%s)", regimen_text(rf), bed_unit(p$ab)),
                                   sprintf("Δ tumor BED DD vs %s (%s)", regimen_text(rf), bed_unit(p$ab)),
                                   sprintf("Acute BED (%s)", bed_unit(p$ab_e)), sprintf("Late BED (%s)", bed_unit(p$ab_l)),
                                   "Dose per fraction outside the tumor data"),
                      do.call(cbind, cols), check.names = FALSE)
    names(out)[-1] <- vapply(regs, `[[`, "", "name")
    out
  }, striped = TRUE, bordered = TRUE)

  output$about <- renderUI({
    txt <- paste(readLines("about.md", encoding = "UTF-8", warn = FALSE), collapse = "\n")
    txt <- gsub("{{VERSIONS}}", sprintf("Version %s of the Custom BED calculator (equations %s, selection rules %s, scenario set %s, uncertainty %s, organs %s).",
                                        APP_VERSION, EQUATION_VERSION, SELECTION_VERSION, SCENARIO_SET_VERSION, UNCERTAINTY_VERSION, OAR_VERSION),
                txt, fixed = TRUE)
    rev <- if (file.exists("deployed_revision.txt")) paste0("Deployed revision ", substr(readLines("deployed_revision.txt", n = 1), 1, 7), ".") else
      "Local copy (not a deployed revision)."
    txt <- gsub("{{REVISION}}", rev, txt, fixed = TRUE)
    withMathJax(HTML(markdown::mark(text = txt, options = "+latex_math")))
  })
}

shinyApp(ui, server)
