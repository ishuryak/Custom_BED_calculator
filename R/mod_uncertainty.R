# Uncertainty tab (Shiny module). Requires R/uncertainty.R and R/oar.R.
# settings(): list(regimen, ref, threshold, early_mode, late_mode, tol_early,
#   tol_late, calendar, d_step, d_max, m_range, min_days).
# oar(): the organ module's value (percents, search_check, library, ab).

UNC_VARIANTS <- c("Primary refit (α/β prior 5-20 Gy, background growth estimated, overdispersion)" = "primary",
                  "2018 paper's restrictions (α/β 10 Gy and g 0.0154/day fixed, overdispersion)" = "paper_restrictions")

unc_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "alert alert-secondary",
        tags$strong("What this tab shows. "),
        "Each tumor model was refitted to the 35 randomized trial arms of Shuryak, Hall and Brenner (2018) with",
        "Bayesian methods (one intercept per trial, so only within-trial contrasts inform the fit, plus",
        "between-arm overdispersion). The refit's joint parameter draws are applied to the evaluated regimen and the",
        "comparison together, so their difference carries the parameter correlations. Results are shares of",
        "simulated parameter draws under the stated models and priors, not probabilities of tumor control or of",
        "safety. DI and DD are shown separately and never averaged."),
    uiOutput(ns("refit_status")),
    layout_columns(col_widths = c(4, 4, 4),
      selectInput(ns("variant"), "Refit used", UNC_VARIANTS),
      selectInput(ns("tumor_crit"), "Tumor criterion for the joint share",
                  c("DI model (refit informed by data)" = "di", "both models" = "both",
                    "DD model (refit prior-dominated)" = "dd")),
      radioButtons(ns("tox_source"), "Toxicity parameters",
                   c("Uncertain (stated priors)" = "prior", "Fixed (the sidebar values)" = "fixed"))),
    uiOutput(ns("context")),
    h5("Evaluated regimen against the comparison"),
    div(style = "overflow-x:auto;", tableOutput(ns("summary"))),
    plotOutput(ns("hist"), height = "320px"),
    helpText("Each bar counts parameter draws. Grey line: zero difference; dashed line: the tumor threshold."),
    h5("Tumor control, conditional on a stated control rate"),
    layout_columns(col_widths = c(4, 8),
      numericInput(ns("tcp0"), "Tumor control of the comparison regimen (%)", 70, 1, 99, 1),
      helpText("With tumor control exp(−K S), a BED difference ΔBED multiplies S by exp(−α ΔBED), so the evaluated",
               "regimen's control is (control of the comparison)^exp(−α ΔBED), with α from the same draw. The refit's",
               "α is population-level (it includes patient-to-patient variation), so these slopes are shallow. The",
               "control rate is your input: the refit's trials did not use concurrent chemotherapy.")),
    div(style = "overflow-x:auto;", tableOutput(ns("tcp"))),
    h5("Search under uncertainty"),
    helpText("The search uses 1000 of the 4000 draws; the listed regimens are then summarized with all 4000. The",
             "selection rule is fixed before the search. The best value over draws is never reported, because it is",
             "optimistic by construction. Each listed regimen is shown with its neighbors on the dose grid."),
    layout_columns(col_widths = c(4, 4, 4),
      selectInput(ns("rule"), "Selection rule, per number of fractions",
                  c("Highest 10th percentile of the smaller of ΔDI and ΔDD" = "min_q10",
                    "Highest share of draws meeting the threshold under both models and every limit" = "share_joint")),
      numericInput(ns("tox_share"), "Required share of draws within all toxicity limits (%)", 90, 50, 100, 5),
      actionButton(ns("run"), "Run the search", class = "btn-primary", style = "margin-top:1.9rem;")),
    uiOutput(ns("search_msg")),
    div(style = "overflow-x:auto;", tableOutput(ns("search_table"))),
    downloadButton(ns("download"), "Download the uncertainty analysis (ZIP)"))
}

unc_server <- function(id, settings, oar) {
  moduleServer(id, function(input, output, session) {
    draws <- reactive(list(di = load_draws("di", input$variant), dd = load_draws("dd", input$variant)))
    tx_all <- reactive({
      s <- settings()
      if (identical(input$tox_source, "fixed")) fixed_toxicity_draws(s$p, nrow(draws()$di)) else toxicity_draws(nrow(draws()$di))
    })
    limits <- reactive({
      s <- settings()
      validate(need(s$early_mode %in% c("hypno", "reference"),
                    "Under uncertainty the early limit must be HYPNO 55 Gy/20 or the comparison plus an allowed excess (chosen on the Search tab): the Fowler and Fenwick limits are tied to fixed parameter values."))
      list(early = list(mode = s$early_mode, tol = s$tol_early), late = list(mode = s$late_mode, tol = s$tol_late))
    })
    oar_fun <- reactive({
      o <- oar()
      if (!isTRUE(o$search_check) || length(o$percents) == 0) return(NULL)
      oar_draw_checker(load_constraints(), o$percents, o$library)
    })

    one <- reactive({
      s <- settings(); r <- s$regimen
      regs <- data.frame(m = r$m, d = r$d, T = r$T)
      td <- tumor_deltas(regs, s$ref, draws()$di, draws()$dd)
      tox <- toxicity_ok_draws(regs, s$ref, tx_all(), limits(), if (is.null(oar_fun())) NULL else oar_fun()(regs, tx_all()))
      list(td = td, tox = tox, sum = summarize_draws(regs, td, tox, s$threshold))
    })

    output$refit_status <- renderUI({
      d <- draws(); sup <- tumor_data_support()
      la <- stats::quantile(d$di$lambda / d$di$alpha, c(.025, .975), names = FALSE)
      div(class = "alert alert-warning",
          tags$strong("Status of the refit (for review, not yet the default). "),
          sprintf(paste("The DI refit is informed by the data (daily repopulation loss λ/α, 95%% %.2f-%.2f Gy₁₀/day).",
                        "The DD refit is not: its rate and onset constant reproduce their priors, so DD intervals here",
                        "mostly reflect assumptions. Tumor α/β is not informed by these data (posterior = prior). The trial",
                        "data cover %g-%g Gy per fraction; larger fractions are extrapolations. See Model and sources."),
                  la[1], la[2], sup$d_min, sup$d_max))
    })

    output$context <- renderUI({
      s <- settings()
      div(class = "alert alert-light", style = "padding:0.5rem 0.75rem;",
          sprintf("Evaluating %s against %s. Tumor threshold: a change of at least %g Gy. Early limit: %s. Late limit: %s.%s",
                  regimen_text(s$regimen), regimen_text(s$ref), s$threshold,
                  if (s$early_mode == "hypno") "HYPNO 55 Gy/20 under the same draw" else sprintf("comparison + %g Gy under the same draw", s$tol_early),
                  if (s$late_mode == "hypno") "HYPNO 55 Gy/20" else sprintf("comparison + %g Gy", s$tol_late),
                  if (is.null(oar_fun())) " Organ limits: not applied." else " Organ limits: applied (Organs at risk tab)."))
    })

    output$summary <- renderTable({
      x <- one()$sum; s <- settings()
      n <- x$n_draws
      ess <- attr(draws()$di, "ess")
      rng <- function(a, b) sprintf("%s to %s", fmt_half(a), fmt_half(b))
      crit <- switch(input$tumor_crit, di = x$share_di, dd = x$share_dd, both = x$share_both)
      joint <- mean(switch(input$tumor_crit,
                           di = !is.na(one()$td$delta_di[1, ]) & one()$td$delta_di[1, ] >= s$threshold - 1e-9,
                           dd = !is.na(one()$td$delta_dd[1, ]) & one()$td$delta_dd[1, ] >= s$threshold - 1e-9,
                           both = !is.na(one()$td$delta_di[1, ]) & !is.na(one()$td$delta_dd[1, ]) &
                             one()$td$delta_di[1, ] >= s$threshold - 1e-9 & one()$td$delta_dd[1, ] >= s$threshold - 1e-9) &
                      one()$tox$ok[1, ])
      u <- bed_unit(10)
      data.frame(
        Quantity = c(sprintf("Δ tumor BED, DI, vs %s (%s): median", regimen_text(s$ref), u), "    80% range", "    95% range",
                     sprintf("Δ tumor BED, DD, vs %s (%s): median", regimen_text(s$ref), u), "    80% range", "    95% range",
                     sprintf("Share of draws with Δ ≥ %g Gy, DI", s$threshold), sprintf("Share of draws with Δ ≥ %g Gy, DD", s$threshold),
                     "Share of draws within every toxicity limit",
                     sprintf("Share meeting the tumor criterion (%s) and every limit, draw by draw", input$tumor_crit),
                     "Draws not evaluable (a BED not positive): DI / DD",
                     "Monte Carlo standard error of a share (at most)"),
        Value = c(fmt_half(x$di_median), rng(x$di_q10, x$di_q90), rng(x$di_q025, x$di_q975),
                  fmt_half(x$dd_median), rng(x$dd_q10, x$dd_q90), rng(x$dd_q025, x$dd_q975),
                  fmt_share(x$share_di), fmt_share(x$share_dd), fmt_share(x$share_tox), fmt_share(joint),
                  sprintf("%d / %d of %d", x$n_invalid_di, x$n_invalid_dd, n),
                  sprintf("%.1f percentage points (effective sample size %.0f)", 100 * share_mcse(0.5, min(n, ess)), min(n, ess))),
        check.names = FALSE)
    }, striped = TRUE, bordered = TRUE)

    output$hist <- renderPlot({
      td <- one()$td; s <- settings()
      long <- rbind(data.frame(delta = td$delta_di[1, ], Model = "DI (refit)"),
                    data.frame(delta = td$delta_dd[1, ], Model = "DD (refit, prior-dominated)"))
      long <- long[is.finite(long$delta), ]
      validate(need(nrow(long) > 0, "No evaluable draws."))
      ggplot(long, aes(delta, fill = Model)) +
        geom_histogram(bins = 60, alpha = 0.6, position = "identity", colour = NA) +
        geom_vline(xintercept = 0, colour = "grey40") +
        geom_vline(xintercept = s$threshold, linetype = "dashed") +
        scale_fill_manual(values = c("DI (refit)" = "#0072B2", "DD (refit, prior-dominated)" = "#D55E00")) +
        labs(x = sprintf("Tumor BED of the evaluated regimen minus %s (Gy₁₀)", regimen_text(s$ref)),
             y = "Number of draws") +
        theme_minimal(base_size = 15) + theme(axis.title = element_text(face = "bold"), legend.position = "top")
    })

    output$tcp <- renderTable({
      validate(need(is.numeric(input$tcp0) && input$tcp0 > 0 && input$tcp0 < 100, "Enter a control rate between 1% and 99%."))
      td <- one()$td; d <- draws()
      q <- function(x) stats::quantile(x, c(.1, .5, .9), na.rm = TRUE, names = FALSE)
      di <- q(100 * tcp_given_control(td$delta_di, d$di$alpha, input$tcp0 / 100)[1, ] - input$tcp0)
      dd <- q(100 * tcp_given_control(td$delta_dd, d$dd$alpha, input$tcp0 / 100)[1, ] - input$tcp0)
      data.frame(Model = c("DI (refit)", "DD (refit, prior-dominated)"),
                 `Change in tumor control, median (percentage points)` = sprintf("%+.1f", c(di[2], dd[2])),
                 `80% range` = sprintf("%+.1f to %+.1f", c(di[1], dd[1]), c(di[3], dd[3])), check.names = FALSE)
    }, striped = TRUE, bordered = TRUE)

    result <- reactiveVal(NULL)
    observeEvent(input$run, {
      s <- settings()
      validate(need(is.numeric(input$tox_share) && input$tox_share >= 50 && input$tox_share <= 100, "Required share must be 50-100%."))
      withProgress(message = "Searching under uncertainty", value = 0.2, {
        d <- draws(); tx <- tx_all()
        sub <- list(di = search_subset(d$di), dd = search_subset(d$dd), tx = search_subset(tx))
        res <- uncertain_search(s$ref, sub$di, sub$dd, sub$tx, limits(), s$threshold, input$rule,
                                m_range = s$m_range, d_step = s$d_step, d_max = s$d_max,
                                min_days_saved = s$min_days, calendar = s$calendar,
                                tox_share_min = input$tox_share / 100, oar_fun = oar_fun())
        incProgress(0.5)
        full <- if (nrow(res) > 0) {
          grid <- attr(res, "grid")
          nb <- do.call(rbind, lapply(seq_len(nrow(res)), function(i) {
            g <- grid_neighbors(grid, res$m[i], res$d[i], s$d_step)
            g$listed <- abs(g$d - res$d[i]) < 1e-9
            g
          }))
          regs <- nb[, c("m", "d", "T")]
          td <- tumor_deltas(regs, s$ref, d$di, d$dd)
          tox <- toxicity_ok_draws(regs, s$ref, tx, limits(), if (is.null(oar_fun())) NULL else oar_fun()(regs, tx))
          cbind(summarize_draws(regs, td, tox, s$threshold), listed = nb$listed,
                outside_data = nb$d > tumor_data_support()$d_max + 1e-9)
        } else NULL
        result(list(res = res, full = full, settings = s, rule = input$rule, tox_share = input$tox_share,
                     variant = input$variant, tox_source = input$tox_source))
      })
    })
    output$search_msg <- renderUI({
      r <- result()
      if (is.null(r)) return(helpText("Press 'Run the search' to start (about 2 seconds)."))
      if (nrow(r$res) == 0)
        return(div(class = "alert alert-info", "No regimen in the grid reaches the required share of draws within every toxicity limit."))
      helpText(sprintf("%d numbers of fractions have an eligible regimen. Listed rows are marked; their neighbors (one dose step either side) follow. Summaries use all %d draws.",
                       nrow(r$res), r$full$n_draws[1]))
    })
    output$search_table <- renderTable({
      r <- result(); req(r, !is.null(r$full))
      f <- r$full
      f <- f[order(f$m, f$d), ]
      data.frame(Listed = ifelse(f$listed, "listed", "neighbor"), Fractions = f$m, `Gy per fraction` = sprintf("%.2f", f$d),
                 Days = f$T, `Days saved` = r$settings$ref$T - f$T,
                 `ΔDI median (80% range)` = sprintf("%s (%s to %s)", fmt_half(f$di_median), fmt_half(f$di_q10), fmt_half(f$di_q90)),
                 `ΔDD median (80% range)` = sprintf("%s (%s to %s)", fmt_half(f$dd_median), fmt_half(f$dd_q10), fmt_half(f$dd_q90)),
                 `Within toxicity limits` = fmt_share(f$share_tox),
                 `Both models ≥ threshold and within limits` = fmt_share(f$share_joint),
                 `Outside the tumor data's dose per fraction` = ifelse(f$outside_data, sprintf("yes (> %g Gy)", tumor_data_support()$d_max), ""),
                 check.names = FALSE)
    }, striped = TRUE, bordered = TRUE)

    output$download <- downloadHandler(
      filename = function() sprintf("bed_v%s_uncertainty_%s.zip", APP_VERSION, format(Sys.Date())),
      content = function(file) {
        r <- result()
        s <- settings()
        sc <- uncertainty_scenario(s, input$variant, input$tox_source, if (is.null(r)) NULL else r)
        write_uncertainty_zip(file, sc, one()$sum, if (is.null(r)) NULL else r$full)
      }, contentType = "application/zip")
  })
}
