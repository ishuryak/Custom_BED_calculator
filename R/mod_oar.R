# Organs-at-risk tab (Shiny module). Requires R/oar.R.
# settings(): list(regimen = list(m, d, T), ref, label_ref).
# Returns a reactive list: table (feasibility for the evaluated regimen),
# percents (organ -> % of prescription, or NULL), search_check (TRUE when the
# entered percentages should also limit the searches), library, ab.

oar_ui <- function(id) {
  ns <- NS(id)
  organ_inputs <- lapply(seq_len(nrow(OAR_ORGANS)), function(i) {
    k <- OAR_ORGANS$key[i]
    div(class = "oar-organ", style = "display:flex; gap:0.5rem; align-items:end; flex-wrap:wrap;",
        div(style = "min-width:15rem;", tags$label(class = "form-label", `for` = ns(paste0("v_", k)),
                                                   sprintf("%s (%s)", OAR_ORGANS$label[i],
                                                           if (OAR_ORGANS$metric[i] == "Dmax") "maximum" else "mean"))),
        div(style = "width:9rem;", numericInput(ns(paste0("v_", k)), NULL, value = NA, min = 0, max = 500, step = 1)),
        conditionalPanel(sprintf("input['%s'] == 'dvh'", ns("mode")),
          div(style = "width:14rem;", selectInput(ns(paste0("s_", k)), NULL, choices = c("(not assessed)" = "")))))
  })
  tagList(
    div(class = "alert alert-secondary",
        "Maximum-dose constraints from conventional fractionation (QUANTEC) or from the HYHOPE protocol's",
        "hypofractionated tables are converted to the evaluated regimen's number of fractions with the linear-quadratic",
        "model and an organ α/β, the method the HYPNO trial describes. Mean-dose limits are physical means from their",
        "source: they are compared with the planned physical mean, and the mean EQD2 is shown for illustration only.",
        "A status is a statement about the source limit, not about safety. Organs left blank are not assessed. These are",
        "selected constraint metrics, not full protocol compliance. No normal-tissue complication probability is calculated."),
    layout_columns(col_widths = c(4, 4, 4),
      radioButtons(ns("library"), "Constraint source",
                   c("QUANTEC 2010 (conventional fractionation: maxima converted, means physical)" = "QUANTEC",
                     "HYHOPE protocol tables (12, 15, 20 fractions; nearest converted)" = "HYHOPE")),
      div(numericInput(ns("ab_cord"), "α/β, spinal cord and brainstem (Gy)", 2, 0.5, 10, 0.1),
          helpText("HYPNO used 2 Gy; QUANTEC (Kirkpatrick 2010) gives 0.87 Gy. The status is also checked from 0.87 to 2 Gy.")),
      div(numericInput(ns("ab_oar"), "α/β, other organs (Gy)", 3, 0.5, 20, 0.5),
          helpText("HYPNO used 3 Gy. The status is also checked from 2 to 4.5 Gy."))),
    radioButtons(ns("mode"), "Enter organ doses as", inline = TRUE,
                 c("% of the prescription dose (from a reference plan)" = "percent",
                   "Gy, in the evaluated regimen" = "gy",
                   "a dose-volume histogram (DVH) file" = "dvh")),
    conditionalPanel(sprintf("input['%s'] == 'percent'", ns("mode")),
      helpText("Assumes the plan scales linearly with the prescription: an organ at 70% of the prescription receives",
               "70% of every fraction. Suited to a planning study, not to a specific patient's plan.")),
    conditionalPanel(sprintf("input['%s'] == 'dvh'", ns("mode")),
      div(class = "alert alert-warning",
          tags$strong("Do not upload identifiable patient data to the public website. "),
          "Files containing patient name or ID fields with values are refused. The upload is saved as a temporary file on the server,",
          "read into memory and deleted at once. It is not kept, logged or included in downloads. For real patient plans",
          "use the local version (run_app.bat).",
          "Use of real plans in a study may need institutional review board approval."),
      fileInput(ns("dvh_file"), "Cumulative DVH: Eclipse-style text export, or a CSV with dose and volume columns (max 2 MB)",
                accept = c(".txt", ".csv", ".dvh")),
      layout_columns(col_widths = c(6, 6),
        numericInput(ns("dvh_total"), "Prescription total dose of the plan the DVH comes from (Gy)", 70, 1, 200, 0.5),
        numericInput(ns("dvh_n"), "Number of fractions of that plan", 35, 1, 100, 1)),
      helpText("DVH doses are scaled to the evaluated regimen by its total dose divided by the plan's, and treated as",
               "given in the evaluated regimen's number of fractions. Mean doses are converted bin by bin."),
      uiOutput(ns("dvh_msg"))),
    h5("Organ doses"),
    helpText(textOutput(ns("unit_hint"), inline = TRUE)),
    do.call(tagList, organ_inputs),
    checkboxInput(ns("use_in_search"),
                  "Also apply the entered organ limits in the Search and Uncertainty tabs (% of prescription only)", FALSE),
    uiOutput(ns("context")),
    h5("Feasibility for the evaluated regimen"),
    div(style = "overflow-x:auto;", tableOutput(ns("table"))),
    uiOutput(ns("binding")),
    h5("The same organ doses under published hypofractionated regimens"),
    helpText("Only for doses entered as % of the prescription."),
    div(style = "overflow-x:auto;", tableOutput(ns("bench"))),
    downloadButton(ns("download"), "Download the feasibility table (CSV, protocol-style)"))
}

oar_server <- function(id, settings) {
  moduleServer(id, function(input, output, session) {
    lib <- load_constraints()
    dvh <- reactiveVal(NULL)
    dvh_error <- reactiveVal(NULL)

    observeEvent(input$dvh_file, {
      f <- input$dvh_file
      res <- tryCatch(read_dvh_upload(f$datapath, f$size), error = function(e) e)
      unlink(f$datapath)  # delete Shiny's temporary upload file at once
      if (inherits(res, "error")) { dvh(NULL); dvh_error(conditionMessage(res)); return() }
      dvh(res); dvh_error(NULL)
      ch <- c("(not assessed)" = "", stats::setNames(names(res), names(res)))
      for (k in OAR_ORGANS$key) updateSelectInput(session, paste0("s_", k), choices = ch)
    })
    output$dvh_msg <- renderUI({
      if (!is.null(dvh_error())) return(div(class = "alert alert-danger", dvh_error()))
      if (!is.null(dvh())) helpText(sprintf("Read %d structure%s: %s. Assign them to organs below.", length(dvh()),
                                            if (length(dvh()) > 1) "s" else "", paste(names(dvh()), collapse = ", ")))
    })
    output$unit_hint <- renderText(switch(input$mode,
      percent = "Enter each organ's maximum or mean dose as a percentage of the prescription total dose.",
      gy = "Enter each organ's physical maximum or mean dose in Gy, for the evaluated regimen.",
      dvh = "Choose the DVH structure for each organ."))

    ab <- reactive({
      validate(need(is.numeric(input$ab_cord) && is.finite(input$ab_cord) && input$ab_cord > 0, "Enter the cord α/β."),
               need(is.numeric(input$ab_oar) && is.finite(input$ab_oar) && input$ab_oar > 0, "Enter the organ α/β."))
      c(cord = input$ab_cord, oar = input$ab_oar)
    })

    inputs_for <- function(m, d) {
      out <- list()
      for (k in OAR_ORGANS$key) {
        v <- input[[paste0("v_", k)]]
        if (input$mode == "dvh") {
          s <- input[[paste0("s_", k)]]
          if (is.null(dvh()) || is.null(s) || !nzchar(s) || is.null(dvh()[[s]])) next
          validate(need(is.numeric(input$dvh_total) && input$dvh_total > 0, "Enter the plan's prescription dose."))
          x <- dvh()[[s]]
          scale <- m * d / input$dvh_total
          out[[k]] <- list(mode = "dvh", value = data.frame(dose = x$dose * scale, volume = x$volume))
        } else if (is.numeric(v) && length(v) == 1 && is.finite(v)) {
          validate(need(v >= 0 && v <= if (input$mode == "percent") 200 else 500,
                        "Organ doses must be between 0 and 200% (or 500 Gy)."))
          out[[k]] <- list(mode = input$mode, value = v)
        }
      }
      out
    }

    table_for <- function(r) oar_table(lib, inputs_for(r$m, r$d), r$m, r$d, input$library, ab(),
                                       ab_ranges = list(cord = c(0.87, 2), oar = c(2, 4.5)))
    tab <- reactive(table_for(settings()$regimen))

    pretty <- function(t, m) {
      num <- function(x) ifelse(is.na(x), "", sprintf("%.1f", x))
      mean_m <- t$metric == "Dmean"
      planned_phys <- ifelse(!is.na(t$planned_max_lo), sprintf("%.1f to %.1f (DVH bin)", t$planned_max_lo, t$planned_max_hi),
                             num(t$planned_physical))
      data.frame(Organ = t$organ, Metric = ifelse(t$metric == "Dmax", "maximum", "mean"),
                 `Source constraint` = ifelse(is.na(t$source), "none verified", t$source),
                 `Source metric and fractionation` = ifelse(is.na(t$source_metric), "",
                                                            sprintf("%s; %s", t$source_metric, t$source_fractionation)),
                 `Compared as` = ifelse(is.na(t$basis), "", t$basis),
                 `α/β (Gy)` = sprintf("%g", t$ab),
                 `Limit, EQD2 (Gy)` = num(t$limit_eqd2),
                 `Limit, physical in this regimen (Gy)` = ifelse(is.na(t$limit_physical), "", sprintf("%s%.1f", t$operator, t$limit_physical)),
                 `Planned, physical (Gy)` = planned_phys,
                 `Planned, EQD2 (Gy)` = ifelse(mean_m & !is.na(t$planned_eqd2), paste(num(t$planned_eqd2), "(illustration)"), num(t$planned_eqd2)),
                 `Margin (Gy, in the compared quantity)` = num(t$margin_gy),
                 Status = t$status, `Over the α/β range` = t$flip, Notes = t$flags,
                 check.names = FALSE)
    }
    coverage <- reactive({
      t <- tab()
      oar_coverage(lib, stats::setNames(as.list(ifelse(t$entered, 1, NA)), OAR_ORGANS$key), input$library)
    })
    output$context <- renderUI({
      r <- settings()$regimen
      div(class = "alert alert-light", style = "padding:0.5rem 0.75rem;",
          sprintf("Evaluated regimen: %d × %g Gy (%g Gy) in %g days. Constraint source: %s. %s α/β: %g Gy (cord, brainstem), %g Gy (other organs), the HYPNO assumptions.",
                  r$m, r$d, r$m * r$d, r$T, OAR_LIBRARY_META[[input$library]]$version,
                  OAR_LIBRARY_META[[input$library]]$scope, ab()[["cord"]], ab()[["oar"]]))
    })
    output$table <- renderTable(pretty(tab(), settings()$regimen$m), striped = TRUE, bordered = TRUE)
    output$binding <- renderUI({
      t <- tab()
      cov <- helpText(oar_coverage_text(coverage()))
      unver <- t$organ[t$entered & is.na(t$source)]
      unver_txt <- if (length(unver)) helpText(sprintf("Not assessed: %s %s no verified constraint in this library.",
                                                       paste(unver, collapse = ", "), if (length(unver) > 1) "have" else "has"))
      if (!any(t$binds_first)) {
        return(tagList(helpText(if (any(t$entered)) "No entered organ has a verified constraint in this library, so nothing is assessed." else
          "No organ dose is entered, so nothing is assessed."), unver_txt, cov))
      }
      b <- t[t$binds_first, ]
      lim <- if (b$metric == "Dmean") b$limit_physical else b$limit_eqd2
      tagList(helpText(sprintf("Binds first: %s, with a margin of %.1f Gy (%s, %.0f%% of its limit).",
                               b$organ, b$margin_gy, if (b$metric == "Dmean") "physical mean" else sprintf("EQD2 under α/β %g Gy", b$ab),
                               100 * b$margin_gy / lim)), unver_txt, cov)
    })
    output$bench <- renderTable({
      validate(need(input$mode == "percent", "Enter organ doses as % of the prescription to compare regimens."))
      bm <- benchmark_regimens()
      rows <- lapply(seq_len(nrow(bm)), function(i) {
        t <- table_for(list(m = bm$m[i], d = bm$d[i], T = bm$T[i]))
        t <- t[t$entered & !is.na(t$source), ]
        if (nrow(t) == 0) return(NULL)
        mean_m <- t$metric == "Dmean"
        data.frame(Regimen = bm$name[i], Organ = t$organ,
                   `Compared as` = ifelse(mean_m, "physical mean", "EQD2"),
                   `Planned (Gy)` = sprintf("%.1f", ifelse(mean_m, t$planned_physical, t$planned_eqd2)),
                   `Limit (Gy)` = sprintf("%s%.1f", t$operator, ifelse(mean_m, t$limit_physical, t$limit_eqd2)),
                   Status = t$status, check.names = FALSE)
      })
      do.call(rbind, rows)
    }, striped = TRUE, bordered = TRUE)

    output$download <- downloadHandler(
      filename = function() sprintf("oar_feasibility_v%s_%s.csv", APP_VERSION, format(Sys.Date())),
      content = function(file) {
        r <- settings()$regimen
        t <- tab()
        hdr <- c(sprintf("# HNSCC BED calculator %s, organ-at-risk feasibility. Regimen %d x %g Gy in %g days.",
                         APP_VERSION, r$m, r$d, r$T),
                 sprintf("# Library %s; alpha/beta cord %g Gy, other %g Gy. Organ rules %s. Maxima are compared as EQD2; mean-dose limits as physical means (EQD2 for illustration only). Status is not a statement of safety.",
                         OAR_LIBRARY_META[[input$library]]$version, ab()[["cord"]], ab()[["oar"]], OAR_VERSION),
                 paste("#", oar_coverage_text(coverage())),
                 "# Uploaded DVH content is not included in this file.")
        writeLines(hdr, file)
        suppressWarnings(utils::write.table(t, file, sep = ",", row.names = FALSE, append = TRUE))
      })

    reactive(list(
      table = tab(), library = input$library, ab = ab(),
      percents = if (input$mode == "percent") {
        p <- lapply(stats::setNames(OAR_ORGANS$key, OAR_ORGANS$key), function(k) input[[paste0("v_", k)]])
        p[vapply(p, function(v) is.numeric(v) && length(v) == 1 && is.finite(v), logical(1))]
      } else NULL,
      search_check = isTRUE(input$use_in_search) && identical(input$mode, "percent")))
  })
}
