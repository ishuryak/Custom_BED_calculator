# Mutation check for the version 2.5 fixes. Each mutation restores a defect
# (or a plausible wrong implementation) in a temporary copy of the app and
# runs the suite there; every mutation must turn at least one test red.
# Run from the app folder: Rscript docs/development/2026_10_19_review_scripts/mutations.R
root <- normalizePath(".")
mutations <- list(
  M1_sensitivity_bypass = list("R/bed_models.R", "reason_di = e_di$reason, d_di = e_di$delta",
                               "reason_di = e_di$reason, d_di = a$bed_di[1] - a$bed_di[2]"),
  M2_pass_counts_invalid = list("R/bed_models.R", "n_pass_both = sum(both & st$d_di",
                                "n_pass_both = sum(st$d_di"),
  M3_agreement_wrong_comparator = list("app.R", "    if (identical(input$agree_cmp, \"sidebar\")) { r <- ref(); return(list(m = r$m, d = r$d, T = r$T)) }\n    r <- comparator_presets()$hypno_test",
                                       "    r <- ref()"),
  M4_candidate_leaks_into_agreement = list("app.R", "min_days_saved = checked(\"agree_days\")), tox_args())",
                                           "min_days_saved = checked(\"agree_days\")), search_args()[setdiff(names(search_args()), c(\"m_range\", \"min_days_saved\", \"min_gain\", \"criterion\", \"rank_by\"))])"),
  M5_replay_ignores_hinge = list("R/export.R", "if (!identical(sc$hinge$type, \"softplus\")", "if (FALSE && !identical(sc$hinge$type, \"softplus\")"),
  M6_replay_ignores_equation = list("R/export.R", "if (!identical(as.character(sc$equation_version), EQUATION_VERSION))",
                                    "if (FALSE)"),
  M7_id_key_order = list("R/export.R", "x <- x[order(names(x), method = \"radix\")]", "x <- x"),
  M8_unit_label = list("app.R", "if (isTRUE(all.equal(ab, 10))) \"Gy\\u2081\\u2080\" else", "if (TRUE) \"Gy\\u2081\\u2080\" else"),
  M9_agreement_sidebar_anchor = list("app.R", "agree_anchor <- reactive(anchor_for(agree_ref()))", "agree_anchor <- reactive(anchor())"),
  M10_calibration_not_checked = list("R/export.R", "replay_fail(\"saved calibration C/alpha does not match", "warning(\"saved calibration C/alpha does not match"),
  M11_export_needs_regimen = list("app.R", "regimen_evaluated = tryCatch(regimen(), error = function(e) NULL))",
                                  "regimen_evaluated = regimen())"),
  M12_search_negative_candidate_credited = list("R/bed_models.R", "  cand_ok <- is.finite(cand_bed) & cand_bed > 0",
                                                "  cand_ok <- !is.na(cand_bed)"),
  M13_onset_inf_printed = list("app.R", "ifelse(is.finite(a$dd_onset), sprintf(\"%.1f\", a$dd_onset), \"none\")",
                               "sprintf(\"%.1f\", a$dd_onset)"))
results <- data.frame(mutation = names(mutations), applied = NA, failed = NA)
for (i in seq_along(mutations)) {
  m <- mutations[[i]]
  tmp <- file.path(tempdir(), paste0("mut_", i)); dir.create(tmp)
  for (f in c("app.R", "about.md", "README.md", "DESCRIPTION")) file.copy(file.path(root, f), tmp)
  for (d in c("R", "tests", "References")) file.copy(file.path(root, d), tmp, recursive = TRUE)
  path <- file.path(tmp, m[[1]])
  src <- paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  hit <- grepl(m[[2]], src, fixed = TRUE)
  results$applied[i] <- hit
  if (hit) {
    writeLines(sub(m[[2]], m[[3]], src, fixed = TRUE), path, useBytes = TRUE)
    r <- as.data.frame(testthat::test_dir(file.path(tmp, "tests", "testthat"), reporter = "silent",
                                          stop_on_failure = FALSE))
    results$failed[i] <- sum(r$failed) + sum(r$error)
  }
  unlink(tmp, recursive = TRUE)
  cat(sprintf("%-42s applied=%s failing tests=%s\n", results$mutation[i], results$applied[i], results$failed[i]))
}
stopifnot(all(results$applied), all(results$failed > 0))
cat("All", nrow(results), "mutations caught.\n")
