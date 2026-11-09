# Review of 2026-11-09, items 2.2 and 2.3: the fitted-trial list joins the fitting data and its reference
# sidecar only when both have unique nonempty keys, the same key set and the same PMID for every key, and the
# fitting data matches the SHA-256 in the draws manifest. The list shows each arm's published outcome.
app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
for (f in c("bed_models.R", "schedule.R", "uncertainty.R", "v4_home.R"))
  source(file.path(app_abs, "R", f), local = TRUE)

# A temporary copy of the fitting data and a manifest that names it, optionally with an updated SHA-256.
v45_copy <- function(edit = identity, update_sha = FALSE) {
  d <- tempfile("v45data"); m <- tempfile("v45draws"); dir.create(d); dir.create(m)
  man <- jsonlite::fromJSON(file.path(app_abs, "draws", "draws_manifest.json"), simplifyVector = FALSE)
  x <- utils::read.csv(file.path(app_abs, "data", man$dataset), stringsAsFactors = FALSE,
                       colClasses = c(pmid = "character", source_value_pct = "character"))
  utils::write.csv(edit(x), file.path(d, man$dataset), row.names = FALSE)
  if (update_sha) {
    f <- file.path(d, man$dataset)
    man$data_sha256 <- digest::digest(readBin(f, "raw", n = file.size(f)), algo = "sha256", serialize = FALSE)
  }
  jsonlite::write_json(man, file.path(m, "draws_manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  list(data_dir = d, dir = m)
}
v45_list <- function(edit) {
  ref <- utils::read.csv(file.path(app_abs, V4_FIT_TRIAL_LIST), stringsAsFactors = FALSE,
                         colClasses = c(pmid = "character"))
  tmp <- tempfile(fileext = ".csv"); utils::write.csv(edit(ref), tmp, row.names = FALSE); tmp
}

test_that("a reordered but valid reference list gives the same table (positive control)", {
  withr::with_dir(app_abs, {
    t0 <- v4_fit_trials()
    t1 <- v4_fit_trials(list_path = v45_list(function(r) r[rev(seq_len(nrow(r))), ]))
    expect_identical(t1, t0)
  })
})

test_that("a repeated arm on either side is refused with the arm named", {
  withr::with_dir(app_abs, {
    tmp <- v45_list(function(r) { r[2, ] <- r[1, ]; r })   # same row count, one arm replaced by another
    expect_error(v4_fit_trials(list_path = tmp),
                 "Repeated arm in .*: CAIR_2.0Gy / CAIR, 2.0 Gy/fraction period\\.", class = "v4_fit_trials_error")
    cp <- v45_copy(function(x) { x[2, ] <- x[1, ]; x }, update_sha = TRUE)
    expect_error(v4_fit_trials(dir = cp$dir, data_dir = cp$data_dir),
                 "Repeated arm in the fitting data .*: CAIR_2.0Gy / CAIR, 2.0 Gy/fraction period\\.",
                 class = "v4_fit_trials_error")
  })
})

test_that("an omitted, renamed or empty arm is refused with the arm named", {
  withr::with_dir(app_abs, {
    expect_error(v4_fit_trials(list_path = v45_list(function(r) r[r$ref_no != "Jackson" | !grepl("Twice", r$arm), ])),
                 "does not match the fitting data .*Arms missing from .*: Jackson / Twice daily \\(b\\.i\\.d\\.\\), accelerated\\.",
                 class = "v4_fit_trials_error")
    expect_error(v4_fit_trials(list_path = v45_list(function(r) { r$arm[1] <- "renamed arm"; r })),
                 "Arms not in the fitting data: CAIR_2.0Gy / renamed arm\\.", class = "v4_fit_trials_error")
    expect_error(v4_fit_trials(list_path = v45_list(function(r) { r$ref_no[3] <- " "; r })),
                 "Empty ref_no or arm in row\\(s\\) 3 of ", class = "v4_fit_trials_error")
  })
})

test_that("a substituted PMID is refused, and a reformatted one is accepted", {
  withr::with_dir(app_abs, {
    tmp <- v45_list(function(r) { r$pmid[r$ref_no == "Jackson"][1] <- "9165136"; r })
    expect_error(v4_fit_trials(list_path = tmp),
                 "PMID differs between the fitting data .* for: Jackson / Once daily \\(o\\.d\\.\\), conventional \\(data 9165135, list 9165136\\)\\.",
                 class = "v4_fit_trials_error")
    tmp <- v45_list(function(r) { r$pmid <- paste0("0", r$pmid, " "); r })
    expect_equal(nrow(v4_fit_trials(list_path = tmp)), 45)
  })
})

test_that("an altered source percentage is refused, by the manifest hash and by the published-fitted check", {
  withr::with_dir(app_abs, {
    alter <- function(x) { x$source_value_pct[x$ref_no == "KROG0201" & x$arm == "CONV"] <- "70"; x }
    cp <- v45_copy(alter)
    expect_error(v4_fit_trials(dir = cp$dir, data_dir = cp$data_dir),
                 "differs from the one the draws manifest records \\(SHA-256", class = "v4_fit_trials_error")
    cp <- v45_copy(alter, update_sha = TRUE)
    expect_error(v4_fit_trials(dir = cp$dir, data_dir = cp$data_dir),
                 "Published and fitted outcome disagree .* for: KROG0201 / CONV \\(published 70, fitted 77.8\\)\\.",
                 class = "v4_fit_trials_error")
    # An unaltered copy with its own correct hash is accepted.
    cp <- v45_copy(identity, update_sha = TRUE)
    expect_equal(nrow(v4_fit_trials(dir = cp$dir, data_dir = cp$data_dir)), 45)
  })
})

test_that("all 45 arms join to the right reference and published value", {
  withr::with_dir(app_abs, {
    t <- v4_fit_trials()
    x <- utils::read.csv(file.path("data", draws_manifest()$dataset), stringsAsFactors = FALSE,
                         colClasses = c(pmid = "character", source_value_pct = "character"))
    ref <- utils::read.csv(V4_FIT_TRIAL_LIST, stringsAsFactors = FALSE, colClasses = c(pmid = "character"))
    expect_equal(nrow(t), 45); expect_equal(length(unique(t$pmid)), 19)
    expect_equal(length(unique(t$ref_no)), 21); expect_equal(sum(t$n), 8279)
    k <- function(a, b) paste(a, b, sep = "\r")
    i <- match(k(t$ref_no, t$arm_key), k(x$ref_no, x$arm)); j <- match(k(t$ref_no, t$arm_key), k(ref$ref_no, ref$arm))
    expect_false(anyNA(i)); expect_false(anyNA(j)); expect_equal(sort(i), 1:45)
    expect_identical(t$pmid, x$pmid[i]); expect_identical(t$trial, ref$trial[j]); expect_identical(t$citation, ref$citation[j])
    fail <- x$source_value_is_failure[i] %in% c(TRUE, "True", "TRUE")
    expect_identical(t$published[!fail], x$source_value_pct[i][!fail])
    expect_identical(t$outcome_locator, x$location[i])
    expect_equal(sum(t$figure_read), 4)   # CHART and EORTC 22811, both arms
    expect_setequal(unique(t$ref_no[t$figure_read]), c("CHART", "EORTC22811"))
    expect_equal(sum(t$estimator == "cumulative incidence"), 6)
  })
})

test_that("KROG conventional, both Yamazaki strata and CAIR show the published values", {
  withr::with_dir(app_abs, {
    t <- v4_fit_trials()
    kr <- t[t$ref_no == "KROG0201" & t$arm_key == "CONV", ]
    expect_identical(c(kr$published, kr$horizon, kr$fit_count, kr$fit_rate), c("77.8", "5 y", "64 / 82", "78.0"))
    ya <- t[t$pmid == "16169681", ]
    expect_identical(setNames(ya$published, ya$arm_key), c(`Arm A-1` = "76", `Arm B-1` = "91", `Arm A-2` = "78", `Arm B-2` = "92"))
    expect_identical(sort(unique(ya$ref_no)), c("Yamazaki_large", "Yamazaki_small"))
    ca <- t[t$pmid == "17011446", ]
    expect_identical(ca$published, c("73", "30", "78", "35"))
    expect_identical(ca$fit_rate, c("73.9", "28.6", "78.6", "35.7"))
    expect_true(all(nzchar(ca$schedule_source)) && all(grepl("PMID 10799721", ca$schedule_source)))
  })
})

test_that("trials without an ordinary conventional arm are not given one, and RTOG 88-09 shows its complement", {
  withr::with_dir(app_abs, {
    t <- v4_fit_trials()
    expect_false(any(grepl("control|comparator|reference", names(t), ignore.case = TRUE)))
    for (r in c("CAIR_vs_CB", "RTOG8313_A", "RTOG8809"))
      expect_false(any(grepl("conventional|control|standard", t$arm[t$ref_no == r], ignore.case = TRUE)), info = r)
    expect_identical(t$published[t$ref_no == "RTOG8313_A"], c("25", "37", "42"))
    expect_identical(t$published[t$ref_no == "CAIR_vs_CB"], c("63", "65"))
    rt <- t[t$ref_no == "RTOG8809", ]
    expect_identical(rt$published, c("55 (from 45% failure)", "44 (from 56% failure)"))
    expect_true(all(startsWith(rt$endpoint, "freedom from local-regional failure")))
    expect_equal(sum(t$published_is_failure), 2)
    h <- as.character(v4_fit_trials_table(t))
    expect_match(h, "which also counts patients who had a competing event such as death without failure", fixed = TRUE)
    expect_equal(lengths(regmatches(h, gregexpr("(read from figure)", h, fixed = TRUE))), 4)
    for (k in c("Published outcome (%)", "Endpoint, time point, estimator", "Count in fit",
                "Schedule: Cummings et al. Radiother Oncol 2007 (PMID 17920715), Table 2, p.12"))
      expect_match(h, k, fixed = TRUE)
  })
})

test_that("the dialog says the list is unavailable instead of failing or showing an unchecked citation", {
  withr::with_dir(app_abs, {
    ok <- v4_fit_trials_dialog("refit_explainer.html")
    expect_true(ok$ok); expect_match(ok$title, "(19 trials, 45 arms or strata)", fixed = TRUE)
    bad <- v4_fit_trials_dialog("refit_explainer.html", list_path = v45_list(function(r) r[-1, ]))
    expect_false(bad$ok); expect_match(bad$title, "list unavailable", fixed = TRUE)
    h <- as.character(bad$body)
    expect_match(h, "Arms missing from", fixed = TRUE)
    expect_false(grepl("pubmed.ncbi.nlm.nih.gov", h, fixed = TRUE))
    expect_false(v4_fit_trials_dialog("x", list_path = "no_such_file.csv")$ok)
  })
})

# Review of 2026-11-09, section 3: corrections to the long refit explainer, with counts derived from the data.
test_that("the refit explainer states counts from the data and the corrected model statements", {
  flat <- function(f) gsub("\\s+", " ", paste(readLines(f, encoding = "UTF-8", warn = FALSE), collapse = "\n"))
  md <- flat(file.path(app_abs, "docs", "refit_explainer", "HNSCC_BED_v4_refit_explainer_2026_11_04.md"))
  html <- flat(file.path(app_abs, "www", "refit_explainer.html"))
  x <- utils::read.csv(file.path(app_abs, "data", "refit_trials_2026_10_25c_main.csv"), stringsAsFactors = FALSE)
  y <- round(x$n_patients * x$tcp_actual_pct / 100); fit <- 100 * y / x$n_patients
  n_diff <- sum(round(fit, 1) != round(x$tcp_actual_pct, 1)); max_diff <- max(abs(fit - x$tcp_actual_pct))
  df <- nrow(x) - length(unique(x$ref_no))
  want <- c("fit generation 2026-10-25",
            sprintf("%d within-unit contrast degrees of freedom", df),
            sprintf("%d of the %d arms, by up to %.2f percentage points", n_diff, nrow(x), max_diff),
            sprintf("8,279 analyzed patients in %d fitting units", length(unique(x$ref_no))),
            "35 × 2.1 Gy in 47 days gains 4.94 Gy over 35 × 2 Gy in 47 days under DI and 3.02 Gy under DD",
            "to −7.98 Gy under both", "0.53 Gy EQD2 per day (0.18–0.99)", "which is always shown",
            "the last ten received 62.8 Gy", "Table 1, p.102")
  for (k in want) { expect_match(md, k, fixed = TRUE, info = k); expect_match(html, k, fixed = TRUE, info = k) }
  expect_match(md, "| Published (%) | Count in fit (%) |", fixed = TRUE)
  expect_equal(sum(x$n_patients), 8279)
  for (k in c("It does not change BED differences", "version 4.3", "Observed (%)", "Under Details on the main page",
              "differ only in how they extrapolate", "About 24 within-trial contrasts", "Every value has a verbatim quote"))
    expect_false(grepl(k, md, fixed = TRUE) || grepl(k, html, fixed = TRUE), info = k)
})
