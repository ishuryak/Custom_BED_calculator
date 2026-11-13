# Version 4.6.2: fixes from the 2026-11-12 external review (B01-B04).
# Record: docs/development/2026_11_12_frequentist_review_response.md
app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
for (f in c("bed_models.R", "schedule.R", "uncertainty.R"))
  source(file.path(app_abs, "R", f), local = TRUE)

# A copy of the installed early-fit package with one manifest field changed (NULL removes it).
mutated_pkg <- function(path, value) {
  d <- tempfile("pkg"); dir.create(d)
  file.copy(file.path(app_abs, "draws", c(EARLY_FIT_MANIFEST, EARLY_FIT_FILE)), d)
  mp <- file.path(d, EARLY_FIT_MANIFEST)
  m <- jsonlite::fromJSON(mp, simplifyVector = FALSE)
  m[[path]] <- value
  writeLines(jsonlite::toJSON(m, auto_unbox = TRUE, digits = NA, null = "null"), mp)
  d
}
mutated_gate <- function(field, value) {
  m <- jsonlite::fromJSON(file.path(app_abs, "draws", EARLY_FIT_MANIFEST), simplifyVector = FALSE)
  g <- m$gate; g[[field]] <- value
  mutated_pkg("gate", g)
}

test_that("B01: the valid package still loads with its point values", {
  f <- load_early_fit(file.path(app_abs, "draws"))
  expect_true(f$ok)
  expect_equal(f$model_id, "early-fit-4747b297ac10")
})

test_that("B01: an unknown evaluator version that only shares a prefix is refused", {
  m <- jsonlite::fromJSON(file.path(app_abs, "draws", EARLY_FIT_MANIFEST), simplifyVector = FALSE)
  f <- load_early_fit(mutated_pkg("evaluator", sub("early-bed-fit-1:", "early-bed-fit-10:", m$evaluator, fixed = TRUE)))
  expect_false(f$ok); expect_match(f$reason, "unsupported evaluator")
  expect_false(load_early_fit(mutated_pkg("evaluator", "early-bed-fit-1x"))$ok)
  expect_true(load_early_fit(mutated_pkg("evaluator", m$evaluator))$ok)       # rewriting it unchanged is fine
})

test_that("B01: mistyped, non-finite or out-of-support diagnostics are refused", {
  bad <- list(list("ess_bulk", "unavailable"), list("ess_tail", "9999"), list("rhat", -1), list("rhat", 0),
              list("rhat", "1.001"), list("ebfmi_min", TRUE), list("divergences", "0"), list("divergences", FALSE),
              list("divergences", 0.5), list("passed", "true"), list("rhat", list(1.001, 1.002)), list("ess_bulk", NULL))
  for (b in bad) {
    f <- load_early_fit(mutated_gate(b[[1]], b[[2]]))
    expect_false(f$ok, info = paste(b[[1]], format(b[[2]])))
    expect_match(f$reason, "diagnostic|do not meet the gate", info = b[[1]])
  }
})

test_that("B01: a missing or malformed model identity is refused", {
  for (v in list(NULL, "", "early-fit-", 12, "early-fit-4747b297ac10 extra")) {
    f <- load_early_fit(mutated_pkg("model_id", v))
    expect_false(f$ok, info = format(v)); expect_match(f$reason, "model identity", info = format(v))
  }
  f <- load_early_fit(mutated_pkg("source_identity", "abc"))
  expect_false(f$ok); expect_match(f$reason, "source identity")
})

test_that("B01: point values must be finite numbers", {
  m <- jsonlite::fromJSON(file.path(app_abs, "draws", EARLY_FIT_MANIFEST), simplifyVector = FALSE)
  pv <- m$point_values; pv$ab_e <- as.character(pv$ab_e)
  expect_false(load_early_fit(mutated_pkg("point_values", pv))$ok)
})

# --- B02: the methods record is checked against the draw files' actual bytes -------------------
source(file.path(app_abs, "R", "methods_page.R"), local = TRUE)
draws_copy <- function() {
  d <- tempfile("draws"); dir.create(d)
  file.copy(list.files(file.path(app_abs, "draws"), full.names = TRUE, pattern = "\\.(json|csv)$"), d)
  d
}

test_that("B02: the unchanged draws give a valid methods record", {
  expect_true(load_model_methods(draws_copy())$ok)
})

test_that("B02: a draw file changed without its manifest makes the record unavailable", {
  for (f in c("tumor_di_primary.csv", "tumor_dd_primary.csv", "early_fit_draws.csv")) {
    d <- draws_copy()
    x <- readLines(file.path(d, f)); x[8] <- sub("^[^,]+", "100", x[8]); writeLines(x, file.path(d, f))
    r <- load_model_methods(d)
    expect_false(r$ok, info = f); expect_match(r$why, "other fitted draws", info = f)
  }
})

test_that("B02: malformed records fail safely instead of throwing", {
  for (mut in list(function(m) { m$tumor <- "x"; m }, function(m) { m$early <- NULL; m },
                   function(m) { m$tumor$di <- list(); m }, function(m) { m$early$draws$sha256 <- 5; m },
                   function(m) { m$tumor$di$priors <- "none"; m })) {
    d <- draws_copy()
    m <- mut(jsonlite::fromJSON(file.path(d, MODEL_METHODS_FILE), simplifyVector = FALSE))
    jsonlite::write_json(m, file.path(d, MODEL_METHODS_FILE), auto_unbox = TRUE, digits = NA)
    r <- expect_no_error(load_model_methods(d))
    expect_false(r$ok)
  }
})

test_that("B02: the release preflight refuses a stale methods record", {
  source(file.path(app_abs, "tools", "deploy_preflight.R"), local = TRUE)
  m <- jsonlite::fromJSON(file.path(app_abs, "draws", "model_methods.json"), simplifyVector = FALSE)
  root <- tempfile("root"); dir.create(root)
  for (f in c(names(m$sources), "draws/model_methods.json", "draws/tumor_di_primary.csv",
              "draws/tumor_dd_primary.csv", "draws/early_fit_draws.csv")) {
    dir.create(dirname(file.path(root, f)), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(app_abs, f), file.path(root, f))
  }
  expect_true(check_methods_record(root))
  x <- readLines(file.path(root, "draws/early_fit_draws.csv")); x[3] <- sub(",[^,]+$", ",9", x[3])
  writeLines(x, file.path(root, "draws/early_fit_draws.csv"))
  expect_error(check_methods_record(root), "out of date for: draws/early_fit_draws.csv")
  cat("# edit\n", file = file.path(root, "R/uncertainty.R"), append = TRUE)
  expect_error(check_methods_record(root), "R/uncertainty.R")
})

# --- B03: correcting HYPNO's description changes the package manifest (its trial-table hash) but not
# the draws; scenarios exported from 4.6.0/4.6.1 must still replay exactly -------------------------
test_that("B03: a 4.6.1 research scenario still replays; an unknown package does not", {
  source(file.path(app_abs, "R", "export.R"), local = TRUE)
  f <- load_early_fit(file.path(app_abs, "draws"))
  p <- fitted_early_params(default_params(), f)
  eb <- early_block(p)
  expect_true(is.list(apply_early_block(default_params(), eb, f)))
  eb$package_sha256 <- "e80ecc45015651cf975d6ae7da43d26080ca7d48df539ecccdca93b4fa1a4355"   # the 4.6.0/4.6.1 manifest
  expect_true(is.list(apply_early_block(default_params(), eb, f)))
  eb$package_sha256 <- strrep("a", 64)
  expect_error(apply_early_block(default_params(), eb, f), "exact replay is impossible")
  eb$package_sha256 <- "e80ecc45015651cf975d6ae7da43d26080ca7d48df539ecccdca93b4fa1a4355"
  eb$draws_sha256 <- strrep("b", 64)                                                     # never with other draws
  expect_error(apply_early_block(default_params(), eb, f), "exact replay is impossible")
})

# --- B04: the tumor equation shows the smooth hinge the fit uses; the early one keeps the exact max ---
test_that("B04: the page shows smax for the tumor model and an exact max for the early model", {
  txt <- paste(readLines(file.path(app_abs, "how_fitted.md"), encoding = "UTF-8"), collapse = "\n")
  tumor_eq <- regmatches(txt, regexpr("\\$\\$\\\\eta_i = [^$]+\\$\\$", txt))
  expect_length(tumor_eq, 1)
  expect_match(tumor_eq, "\\lambda\\,\\mathrm{smax}(T_i - T_k)", fixed = TRUE)
  expect_false(grepl("\\max(0, T_i", tumor_eq, fixed = TRUE))
  expect_match(txt, "\\mathrm{smax}(x) = \\log(1 + e^{3x})/3", fixed = TRUE)
  expect_equal(SMOOTH_K, 3)                                            # the constant the text states
  expect_equal(smooth_max0(0), log(2) / 3, tolerance = 1e-12)          # 0.23 days at onset, as stated
  expect_match(txt, "r\\,\\max(0, T_{\\mathrm{elapsed}} - T_k)", fixed = TRUE)
})
