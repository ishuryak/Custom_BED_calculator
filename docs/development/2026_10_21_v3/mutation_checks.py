"""Deliberate-break checks for the version 3 guards. Each mutation is applied to
a temporary copy of the repository; the named test file must then FAIL.
Run from the repository root:  python3 docs/development/2026_10_21_v3/mutation_checks.py"""
import os, shutil, subprocess, sys, tempfile
RS = "/mnt/c/Program Files/R/R-4.5.1/bin/Rscript.exe"
ROOT = os.getcwd()
TMP = "/mnt/c/Users/ishur/AppData/Local/Temp/bed_v3_mutations"
M = [
 ("M1 draws attributes dropped (the original bug)", "R/uncertainty.R", "  attr(x, \"ess\") <- man$files[[key]]$min_ess_bulk\n  x\n}", "  attr(x, \"ess\") <- man$files[[key]]$min_ess_bulk\n  x[names(x)]\n}", "test_v3_uncertainty.R"),
 ("M2 identifier check disabled", "R/oar.R", "identifier_lines <- function(lines) which(grepl(ID_PATTERN, lines, ignore.case = TRUE, perl = TRUE))", "identifier_lines <- function(lines) integer()", "test_v3_oar_schedule.R"),
 ("M3 invalid draws count as meeting", "R/uncertainty.R", "  ok <- function(b, r) is.finite(b) & b > 0 & matrix(", "  ok <- function(b, r) TRUE | matrix(", "test_v3_uncertainty.R"),
 ("M4 conversion to the wrong fraction number", "R/oar.R", "  physical_from_eqd2(constraint_eqd2(limit_gy, n_fractions, ab), m, ab)", "  physical_from_eqd2(constraint_eqd2(limit_gy, n_fractions, ab), m + 1, ab)", "test_v3_oar_schedule.R"),
 ("M5 comparison not on the same draws", "R/uncertainty.R", "  r_di <- bed_di_draws(ref$m, ref$d, ref$T, di)[1, ]", "  r_di <- rev(bed_di_draws(ref$m, ref$d, ref$T, di)[1, ])", "test_v3_uncertainty.R"),
 ("M6 organ check ignored in the search", "R/bed_models.R", "oar_check(m, d)", "rep(TRUE, length(m))", "test_v3_app.R"),
 ("M7 fitted DD ignored", "R/bed_models.R", "  if (identical(p$dd_mode, \"fitted\")) {", "  if (FALSE) {", "test_v3_app.R"),
 ("M8 both-model screen ignores the global threshold", "app.R", "agree_threshold <- reactive(if (identical(input$agree_threshold_mode, \"own\")) checked(\"agree_margin\") else checked(\"min_gain\"))", "agree_threshold <- reactive(checked(\"agree_margin\"))", "test_v3_app.R"),
 ("M9 DVH bins at their lower edge", "R/oar.R", "  mid <- c((dose[-k] + dose[-1]) / 2, dose[k])", "  mid <- c(dose[-k], dose[k])", "test_v3_oar_schedule.R"),
 ("M10 skipped days ignored", "R/schedule.R", "if (wd <= fx_per_week && !(day %in% skip_days))", "if (wd <= fx_per_week)", "test_v3_oar_schedule.R"),
 ("M11 bin-by-bin mean replaced by mean-only conversion", "R/oar.R", "       mean_eqd2 = sum(w * eqd2_physical(diff_dvh$dose, n, ab)),", "       mean_eqd2 = eqd2_physical(sum(w * diff_dvh$dose), n, ab),", "test_v3_oar_schedule.R"),
 ("M12 verdict tolerance removed", "app.R", "VERDICT_TOL <- 0.5", "VERDICT_TOL <- 0", "test_v3_app.R"),
]
log = []
for name, f, old, new, test in M:
    if os.path.exists(TMP): shutil.rmtree(TMP)
    shutil.copytree(ROOT, TMP, ignore=shutil.ignore_patterns(".git", "fit", "docs", "References"))
    src = open(os.path.join(TMP, f), encoding="utf-8").read()
    if src.count(old) != 1:
        log.append(f"SETUP-ERROR {name}: pattern found {src.count(old)} times"); continue
    open(os.path.join(TMP, f), "w", encoding="utf-8").write(src.replace(old, new))
    win = "C:\\Users\\ishur\\AppData\\Local\\Temp\\bed_v3_mutations"
    r = subprocess.run([RS, "-e", f'setwd("{win.replace(chr(92), "/")}"); source("tools/env.R"); invisible(use_bed_library(TRUE)); '
                        f'r <- as.data.frame(testthat::test_file("tests/testthat/{test}", reporter = "silent")); '
                        f'cat("FAILED", sum(r$failed) + sum(r$error), "\\n")'], capture_output=True, text=True)
    out = [l for l in r.stdout.splitlines() if l.startswith("FAILED")]
    nfail = int(out[0].split()[1]) if out else -1
    log.append(f"{'CAUGHT' if nfail > 0 else 'MISSED'}  {name}  ({test}: {nfail} failing)")
    print(log[-1], flush=True)
shutil.rmtree(TMP, ignore_errors=True)
open(os.path.join(ROOT, "docs/development/2026_10_21_v3/mutation_checks.log"), "w").write("\n".join(log) + "\n")
sys.exit(0 if all(l.startswith("CAUGHT") for l in log) else 1)
