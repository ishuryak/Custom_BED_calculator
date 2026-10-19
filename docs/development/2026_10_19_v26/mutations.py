# Mutation check for version 2.6: restore each old line (one at a time) and
# confirm tests/testthat/test_v26_review_fixes.R turns red. Run from the app
# folder: python3 docs/development/2026_10_19_v26/mutations.py
import subprocess, shutil, sys
RS = "/mnt/c/Program Files/R/R-4.5.1/bin/Rscript.exe"
M = [
 ("M1 conservative score ignores missing DD", "R/bed_models.R",
  "worst = ifelse(ok_di & ok_dd, pmin(a$d_di, a$d_dd), NA_real_),", "worst = pmin(a$d_di, a$d_dd, na.rm = TRUE),"),
 ("M2 unrankable rows listed", "R/bed_models.R",
  "  a <- a[a$rank_eligible, , drop = FALSE]\n", "  a <- a\n"),
 ("M3 per-model range masks", "R/bed_models.R",
  "di <- rng(st$d_di, both); dd <- rng(st$d_dd, both)", "di <- rng(st$d_di, st$valid_di); dd <- rng(st$d_dd, st$valid_dd)"),
 ("M4 d_min from one model", "R/bed_models.R",
  "a$d_min <- ifelse(is.finite(a$d_di) & is.finite(a$d_dd), pmin(a$d_di, a$d_dd), NA_real_)", "a$d_min <- pmin(a$d_di, a$d_dd, na.rm = TRUE)"),
 ("M5 next step ignores ranking", "R/bed_models.R",
  'if (length(parts) == 0 && !r) return("ranking not evaluable (a tumor BED is not positive)")', "NULL"),
 ("M6 fenwick coerced with isTRUE", "R/export.R",
  'fenwick = check_flag(s$fenwick, "fenwick"))', "fenwick = isTRUE(s$fenwick))"),
 ("M7 agreement contradictions accepted", "R/export.R",
  'if (!identical(args$criterion, "both") || !identical(args$rank_by, "worst") ||', "if (FALSE &&"),
 ("M8 onset_anchor unchecked", "R/export.R",
  "if (is.finite(cal$onset_ref) && !is.null(saved$onset_anchor)", "if (FALSE"),
 ("M9 old selection replayed silently", "R/export.R",
  'if (schema != SCHEMA_VERSION || !identical(as.character(sc$selection_version), SELECTION_VERSION)) {', "if (FALSE) {"),
 ("M10 grid size unbounded", "R/export.R",
  "if (n_d * length(mr) > MAX_GRID_ROWS)", "if (FALSE)"),
 ("M11 hash via digest(file=)", "R/export.R",
  'digest::digest(bytes, algo = "sha256", serialize = FALSE)', 'digest::digest(bytes, algo = "md5", serialize = FALSE)'),
 ("M12 code identity taken at export", "R/export.R",
  "  if (is.null(code)) code <- now\n", "  code <- now\n"),
 ("M13 old sidebar help restored", "app.R",
  '"Used by the Results, Candidate regimens and Sensitivity tabs: their differences, and",', '"All differences and toxicity limits are relative to this regimen.",'),
 ("M14 version footer hard-coded", "app.R",
  'txt <- gsub("{{VERSIONS}}", sprintf("Version %s of', 'txt <- gsub("{{VERSIONS}}", sprintf("Version 2.5 %s of'),
]
res = []
for name, f, old, new in M:
    src = open(f, encoding="utf-8").read()
    assert src.count(old) == 1, (name, src.count(old))
    shutil.copy(f, f + ".bak")
    open(f, "w", encoding="utf-8").write(src.replace(old, new))
    try:
        out = subprocess.run([RS, "-e", 'r <- as.data.frame(testthat::test_file("tests/testthat/test_v26_review_fixes.R", reporter = "silent")); cat("NFAIL", sum(r$failed) + sum(r$error), "\\n")'],
                             capture_output=True, text=True, timeout=900).stdout
    finally:
        shutil.move(f + ".bak", f)
    n = int(out.split("NFAIL")[1].split()[0]) if "NFAIL" in out else -1
    res.append((name, n)); print(name, "->", "CAUGHT" if n > 0 else "SURVIVED", n, flush=True)
sys.exit(0 if all(n > 0 for _, n in res) else 1)
