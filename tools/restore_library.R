# Restores the package versions in renv.lock into the release library (see
# tools/env.R). Run from the app folder:  Rscript tools/restore_library.R
# Binaries come from the dated Posit Package Manager snapshot named in
# renv.lock (Windows and macOS). On Linux set RENV_CONFIG_REPOS_OVERRIDE to
# the matching __linux__ URL for binaries (see .github/workflows/tests.yml).
source("tools/env.R")
lib <- bed_library()
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
repo <- Sys.getenv("RENV_CONFIG_REPOS_OVERRIDE", LOCK_REPO)
.libPaths(c(lib, .Library))
if (!requireNamespace("renv", quietly = TRUE))
  install.packages("renv", lib = lib, repos = repo,
                   type = if (.Platform$OS.type == "windows") "win.binary" else getOption("pkgType"))
renv::restore(lockfile = "renv.lock", library = lib, prompt = FALSE)
cat("Release library:", lib, "\n")
environment_report()
