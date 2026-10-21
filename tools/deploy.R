# Deploys version 2 to shinyapps.io as a NEW app, "custom_bed_calculator_v2",
# so the version 1 site (custom_bed_calculator) is never overwritten.
# Run from the app folder, after tools/restore_library.R and after the
# shinyapps.io account has been set once with rsconnect::setAccountInfo():
#   Rscript tools/deploy.R
# rsconnect lives in a separate tools library (BED_DEPLOY_LIBRARY), so the
# release library holds only what the app and tests use.
source("tools/env.R")
source("tools/deploy_preflight.R")
APP_NAME <- "custom_bed_calculator_v2"
ACCOUNT <- "ishuryak"  # pinned: never deploy to whichever account happens to be first
deploy_lib <- Sys.getenv("BED_DEPLOY_LIBRARY",
                         file.path(tools::R_user_dir("HNSCC_BED_calculator", "cache"), "deploy-tools"))
invisible(use_bed_library(required = TRUE))
.libPaths(c(.libPaths()[1], deploy_lib, .Library))
if (!requireNamespace("rsconnect", quietly = TRUE)) {
  dir.create(deploy_lib, recursive = TRUE, showWarnings = FALSE)
  install.packages("rsconnect", lib = deploy_lib, repos = LOCK_REPO,
                   type = if (.Platform$OS.type == "windows") "win.binary" else getOption("pkgType"))
}
dirty <- system("git status --porcelain --untracked-files=no", intern = TRUE)
deploy_preflight(APP_NAME, "custom_bed_calculator_v2", rsconnect::accounts(), ACCOUNT, dirty)
cat("Deploying revision", system("git rev-parse HEAD", intern = TRUE), "\n")
environment_report()
cat("Files:", rsconnect::listDeploymentFiles("."), "\n")
rsconnect::deployApp(appDir = ".", appName = APP_NAME,
                     appTitle = "HNSCC BED calculator, version 2",
                     account = ACCOUNT, server = "shinyapps.io",
                     launch.browser = FALSE, forceUpdate = TRUE)
