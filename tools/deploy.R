# Deploys version 3 to shinyapps.io as its own app, "custom_bed_calculator_v3".
# The version 1 site (custom_bed_calculator) and the version 2 site
# (custom_bed_calculator_v2) are never touched: deploy_preflight() refuses them.
# Run from the app folder, after tools/restore_library.R and a green
# tools/run_tests.R on the same commit:
#   Rscript tools/deploy.R
# rsconnect lives in a separate tools library (BED_DEPLOY_LIBRARY), so the
# release library holds only what the app and tests use.
source("tools/env.R")
source("tools/deploy_preflight.R")
APP_NAME <- "custom_bed_calculator_v3"
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
# Checked git calls only (review M01): stops before anything is written or uploaded.
id <- deploy_identity()
head <- id$rev; passed <- id$passed
deploy_preflight(APP_NAME, "custom_bed_calculator_v3", rsconnect::accounts(), ACCOUNT, character(),
                 head = head, passed = passed)
writeLines(head, "deployed_revision.txt")   # shown on the Model and sources page
cat("Deploying revision", head, "\n")
environment_report()
files <- rsconnect::listDeploymentFiles(".")
check_bundle(files, tracked = id$tracked)
cat("Files:", files, "\n")
rsconnect::deployApp(appDir = ".", appName = APP_NAME,
                     appTitle = "HNSCC BED calculator, version 3",
                     account = ACCOUNT, server = "shinyapps.io",
                     launch.browser = FALSE, forceUpdate = TRUE)
