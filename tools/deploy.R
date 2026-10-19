# Deploys version 2 to shinyapps.io as a NEW app, "custom_bed_calculator_v2",
# so the version 1 site (custom_bed_calculator) is never overwritten.
# Run from the app folder, after tools/restore_library.R and after the
# shinyapps.io account has been set once with rsconnect::setAccountInfo():
#   Rscript tools/deploy.R
# rsconnect lives in a separate tools library (BED_DEPLOY_LIBRARY), so the
# release library holds only what the app and tests use.
source("tools/env.R")
APP_NAME <- "custom_bed_calculator_v2"
stopifnot(APP_NAME != "custom_bed_calculator")
deploy_lib <- Sys.getenv("BED_DEPLOY_LIBRARY",
                         file.path(tools::R_user_dir("HNSCC_BED_calculator", "cache"), "deploy-tools"))
invisible(use_bed_library(required = TRUE))
.libPaths(c(.libPaths()[1], deploy_lib, .Library))
if (!requireNamespace("rsconnect", quietly = TRUE)) {
  dir.create(deploy_lib, recursive = TRUE, showWarnings = FALSE)
  install.packages("rsconnect", lib = deploy_lib, repos = LOCK_REPO,
                   type = if (.Platform$OS.type == "windows") "win.binary" else getOption("pkgType"))
}
acc <- rsconnect::accounts()
acc <- acc[acc$server == "shinyapps.io", , drop = FALSE]
if (nrow(acc) == 0)
  stop("No shinyapps.io account is set. In R, run rsconnect::setAccountInfo(name, token, secret) ",
       "with the values from shinyapps.io > Account > Tokens, then rerun this script.")
environment_report()
cat("Files:", rsconnect::listDeploymentFiles("."), "\n")
rsconnect::deployApp(appDir = ".", appName = APP_NAME,
                     appTitle = "HNSCC BED calculator, version 2",
                     account = acc$name[1], server = "shinyapps.io",
                     launch.browser = FALSE, forceUpdate = TRUE)
