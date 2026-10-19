# Starts the calculator with the release library when it is restored
# (tools/restore_library.R), otherwise with the packages installed in R.
source("tools/env.R")
if (!use_bed_library()) cat("Note: release library not restored; using the packages installed in R.\n")
environment_report()
shiny::runApp(".", launch.browser = TRUE)
