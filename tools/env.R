# The release package library: the versions in renv.lock, restored by
# tools/restore_library.R into a per-user folder outside this project (so it
# is not synced by OneDrive or deployed). BED_LIBRARY overrides the location.
LOCK_REPO <- local({
  x <- readLines("renv.lock", warn = FALSE)
  sub('.*"URL": "([^"]+)".*', "\\1", grep('"URL":', x, value = TRUE)[1])
})
bed_library <- function() {
  Sys.getenv("BED_LIBRARY", file.path(tools::R_user_dir("HNSCC_BED_calculator", "cache"),
                                      paste0("library-", basename(LOCK_REPO), "-R", getRversion())))
}
# Put the release library first and drop the user library, so nothing else
# installed on this computer is used. Returns FALSE if it is not restored.
use_bed_library <- function(required = FALSE) {
  lib <- bed_library()
  ok <- dir.exists(file.path(lib, "shiny"))
  if (ok) .libPaths(c(lib, .Library))
  else if (required) stop("The release library is not restored. Run: Rscript tools/restore_library.R")
  invisible(ok)
}
environment_report <- function() {
  pk <- c("shiny", "bslib", "ggplot2", "DT", "markdown", "jsonlite", "zip", "digest", "testthat", "callr")
  v <- vapply(pk, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "missing"), "")
  cat(R.version.string, "\nlibraries:", .libPaths(), "\n", paste(pk, v, collapse = ", "), "\n")
}
