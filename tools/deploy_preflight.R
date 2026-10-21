# Checks run before any deployment. Each failure stops the deploy with a reason.
#   app_name   must be the app this branch deploys (never the version 1 site)
#   accounts   data frame from rsconnect::accounts() (columns name, server)
#   account    the shinyapps.io account to use (pinned, never "the first one")
#   dirty      output of `git status --porcelain` for tracked files
deploy_preflight <- function(app_name, expected_app, accounts, account, dirty) {
  if (!identical(app_name, expected_app))
    stop("This branch deploys only to '", expected_app, "', not '", app_name, "'.", call. = FALSE)
  if (app_name %in% c("custom_bed_calculator"))
    stop("Refusing to deploy over the version 1 site.", call. = FALSE)
  ok <- accounts[accounts$server == "shinyapps.io" & accounts$name == account, , drop = FALSE]
  if (nrow(ok) != 1)
    stop("shinyapps.io account '", account, "' is not configured on this computer. ",
         "Set it once with rsconnect::setAccountInfo() typed in the R console (never saved in a file).",
         call. = FALSE)
  if (length(dirty) > 0)
    stop("Tracked files have uncommitted changes; commit them first so the deployed revision is known:\n",
         paste(dirty, collapse = "\n"), call. = FALSE)
  invisible(TRUE)
}
