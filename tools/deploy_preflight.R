# Checks run before any deployment. Each failure stops the deploy with a reason.
#   app_name     the app to deploy; must equal expected_app
#   protected    app names this branch must never deploy to (live versions 1 and 2)
#   accounts     data frame from rsconnect::accounts() (columns name, server)
#   account      the shinyapps.io account to use (pinned, never "the first one")
#   dirty        `git status --porcelain --untracked-files=no` output
#   head, passed the current revision, and the revision the full suite last passed on
deploy_preflight <- function(app_name, expected_app, accounts, account, dirty,
                             protected = c("custom_bed_calculator", "custom_bed_calculator_v2"),
                             head = NULL, passed = NULL) {
  if (app_name %in% protected)
    stop("Refusing to deploy over a live earlier version ('", app_name, "').", call. = FALSE)
  if (!identical(app_name, expected_app))
    stop("This branch deploys only to '", expected_app, "', not '", app_name, "'.", call. = FALSE)
  ok <- accounts[accounts$server == "shinyapps.io" & accounts$name == account, , drop = FALSE]
  if (nrow(ok) != 1)
    stop("shinyapps.io account '", account, "' is not configured on this computer. ",
         "Set it once with rsconnect::setAccountInfo() typed in the R console (never saved in a file).",
         call. = FALSE)
  if (length(dirty) > 0)
    stop("Tracked files have uncommitted changes; commit them first so the deployed revision is known:\n",
         paste(dirty, collapse = "\n"), call. = FALSE)
  if (!is.null(head) && !identical(head, passed))
    stop("The full test suite has not passed on this revision (", substr(head, 1, 7), "). Run Rscript tools/run_tests.R first.",
         call. = FALSE)
  invisible(TRUE)
}

# The bundle rsconnect will upload must hold every file the app reads and none
# of the local-only files. `files`: rsconnect::listDeploymentFiles(".").
# (rsconnect silently skips files named manifest.json, which is how the draws
# manifest was once left out of a bundle.)
# `tracked`: `git ls-files` output. Every bundled file must be tracked, so the
# deployed code is exactly the recorded revision (an untracked file would ship
# without being part of any commit).
check_bundle <- function(files, root = ".", tracked = NULL) {
  need <- c("app.R", "about.md",
            file.path("R", list.files(file.path(root, "R"), pattern = "\\.R$")),
            file.path("data", list.files(file.path(root, "data"))),
            file.path("draws", list.files(file.path(root, "draws"))))
  missing <- setdiff(need, files)
  forbidden <- intersect(files, c(".Rhistory", ".Renviron", "external_proofs.yaml", ".suite_passed"))
  if (length(missing)) stop("Files the app reads are missing from the bundle: ", paste(missing, collapse = ", "), call. = FALSE)
  if (length(forbidden)) stop("Local-only files would be deployed: ", paste(forbidden, collapse = ", "), call. = FALSE)
  untracked <- if (is.null(tracked)) character() else setdiff(files, tracked)
  if (length(untracked)) stop("Files not committed to git would be deployed: ", paste(untracked, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}
