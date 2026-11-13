# Checks run before any deployment. Each failure stops the deploy with a reason.
#   app_name     the app to deploy; must equal expected_app
#   protected    app names this branch must never deploy to (live versions 1, 2 and 3)
#   accounts     data frame from rsconnect::accounts() (columns name, server)
#   account      the shinyapps.io account to use (pinned, never "the first one")
#   dirty        `git status --porcelain --untracked-files=no` output
#   head, passed the current revision, and the revision the full suite last passed on
deploy_preflight <- function(app_name, expected_app, accounts, account, dirty,
                             protected = c("custom_bed_calculator", "custom_bed_calculator_v2", "custom_bed_calculator_v3"),
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
  if (!is.null(head) && (length(head) != 1 || is.na(head) || !nzchar(head) || !identical(head, passed)))
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
# without being part of any commit). deployed_revision.txt is the one exception:
# deploy.R writes it just before upload to name the revision (it is gitignored).
# Documents the pages link to and serve (review 6.4, 2026-11-09: a bundle missing both refit
# explainer copies was accepted). Each must be in the bundle and not empty.
SERVED_ASSETS <- c("www/refit_explainer.html", "www/refit_explainer.docx")
# Text the pages are built from (version 4.6.1): an empty file would render an empty page section,
# so each must be in the bundle and not empty, like the served documents.
PAGE_CONTENT <- c("about.md", "refit_in_brief.md", "early_fit_in_brief.md", "how_fitted.md", "draws/model_methods.json")
# Files (not folders) directly in a folder: a bundle lists files, so a subfolder such as
# data/early_tox (research data the app does not read) must not be required by name (4.6.0).
top_files <- function(root, dir) {
  x <- list.files(file.path(root, dir))
  file.path(dir, x[!dir.exists(file.path(root, dir, x))])
}
check_bundle <- function(files, root = ".", tracked = NULL) {
  need <- c("app.R", PAGE_CONTENT, SERVED_ASSETS,
            file.path("R", list.files(file.path(root, "R"), pattern = "\\.R$")),
            top_files(root, "data"), top_files(root, "draws"))
  missing <- setdiff(need, files)
  content <- c(SERVED_ASSETS, PAGE_CONTENT)
  empty <- content[content %in% files & !(file.exists(file.path(root, content)) & file.size(file.path(root, content)) > 0)]
  if (length(empty)) stop("Linked documents are empty or unreadable: ", paste(empty, collapse = ", "), call. = FALSE)
  forbidden <- intersect(files, c(".Rhistory", ".Renviron", "external_proofs.yaml", ".suite_passed"))
  if (length(missing)) stop("Files the app reads are missing from the bundle: ", paste(missing, collapse = ", "), call. = FALSE)
  if (length(forbidden)) stop("Local-only files would be deployed: ", paste(forbidden, collapse = ", "), call. = FALSE)
  untracked <- if (is.null(tracked)) character() else setdiff(files, c(tracked, "deployed_revision.txt"))
  if (length(untracked)) stop("Files not committed to git would be deployed: ", paste(untracked, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

# The methods record (draws/model_methods.json) must be current: every source file it names and
# every draw file it describes must still have the recorded bytes (review B02, 2026-11-12). The full
# rebuild check is `python3 tools/build_model_methods.py --check`; this is its dependency-free
# counterpart, run by tools/deploy.R and by the test suite.
check_methods_record <- function(root = ".") {
  sha <- function(f) digest::digest(readBin(f, "raw", n = file.size(f)), algo = "sha256", serialize = FALSE)
  p <- file.path(root, "draws", "model_methods.json")
  m <- tryCatch(jsonlite::fromJSON(p, simplifyVector = FALSE), error = function(e) NULL)
  if (!is.list(m) || !is.list(m$sources) || !length(m$sources))
    stop("The methods record is missing or has no source list.", call. = FALSE)
  want <- c(m$sources, list("draws/tumor_di_primary.csv" = m$tumor$di$draws_sha256,
                            "draws/tumor_dd_primary.csv" = m$tumor$dd$draws_sha256,
                            "draws/early_fit_draws.csv" = m$early$draws$sha256))
  stale <- names(want)[vapply(names(want), function(f) {
    path <- file.path(root, f)
    !(is.character(want[[f]]) && file.exists(path) && identical(sha(path), want[[f]]))
  }, TRUE)]
  if (length(stale)) stop("The methods record is out of date for: ", paste(stale, collapse = ", "),
                          ". Run python3 tools/build_model_methods.py.", call. = FALSE)
  invisible(TRUE)
}

# The revision a passing suite may certify (review R09, 2026-10-26). Fails closed: a git that
# cannot run, a non-zero exit, output that is not one 40-character hash, or uncommitted tracked
# changes all give ok = FALSE with the reason, and no suite marker is written.
# git: function(args) -> list(status, out); replaceable in tests.
git_run <- function(args) {
  out <- tryCatch(suppressWarnings(system2("git", args, stdout = TRUE, stderr = TRUE)),
                  error = function(e) structure(conditionMessage(e), status = 127L))
  list(status = if (is.null(attr(out, "status"))) 0L else as.integer(attr(out, "status")), out = as.character(out))
}
release_revision <- function(git = git_run) {
  h <- git(c("rev-parse", "HEAD"))
  if (h$status != 0) return(list(ok = FALSE, why = "git rev-parse failed"))
  rev <- trimws(h$out)
  if (length(rev) != 1 || !grepl("^[0-9a-f]{40}$", rev)) return(list(ok = FALSE, why = "git did not return one revision hash"))
  d <- git(c("status", "--porcelain", "--untracked-files=no"))
  if (d$status != 0) return(list(ok = FALSE, why = "git status failed"))
  if (length(d$out[nzchar(d$out)]) > 0) return(list(ok = FALSE, why = "tracked files have uncommitted changes"))
  list(ok = TRUE, rev = rev)
}

# Everything deploy.R needs from git, checked before anything is written or uploaded (review M01
# of the 2026-10-26 v3.2 review: deploy.R used unchecked system() calls). Returns list(rev, tracked)
# or stops: HEAD must be one 40-character hash with a clean tracked tree (release_revision()), the
# suite marker must name exactly that hash, and `git ls-files` must succeed and list files.
deploy_identity <- function(git = git_run, marker = ".suite_passed") {
  rr <- release_revision(git)
  if (!rr$ok) stop("Cannot identify the revision to deploy: ", rr$why, ".", call. = FALSE)
  passed <- if (file.exists(marker)) readLines(marker, warn = FALSE) else character()
  if (length(passed) != 1 || !grepl("^[0-9a-f]{40}$", passed) || !identical(passed, rr$rev))
    stop("The full test suite has not passed on this revision (", substr(rr$rev, 1, 7),
         "). Run Rscript tools/run_tests.R first.", call. = FALSE)
  ls <- git(c("ls-files"))
  if (ls$status != 0 || !length(ls$out[nzchar(ls$out)])) stop("git ls-files failed or listed nothing.", call. = FALSE)
  list(rev = rr$rev, passed = passed, tracked = ls$out[nzchar(ls$out)])
}

# The tools library for rsconnect (review 6.4): create it BEFORE putting it on the library path.
# R drops a directory that does not exist from .libPaths(), so a first-time install there was not
# found in the same process. Then install if needed and confirm rsconnect loads from it.
# install: function(pkg, lib) -> installs; replaceable in tests.
ensure_deploy_tools <- function(lib, install = function(pkg, lib)
                                  utils::install.packages(pkg, lib = lib, repos = LOCK_REPO,
                                    type = if (.Platform$OS.type == "windows") "win.binary" else getOption("pkgType"))) {
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(lib)) stop("Cannot create the deploy tools library ", lib, call. = FALSE)
  .libPaths(c(.libPaths()[1], lib, .Library))
  if (!normalizePath(lib, winslash = "/") %in% normalizePath(.libPaths(), winslash = "/"))
    stop("The deploy tools library is not on the library path: ", lib, call. = FALSE)
  if (!requireNamespace("rsconnect", quietly = TRUE)) install("rsconnect", lib)
  if (!requireNamespace("rsconnect", quietly = TRUE))
    stop("rsconnect could not be installed or loaded from ", lib, call. = FALSE)
  invisible(lib)
}
