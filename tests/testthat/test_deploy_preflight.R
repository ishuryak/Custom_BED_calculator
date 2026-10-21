# Deployment preflight (2026-10-21, PLAN_v3 section 3): pinned account, the
# right app name only, and no uncommitted tracked changes.
app_abs <- normalizePath(testthat::test_path("..", ".."))
source(file.path(app_abs, "tools", "deploy_preflight.R"), local = TRUE)
acc <- data.frame(name = c("other", "ishuryak"), server = c("shinyapps.io", "shinyapps.io"))

test_that("a clean tree, the pinned account and the branch's app name pass", {
  expect_true(deploy_preflight("custom_bed_calculator_v2", "custom_bed_calculator_v2", acc, "ishuryak", character()))
})
test_that("the first configured account is never used in place of the pinned one", {
  expect_error(deploy_preflight("custom_bed_calculator_v2", "custom_bed_calculator_v2", acc[1, ], "ishuryak", character()),
               "not configured")
})
test_that("the version 1 site and other app names are refused", {
  expect_error(deploy_preflight("custom_bed_calculator", "custom_bed_calculator_v2", acc, "ishuryak", character()),
               "deploys only")
  expect_error(deploy_preflight("custom_bed_calculator", "custom_bed_calculator", acc, "ishuryak", character()),
               "version 1")
})
test_that("uncommitted tracked changes stop the deploy", {
  expect_error(deploy_preflight("custom_bed_calculator_v2", "custom_bed_calculator_v2", acc, "ishuryak", " M app.R"),
               "uncommitted")
})
test_that("deploy.R uses the preflight and the pinned account, not the first account", {
  s <- readLines(file.path(app_abs, "tools", "deploy.R"))
  expect_true(any(grepl("deploy_preflight(", s, fixed = TRUE)))
  expect_false(any(grepl("acc$name[1]", s, fixed = TRUE)))
  expect_true(any(grepl("account = ACCOUNT", s, fixed = TRUE)))
})
test_that(".Rhistory and .Renviron are neither tracked nor deployed", {
  expect_true(all(c(".Rhistory", ".Renviron") %in% readLines(file.path(app_abs, ".rscignore"))))
  expect_true(all(c(".Rhistory", ".Renviron") %in% readLines(file.path(app_abs, ".gitignore"))))
})
