# Version 3.2.1, review M01 (2026-10-26 v3.2 review): deploy identity from checked git calls only.
# Mocked git; nothing is deployed.

app_abs <- normalizePath(testthat::test_path("..", ".."))
source(file.path(app_abs, "tools", "deploy_preflight.R"), local = TRUE)
H <- paste(rep("a1", 20), collapse = "")
fake_git <- function(head = list(status = 0L, out = H), status = list(status = 0L, out = character()),
                     ls = list(status = 0L, out = c("app.R", "about.md")))
  function(args) switch(args[1], "rev-parse" = head, status = status, "ls-files" = ls)
marker <- function(x) { f <- tempfile(); if (!is.null(x)) writeLines(x, f); f }

test_that("M01: a clean, certified revision gives its hash and tracked files", {
  id <- deploy_identity(fake_git(), marker(H))
  expect_identical(id$rev, H); expect_identical(id$tracked, c("app.R", "about.md"))
})

test_that("M01: every failed or empty git answer, and every bad marker, stops before deployment", {
  expect_error(deploy_identity(fake_git(head = list(status = 128L, out = "fatal")), marker(H)), "revision")
  expect_error(deploy_identity(fake_git(head = list(status = 0L, out = character())), marker(H)), "revision")
  expect_error(deploy_identity(fake_git(status = list(status = 1L, out = character())), marker(H)), "revision")
  expect_error(deploy_identity(fake_git(status = list(status = 0L, out = " M app.R")), marker(H)), "uncommitted")
  expect_error(deploy_identity(fake_git(), marker(NULL)), "has not passed")
  expect_error(deploy_identity(fake_git(), marker(character())), "has not passed")
  expect_error(deploy_identity(fake_git(), marker("abc1234")), "has not passed")
  expect_error(deploy_identity(fake_git(), marker(sub("a1$", "b2", H))), "has not passed")
  expect_error(deploy_identity(fake_git(ls = list(status = 1L, out = character())), marker(H)), "ls-files")
  expect_error(deploy_identity(fake_git(ls = list(status = 0L, out = character())), marker(H)), "ls-files")
})

test_that("M01: preflight refuses matching empty identities, and deploy.R uses only checked git calls", {
  acc <- data.frame(name = "ishuryak", server = "shinyapps.io")
  expect_error(deploy_preflight("custom_bed_calculator_v3", "custom_bed_calculator_v3", acc, "ishuryak", character(),
                                head = character(), passed = character()), "has not passed")
  s <- readLines(file.path(app_abs, "tools", "deploy.R"))
  expect_false(any(grepl("system(\"git", s, fixed = TRUE)))
  expect_true(any(grepl("deploy_identity()", s, fixed = TRUE)))
  # the identity is established before deployed_revision.txt is written
  expect_lt(grep("deploy_identity()", s, fixed = TRUE)[1], grep("deployed_revision.txt", s, fixed = TRUE)[1])
})
