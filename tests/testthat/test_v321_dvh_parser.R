# Version 3.2.1, review F01 (2026-10-26 v3.2 review): a DVH row with a missing or extra field,
# or two structures with the same name, must be refused, never recycled or overwritten.
# Every file here is synthetic.

app_abs <- normalizePath(testthat::test_path("..", ".."))
source(file.path(app_abs, "R", "oar.R"), local = TRUE)

hdr3 <- "Dose [Gy] Relative dose [%] Ratio of Total Structure Volume [%]"
hdr2 <- "Dose [Gy] Ratio of Total Structure Volume [%]"
blk <- function(name, hdr, rows) c(paste("Structure:", name), "", hdr, rows)
ecl <- function(...) parse_dvh_eclipse(c(...))
csv <- function(...) parse_dvh_csv(c(...))

test_that("F01-A: a three-column row missing its relative-dose field is refused (review false parotid pass)", {
  good <- ecl(blk("Parotid", hdr3, c("0 0 100", "10 20 90", "50 100 0")))
  expect_equal(good$Parotid$volume, c(100, 90, 0))
  expect_equal(good$Parotid$dose, c(0, 10, 50))
  expect_error(ecl(blk("Parotid", hdr3, c("0 0 100", "10 90", "50 100 0"))), "line 5.*3 fields")
})

test_that("F01-A: removing any one field, or adding one, is refused in both layouts", {
  rows3 <- c("0 0 100", "10 20 90", "50 100 0")
  for (j in 1:3) {
    r <- rows3; tok <- strsplit(r[2], " ")[[1]]; r[2] <- paste(tok[-j], collapse = " ")
    expect_error(ecl(blk("P", hdr3, r)), "fields", info = j)
  }
  expect_error(ecl(blk("P", hdr3, c("0 0 100", "10 20 90 7", "50 100 0"))), "fields")
  expect_error(ecl(blk("P", hdr2, c("0 100", "10 20 90", "50 0"))), "fields")
  ok2 <- ecl(blk("P", hdr2, c("0 100", "10 90", "50 0")))
  expect_equal(ok2$P$volume, c(100, 90, 0))
})

test_that("F01-A: the volume is read from the column the header names, not the last one", {
  h <- "Dose [Gy] Ratio of Total Structure Volume [%] Relative dose [%]"
  x <- ecl(blk("P", h, c("0 100 0", "10 90 20", "50 0 100")))
  expect_equal(x$P$volume, c(100, 90, 0))
})

test_that("F01-B: two structures with the same name are refused (review false cord pass)", {
  a <- blk("Cord", hdr3, c("0 0 100", "50 100 1", "60 120 0"))
  b <- blk("Cord", hdr3, c("0 0 100", "10 20 1", "20 40 0"))
  expect_error(ecl(a, b), "Cord.*more than once")
  expect_error(ecl(a, blk("  Cord ", hdr3, c("0 0 100", "10 20 1", "20 40 0"))), "more than once")
  expect_error(ecl(a, a), "more than once")
  two <- ecl(a, blk("Brainstem", hdr3, c("0 0 100", "10 20 1", "20 40 0")))
  expect_equal(names(two), c("Cord", "Brainstem"))
  expect_equal(two$Cord$dose, c(0, 50, 60))
})

test_that("F01-B: a blank structure name is refused", {
  expect_error(ecl(blk("", hdr3, c("0 0 100", "10 20 1", "20 40 0"))), "no name")
})

test_that("F01: CSV rows with a missing or extra field, blank or split structures are refused", {
  expect_equal(csv("Dose (Gy),Volume (%)", "0,100", "10,90", "50,0")[[1]]$volume, c(100, 90, 0))
  expect_error(csv("Dose (Gy),Volume (%)", "0,100", "10,90,5", "50,0"), "fields")
  expect_error(csv("Dose (Gy),Volume (%)", "0,100", "10", "50,0"), "fields|non-numeric")
  s <- "Structure,Dose (Gy),Volume (%)"
  expect_error(csv(s, "Cord,0,100", "Cord,50,1", "Cord,60,0", "Lips,0,100", "Lips,5,0",
                   "Cord,0,100", "Cord,10,1", "Cord,20,0"), "more than once")
  expect_error(csv(s, ",0,100", ",50,1", ",60,0"), "no name")
  ok <- csv(s, "Cord,0,100", "Cord,50,1", "Cord,60,0", "Lips,0,100", "Lips,5,0")
  expect_equal(names(ok), c("Cord", "Lips"))
})
