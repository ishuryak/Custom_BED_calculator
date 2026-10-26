# Version 3: treatment calendars (R/schedule.R) and the organ-at-risk module
# (R/oar.R). Analytic answers and published values only; no patient data.

app_abs <- normalizePath(testthat::test_path("..", ".."))
source(file.path(app_abs, "R", "bed_models.R"), local = TRUE)
source(file.path(app_abs, "R", "schedule.R"), local = TRUE)
source(file.path(app_abs, "R", "oar.R"), local = TRUE)
lib <- load_constraints(file.path(app_abs, "data", "oar_constraints.csv"))
dvh <- function(f) file.path(app_abs, "tests", "testthat", "fixtures", "dvh", f)

# --- calendars ----------------------------------------------------------------

test_that("the default calendar reproduces the version 1/2 closed form for 1-80 fractions", {
  expect_equal(schedule_time(1:80), overall_time(1:80))
})
test_that("6 fractions per week gives HYPNO's control arm, 33 fractions in 38 days", {
  expect_equal(schedule_time(33, fx_per_week = 6), 38)
})
test_that("a Thursday start puts a weekend inside short courses", {
  expect_equal(schedule_time(5, 5, start_weekday = 1), 5)
  expect_equal(schedule_time(5, 5, start_weekday = 4), 7)
  expect_equal(treatment_days(5, 5, 4), c(1L, 2L, 5L, 6L, 7L))
})
test_that("skipped days lengthen the course by the treatment days they remove", {
  expect_equal(schedule_time(20, skip_days = c(3, 4)), schedule_time(22))  # two weekdays = two more treatment days
  expect_equal(schedule_time(20, skip_days = 6), schedule_time(20))            # a Saturday changes nothing
  expect_error(treatment_days(10, skip_days = 1), "first treatment day")
})
test_that("typed skipped days are parsed and checked", {
  expect_equal(parse_skip_days("12, 13; 30 13"), c(12L, 13L, 30L))
  expect_equal(parse_skip_days(""), integer())
  expect_error(parse_skip_days("1"), "2 to 400")
  expect_error(parse_skip_days("abc"), "whole numbers")
})

# --- LQ conversion ---------------------------------------------------------------

test_that("EQD2 and its inverse are consistent and fixed at 2 Gy per fraction", {
  expect_equal(eqd2_physical(70, 35, 3), 70)
  for (ab in c(0.87, 2, 3, 10)) for (n in c(5, 15, 20, 35))
    expect_equal(eqd2_physical(physical_from_eqd2(50, n, ab), n, ab), 50, tolerance = 1e-12)
})
test_that("HYHOPE's 15-fraction cord limit of 37 Gy is EQD2 41.3 Gy at α/β 2 and 43.0 at 0.87 (plan section B4)", {
  expect_equal(round(eqd2_physical(37, 15, 2), 1), 41.3)
  expect_equal(round(eqd2_physical(37, 15, 0.87), 1), 43.0)
})
test_that("fewer, larger fractions lower the physical limit of a 2-Gy-fraction constraint", {
  lims <- convert_limit(50, 0, c(35, 25, 20, 15), 2)
  expect_true(all(diff(lims) < 0))
  expect_equal(convert_limit(50, 0, 25, 2), 50)   # 25 x 2 Gy is the constraint itself
})

# HYHOPE validation fixture: the per-protocol limits of the 12/15/20-fraction
# tables are not one exact LQ conversion. Each organ's EQD2 (alpha/beta 2 for
# cord and brainstem, 3 otherwise) agrees across the three tables within 7%;
# the spreads are recorded so a change in the library is noticed.
test_that("HYHOPE's three tables agree with one LQ EQD2 per organ within 7%", {
  h <- lib[lib$library == "HYHOPE" & lib$level == "per_protocol", ]
  h$eqd2 <- eqd2_physical(h$limit_gy, h$n_fractions, ifelse(h$ab_class == "cord", 2, 3))
  spread <- tapply(h$eqd2, h$organ, function(e) diff(range(e)) / mean(e))
  expect_equal(length(spread), 5)
  expect_true(all(spread < 0.07))
  expect_equal(round(as.vector(spread[c("Spinal cord", "Mandible")]), 3), c(0.032, 0.047))
})
test_that("every HYHOPE mandatory value matches the protocol tables (pages 13-14)", {
  pp <- lib[lib$library == "HYHOPE" & lib$level == "per_protocol", ]
  key <- setNames(pp$limit_gy, paste(pp$organ, pp$n_fractions))
  expect_equal(unname(key[c("Spinal cord 12", "Spinal cord 15", "Spinal cord 20")]), c(34, 37, 40))
  expect_equal(unname(key[c("Brainstem 12", "Brainstem 15", "Brainstem 20")]), c(38, 43, 46))
  expect_equal(unname(key[c("Mandible 12", "Mandible 15", "Mandible 20")]), c(46.6, 48.9, 54.6))
})
test_that("every library row carries a source, page and quote", {
  expect_true(all(nzchar(lib$source_short) & nzchar(lib$page) & nzchar(lib$quote)))
  expect_true(all(lib$limit_lo_gy <= lib$limit_gy))
})

# --- DVH ---------------------------------------------------------------------------

test_that("bin-by-bin mean EQD2 of a linear differential DVH matches the analytic integral", {
  x <- read_dvh_upload(dvh("eclipse_synthetic.txt"), file.size(dvh("eclipse_synthetic.txt")))
  dd <- cumulative_to_differential(x$Constrictors$dose, x$Constrictors$volume)
  mm <- dvh_means(dd, 20, 3)
  a <- 10; b <- 50; n <- 20; ab <- 3
  exact <- ((b^3 - a^3) / (3 * n) + ab * (b^2 - a^2) / 2) / (2 + ab) / (b - a)
  expect_equal(mm$mean_dose, 30, tolerance = 1e-4)
  expect_equal(mm$mean_eqd2, exact, tolerance = 1e-4)
  # converting the mean alone understates it (LQ is convex in dose per fraction)
  expect_lt(eqd2_physical(30, n, ab), mm$mean_eqd2)
})
test_that("cGy and Gy, % and cc are read; a uniform dose gives its own mean", {
  x <- read_dvh_upload(dvh("eclipse_synthetic.txt"), 1000)
  expect_equal(dvh_means(cumulative_to_differential(x$Parotid_R$dose, x$Parotid_R$volume), 20, 3)$mean_dose, 20, tolerance = 1e-3)
  y <- read_dvh_upload(dvh("eclipse_gy_cc.txt"), 1000)
  expect_equal(attr(y$Cord, "volume_unit"), "cc")
  mx <- dvh_means(cumulative_to_differential(y$Cord$dose, y$Cord$volume), 20, 2)   # D1 (stage A): the top occupied bin, not a midpoint
  expect_equal(c(mx$max_lo, mx$max_hi), c(29.99, 30), tolerance = 1e-9)
  z <- read_dvh_upload(dvh("generic.csv"), 1000)
  expect_equal(dvh_means(cumulative_to_differential(z$Constrictors$dose, z$Constrictors$volume), 20, 3)$mean_dose, 30, tolerance = 1e-4)
})
test_that("files with patient identifiers are refused without echoing them", {
  err <- tryCatch(read_dvh_upload(dvh("eclipse_with_identifiers.txt"), 1000), error = conditionMessage)
  expect_match(err, "patient identifier fields")
  expect_false(grepl("Synthetic|000000", err))
  expect_length(identifier_lines(c("Patient Name         : ", "Patient ID           : ")), 0)
})
test_that("oversize files and CSVs without units are refused", {
  expect_error(read_dvh_upload(dvh("generic.csv"), 3 * 1024^2), "2 MB")
  expect_error(read_dvh_upload(dvh("generic_no_units.csv"), 100), "unit")
  expect_error(cumulative_to_differential(c(0, 1, 2), c(50, 60, 0)), "must not increase")
})

# --- feasibility --------------------------------------------------------------------

test_that("blank inputs are 'not assessed', never a pass", {
  t <- oar_table(lib, list(), 20, 2.75)
  expect_true(all(t$status %in% c("not assessed", "no verified constraint in this library")))
  expect_false(any(t$binds_first))
})
test_that("status words follow the converted limit, with a 5% band", {
  lim <- convert_limit(50, 0, 20, 2)
  st <- function(D) oar_row(lib, "cord", list(mode = "gy", value = D), 20, 2.75, ab = 2)$status
  expect_equal(st(lim * 0.90), "within")
  expect_equal(st(lim * 0.99), "within 5% of the converted constraint")
  expect_equal(st(lim * 1.01), "exceeds the converted constraint")
})
test_that("organs without a verified QUANTEC constraint say so, and HYHOPE fills them", {
  q <- oar_row(lib, "mandible", list(mode = "gy", value = 50), 20, 2.75, "QUANTEC")
  expect_equal(q$status, "not assessed: no verified constraint in this library")   # S9 (stage A): an entered organ is named as unassessed
  h <- oar_row(lib, "mandible", list(mode = "gy", value = 50), 20, 2.75, "HYHOPE")
  expect_equal(h$limit_physical, 54.6, tolerance = 1e-9)   # 20 fractions: HYHOPE's own value
  h17 <- oar_row(lib, "mandible", list(mode = "gy", value = 50), 17, 2.75, "HYHOPE")
  expect_match(h17$flags, "converted from HYHOPE's 15-fraction table")
})
test_that("a mean dose entered as a number is flagged as converted as if uniform", {
  r <- oar_row(lib, "parotid", list(mode = "gy", value = 20), 20, 2.75)
  expect_match(r$flags, "as if uniform")
})
test_that("the alpha/beta flip check finds where the status changes", {
  # 47 Gy cord max in 20 fractions: EQD2 limit 50 Gy is exceeded at low alpha/beta only
  r <- oar_row(lib, "cord", list(mode = "gy", value = 47), 20, 2.35, ab = 2, ab_range = c(0.87, 4))
  expect_match(r$flip, "^exceeds for α/β 0.87-")
  expect_match(r$flip, "within otherwise")
})
test_that("the binding constraint is the one with the smallest relative margin", {
  t <- oar_table(lib, list(cord = list(mode = "gy", value = 30), parotid = list(mode = "gy", value = 22)), 20, 2.75)
  expect_equal(t$organ[t$binds_first], "Parotid (spared side)")
})

test_that("the 2018 Table 1 file has 35 arms, 16 trials and 7283 patients, as the paper states", {
  x <- read.csv(file.path(app_abs, "data", "shuryak2018_table1.csv"))
  expect_equal(nrow(x), 35)
  expect_equal(length(unique(x$ref_no)), 16)
  expect_equal(sum(x$n_patients), 7283)
})
