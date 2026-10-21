# Treatment calendars (version 3). The models use only the overall time T
# (inclusive days, first to last treatment day), so a calendar matters through T.
#
# A course is m fractions, one per treatment day, on the allowed weekdays
# (5 per week: Monday-Friday; 6 per week: Monday-Saturday), starting on
# `start_weekday` (1 = Monday ... 6 = Saturday). `skip_days` lists calendar
# days (1 = first treatment day) on which no fraction is given: holidays,
# machine down time or a treatment gap. Twice-daily treatment is not covered.

WEEKDAY_NAMES <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")

# Calendar day numbers (1 = first treatment day) of the m treatment days.
treatment_days <- function(m, fx_per_week = 5, start_weekday = 1, skip_days = integer()) {
  stopifnot(length(m) == 1, m >= 1, m == round(m), fx_per_week %in% c(5, 6),
            start_weekday %in% 1:fx_per_week)
  skip_days <- unique(as.integer(skip_days))
  if (1L %in% skip_days) stop("Day 1 is the first treatment day by definition and cannot be skipped.")
  days <- integer(m); k <- 0L; day <- 0L
  while (k < m) {
    day <- day + 1L
    wd <- (start_weekday - 1L + day - 1L) %% 7L + 1L
    if (wd <= fx_per_week && !(day %in% skip_days)) { k <- k + 1L; days[k] <- day }
  }
  days
}

# Overall time (inclusive days) for a calendar. Vectorized over m.
schedule_time <- function(m, fx_per_week = 5, start_weekday = 1, skip_days = integer()) {
  vapply(m, function(mm) {
    d <- treatment_days(mm, fx_per_week, start_weekday, skip_days)
    as.numeric(d[length(d)])
  }, numeric(1))
}

# The default calendar of versions 1 and 2 (5 per week from a Monday) has a
# closed form; overall_time() in bed_models.R must agree with schedule_time().
default_calendar <- function() list(fx_per_week = 5, start_weekday = 1, skip_days = integer())

# Readable description of a calendar.
calendar_text <- function(cal) {
  s <- sprintf("%d fractions/week, starting %s", cal$fx_per_week, WEEKDAY_NAMES[cal$start_weekday])
  if (length(cal$skip_days)) s <- paste0(s, sprintf(", no treatment on day%s %s",
                                                     if (length(cal$skip_days) > 1) "s" else "",
                                                     paste(sort(cal$skip_days), collapse = ", ")))
  s
}

# Parse a typed list of skipped days ("12, 13 30") into integers; errors name the problem.
parse_skip_days <- function(txt) {
  txt <- trimws(if (is.null(txt)) "" else txt)
  if (!nzchar(txt)) return(integer())
  parts <- strsplit(txt, "[,; ]+")[[1]]
  parts <- parts[nzchar(parts)]
  v <- suppressWarnings(as.numeric(parts))
  if (any(is.na(v)) || any(v != round(v)) || any(v < 2) || any(v > 400))
    stop("Skipped days must be whole numbers from 2 to 400 (day 1 is the first treatment day).")
  sort(unique(as.integer(v)))
}

# "m × d Gy in T days" for messages and table headers.
regimen_text <- function(r) sprintf("%d × %g Gy in %g days", as.integer(r$m), r$d, r$T)
