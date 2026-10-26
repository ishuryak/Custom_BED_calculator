# Treatment calendars (version 3). The models use only the overall time T
# (inclusive days, first to last treatment day), so a calendar matters through T.
#
# A course is m fractions, one per treatment day, on the allowed weekdays
# (5 per week: Monday-Friday; 6 per week: Monday-Saturday), starting on
# `start_weekday` (1 = Monday ... 6 = Saturday). `skip_days` lists calendar
# days (1 = first treatment day) on which no fraction is given: holidays,
# machine down time or a treatment gap. Twice-daily treatment is not covered.

WEEKDAY_NAMES <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")

# One validation contract for skipped days, shared by the interface parser,
# the direct API and replayed files: finite whole numbers from 2 to 400,
# duplicates removed. An empty value (NULL, integer(0), an empty JSON array)
# means no skipped days.
check_skip_days <- function(x) {
  if (is.null(x) || length(x) == 0) return(integer())
  if (is.list(x)) x <- unlist(x)
  if (!is.numeric(x) || any(!is.finite(x)) || any(x != round(x)))
    stop("Skipped days must be whole numbers from 2 to 400 (day 1 is the first treatment day).")
  if (any(x == 1)) stop("Day 1 is the first treatment day by definition and cannot be skipped.")
  if (any(x < 2) || any(x > 400))
    stop("Skipped days must be whole numbers from 2 to 400 (day 1 is the first treatment day).")
  sort(unique(as.integer(x)))
}

# Calendar day numbers (1 = first treatment day) of the m treatment days.
treatment_days <- function(m, fx_per_week = 5, start_weekday = 1, skip_days = integer()) {
  one_whole <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x) && x == round(x)
  if (!one_whole(m) || m < 1) stop("The number of fractions must be one whole number of at least 1.")
  if (!one_whole(fx_per_week) || !fx_per_week %in% c(5, 6)) stop("Fractions per week must be 5 or 6.")
  if (!one_whole(start_weekday) || !start_weekday %in% seq_len(fx_per_week))
    stop("The first treatment day must be a weekday number from 1 to ", fx_per_week, ".")
  skip_days <- check_skip_days(skip_days)
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

# Candidate calendars of both searches (deterministic and uncertainty): the
# regular calendar of the sidebar (fractions per week and start day). Skipped
# days describe the evaluated regimen only and are never applied to candidates.
candidate_calendar <- function(cal) {
  list(fx_per_week = cal$fx_per_week, start_weekday = cal$start_weekday, skip_days = integer())
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
  if (any(is.na(v)) || any(v == 1))
    stop("Skipped days must be whole numbers from 2 to 400 (day 1 is the first treatment day).")
  check_skip_days(v)
}

# "m × d Gy in T days" for messages and table headers.
regimen_text <- function(r) sprintf("%d × %g Gy in %g days", as.integer(r$m), r$d, r$T)
