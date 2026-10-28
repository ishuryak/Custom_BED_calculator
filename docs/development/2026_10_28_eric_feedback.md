# Eric's calculator feedback, 2026-10-28

Base: `4568819dbcc30dc2fbd9ed5db683acd4d3bd50e8` on
`fixes_v32_review_2026_10_26` (pending version 3.2.1).

- Report-card tumor threshold status now uses `tumor_endpoint()`, matching the displayed
  difference. Nonpositive or unavailable BEDs read `not evaluable`, even when their raw
  subtraction would meet the threshold. Positive evaluable values retain the same threshold
  comparison and numerical tolerance.
- The uncertainty joint-share display defaults to both paired models. This changes the
  initially displayed joint share, not the equations, draws, or search selection rules.
- The uncertainty ranking has an explicit heading and a note immediately above the results
  explaining that listed doses meet the toxicity share but may miss the tumor threshold.
- The report card explains that its refit interval is separate from the sidebar preset.
- A full deployed commit ID links to its exact GitHub source from Model and sources.
  A missing or malformed revision file does not produce a source link. The moving v3 branch
  is labeled development source. README points to the current v3 website.
- The base already includes the DVH row-count and duplicate-name fixes; they are retained.

## Verification

New tests: `tests/testthat/test_eric_feedback.R` (28 expectations including 3 package warnings).
They exercise the negative-BED report-card case, valid pass/fail/boundary statuses, unavailable
DD, the rendered refit note, both-model UI default, search note placement, and exact source links.

Related suite:

```r
testthat::test_dir("tests/testthat",
  filter = "eric_feedback|v321_dvh_parser|v321_golden|v3_app|v3_uncertainty|v321_ui_items",
  reporter = "summary", stop_on_failure = FALSE)
```

Result: 2204 expectations, 0 failures, 0 errors, 0 skipped, 4 warnings.
Three warnings concern installed packages built under R 4.4.3 (host R 4.4.1);
one concerns an unwritable existing Sass cache during UI rendering. The rendered UI assertions
still pass. The locked release library is absent, so this is an installed-library check,
not release certification. Default numerical BED checks pass unchanged.

At initial preparation, the connected account had read-only access and a patch handoff was
prepared. Eric subsequently confirmed the invitation was accepted; write access was verified
on 2026-10-28. These changes are submitted on a separate review branch. The landing-page
README update targets `main` separately, since it has a different history from this v3 branch.
No website deployment is part of this change; the locked release checks remain necessary.
