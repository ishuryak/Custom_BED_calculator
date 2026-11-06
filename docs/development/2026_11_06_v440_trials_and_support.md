# Version 4.4.0: refit support always shown, and a list of the refitted trials (2026-11-06)

Branch `v4-trials-support`, from `v4` at `ba9fe51` (4.3.1). Not deployed: the live site stays on
its current release until this branch is approved.

## Requests (Igor Shuryak, 2026-11-06)

1. Show the Bayesian refit support column on the main page, so it is visible without Details.
2. Under Details, add a page that lists the trials the current DI and DD tumor models were fitted to:
   each trial's reference, the compared arms (dose, fractionation, time), and the chemotherapy used,
   opened by a click.
3. Build it on a new branch, and make it the live site only if it looks good.

## What was done

- **Refit support column.** `Refit support` is now added to the main table before the Details
  check. The model column ("Tumor condition holds in") and Notes stay under Details. The
  "Refit support." footnote moved out of the Details block, since the column it explains is now
  always shown. It now says that the low-support note appears in Notes (with Details) and that the
  trials are listed under Details.
- **Trial list.** With Details on, a button "Trials the tumor models were refitted to" opens a
  large dialog (`modalDialog`, size xl). One row per fitted arm: trial and reference (spanning the
  trial's arms, with PubMed and DOI links), arm, total dose, dose per fraction, fractions, overall
  time, time basis, chemotherapy and other drugs.
- **Where the numbers come from.** Doses, fraction numbers and times are read at run time from the
  fitting file named in `draws/draws_manifest.json` (`refit_trials_2026_10_25c_main.csv`), the same
  file `tumor_data_support()` reads, so the list cannot differ from what was fitted and moves with
  the data. Only what that file lacks is stored separately, in `data/fit_trial_list.csv`, built by
  `tools/build_fit_trial_list.py`:
  - citation and DOI from the PubMed records saved by the trial search (`pubmed_records.json`,
    `supp_records.json`), not typed.
  - chemotherapy from the `chemo` field of the full-text extraction rows of each arm
    (`all_extract.json`), matched on PubMed ID and exact arm name. The builder refuses an unmatched
    arm, matched rows that disagree, and a trial without a PubMed record.
  - other drugs: nimorazole in DAHANCA 6&7, from the population note that
    `fit/build_provenance.py` records (read in the publication on 2026-10-26).
  The app refuses a list whose arms do not match the fitting file one to one.

## Decisions and judgement calls (made by Claude under the request above, open to change)

- **What "fitted to" means.** The default tumor values (scenario D18) were chosen by the developer
  and were not fitted to trial data. The 19 trials are those of the Bayesian refit of DI and DD,
  which gives the Refit support column and the refit tumor preset. The dialog says this, rather
  than implying the default numbers came from these trials. The button label says "refitted".
- **A dialog, not a menu page.** The request placed the list under Details, which is a check box on
  the main page, not a menu. A button there opening a dialog keeps the list one click from the
  table it explains. Alternative considered: a new side-menu page beside "Model and sources",
  which would be reachable without Details but further from the table.
- **Values as fitted.** Non-whole fraction numbers (DAHANCA 6&7, EORTC 22811, KROG 0201 and the
  CAIR vs concomitant boost trial) are shown with ≈ and explained as averages, rather than replaced
  by the prescriptions. Mixed fraction sizes are shown as their components (for example
  25 × 1.8 + 15 × 1.5 for Ghoshal's concomitant boost arm).
- **Chemotherapy.** The extraction records no chemotherapy in any of the 45 arms (the trials were
  selected as radiotherapy-only). RTOG 83-13 allowed chemotherapy only more than 6 weeks before
  randomization, which is stated. Nimorazole (a radiosensitizer, not chemotherapy) is given in a
  separate column, so the chemotherapy column stays literal.
- **Extraction folder.** The builder reads the extraction and PubMed records outside the repository
  (`../References_refit_trials/search`). The committed `data/fit_trial_list.csv` is the output, and
  `tests/python/test_fit_trial_list_builder.py::test_committed_list_is_what_the_builder_makes`
  checks it is reproduced exactly where that folder is present (skipped elsewhere, including CI).
  The repo's own `fit/build_refit_dataset.py` could not be reused for the matching because the
  local extraction copy lacks `supp_extract_2026_10_24.json`. All 45 arms are in `all_extract.json`.

## Checks

- Two defects were caught by the new tests before commit: `formatC` padded the dose and time cells
  with spaces, and the builder first wrote "earlier chemotherapy allowed earlier".
- Red without the change: with `app.R` and `R/v4_home.R` reverted to 4.3.1 and the new tests kept,
  `test_v44_trials_support.R` gave 5 failures (support column and footnote hidden without Details,
  no button) and 3 errors (`v4_fit_trials` absent). Builder mutations (dropping the nimorazole entry,
  disabling the disagreement check) each turned one builder test red.
- Suite counts and revision: `docs/FIXES.md`, 2026-11-06 entry.

## Not done

- Not deployed. Deployment needs the shinyapps.io token, which the 2026-11-06 audit found exposed in
  `.Rhistory` on `main`. Rotate it first.
- The user guide and About page do not yet describe the dialog.
