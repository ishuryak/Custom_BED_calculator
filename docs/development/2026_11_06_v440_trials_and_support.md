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
- Certifying run on the committed revision `ab9152e` (2026-11-06, Windows, R 4.5.1, pinned library
  `library-2026-10-15-R4.5.1`): `Rscript tools/run_tests.R` gave 6188 expectations, 0 failed, 0 errors,
  0 skipped, 0 warnings, and wrote the suite marker for `ab9152e`. A deploy from a later revision
  (this record, or a merge into `v4`) needs the suite rerun on that revision, as `tools/deploy.R` requires.

## Not done

- Not deployed. Deployment needs the shinyapps.io token, which the 2026-11-06 audit found exposed in
  `.Rhistory` on `main`. Rotate it first.
- The user guide and About page do not yet describe the dialog.

## Second round (Igor Shuryak, 2026-11-06)

Requests: drop the Chemotherapy and Other drugs columns and use footnotes for the few trials that
need one; footnote the trials that are in the refit but were not in the 2018 and 2019 papers with
Hall and Brenner; make the three 5 Gy limits of the main page changeable, each on its own, if easy.

- **Footnotes.** The dialog's table has seven columns. Mark a: DAHANCA 6&7 (nimorazole for all
  patients except glottic cancer, in both arms). Mark b: RTOG 83-13 (no chemotherapy within 6 weeks
  of randomization, earlier chemotherapy allowed). The notes are generated from
  `data/fit_trial_list.csv`, so a future trial with a drug gets a note without code changes.
- **Not in the 2018 and 2019 papers.** `shuryak2018_ref` in `data/fit_trial_list.csv` holds each
  trial's reference number in Shuryak, Hall and Brenner 2018 (Radiother Oncol 127:20-26, PMID
  29534828), read from its full text (Table 1 and references [23]-[55]). The 2019 paper (Int J
  Radiat Oncol Biol Phys 104:272-278, PMID 30776451; full text supplied by Igor Shuryak,
  `Publications/Articles/hypofrac.pdf`) used the same data set: its Methods name 16 randomized
  trials from 16 publications, 7283 patients and 35 regimens, and its references 5, 6 and 12-25 are
  the 2018 Table 1 publications. Eight refit trials are in neither paper and are starred: EORTC 22811,
  Ghoshal, IAEA-ACC, KROG 0201, Princess Margaret, RTOG 83-13, RTOG 88-09 and RTOG 9512. This
  agrees with Section 2.3 of the refit explainer. CHART was cited in 2018 through Dische 1997 and is
  counted as the same trial as Saunders 2010.
- **Changeable limits.** Three boxes under the conditions set the tumor loss and the early and late
  increases (Gy of BED), each 0 to 20 in 0.5 steps, default 5, with "Reset all to 5". The value
  reaches `v4_regimens()`, the refit support (tumor margin), the intro and criteria text, the row
  card, the plot's dashed lines, the Report card transfer (`v4_report_inputs(m, d, margins)`,
  `v4_transfer_scope(margins)`) and the ZIP export. An out-of-range value stops the list with a
  message.
- **Export compatibility (judgement call).** The margins were already recorded and hashed in every
  `main-1` scenario, so the schema and analysis version stay `main-1` / `1.0`. Replay now checks the
  recorded margins against the 0 to 20 range and requires the support margin to equal the tumor
  margin, instead of requiring 5/5/5. Exports made with 4.3.x (5/5/5) replay unchanged. A 4.3.0 test
  that expected an edited tumor margin to be refused as a "fixed rule" now expects the refusal for
  the mismatched support margin.
- **Range 0 to 20 Gy (judgement call).** Negative margins (demanding a gain) and values above 20 are
  refused. Either bound can be widened in `V4_MARGIN_RANGE` if wanted.
- Docs: About, the user guide (Markdown and its Word copy, rebuilt with pandoc) and the README say
  the limits default to 5 and can be changed.

Red without the change: with `app.R`, `R/v4_home.R` and `R/v4_export.R` from `ce41191`, the new and
changed tests in `test_v44_trials_support.R` gave 9 failures and 4 errors (drug columns present, no
footnotes, no `v4_check_margins`, the page ignoring `home_m_*`). Full suite on the working tree:
6222 expectations, 1 failed (the 4.3.0 expectation above, updated), 0 errors. The certifying run on
the committed revision is recorded below.

## Third round: the plot went blank while a limit was retyped (reported by Igor Shuryak, 2026-11-06)

**Symptom.** Changing a limit could leave the plot section empty. Reproduced in headless Chrome
(`chromote`) against the local app at `8da8c18`: whenever a limit box was empty or out of range, for
example after Backspace while retyping, or after typing 25, `margins()` failed validation. The plot
section (`#home_plot_ui`) shrank to an empty error element, and the table kept showing the previous
list with no explanation. On the server all 150 combinations of valid limits tried (0 to 20) drew
the plot, so the cause was the transient invalid value, not a limit combination.

**Fix.** An observer checks the three boxes. A valid set becomes the limits in use. An empty or
out-of-range box keeps the last valid limits, so the list and the plot stay, and a red note beside
the boxes names the problem and the limits still used. The check `v4_check_margins()` is unchanged.

**Evidence.** `test_v44_trials_support.R::the main page follows the limits the user sets, and passes
them on` now sets the tumor box to NA, −1 and 25 and expects the previous limits, the same rows, the
plot output and the message. With `app.R` from `8da8c18` it errored (`home()` failed on the invalid
value). With the fix it passes. The browser probe after the fix: plot image present (1100 × 460) in
all 17 states, including the three invalid ones, with the note shown only in those.
