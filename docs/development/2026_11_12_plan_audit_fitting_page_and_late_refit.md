# Audit of the 2026-11-12 plan: a findable Bayesian-method page (Issue 1) and a late-toxicity refit (Issue 2)

> **Superseded in part, 2026-11-12 (later).** The independent review
> `2026_11_12_bayesian_methods_and_late_toxicity_review.md` (section 8) corrects F1, F2, F4, F5,
> L2, L3, L4, L6 and decision D4 of this audit. Where the two differ, the review holds; the
> corrections were checked against the code and are recorded in `2026_11_12_v461_methods_page.md`.
> In particular, the research-mode rows of L2 were not invalid (`early_cal` is calendar information):
> the research list has 9 doses at every late α/β tested.

Written 2026-11-12 against revision `52f0d3b` (4.6.0, live at `c80f63b`; the two differ only in
walkthrough screenshots and a log). Status: AUDIT ONLY. No app file was changed. Nothing is committed.

## 0. The plan as proposed (version 1, in chat on 2026-11-12)

- **P1.** A new menu item "How the models were fitted (Bayesian)". It opens one page covering what a
  Bayesian fit is and why it was used, the tumor fit (described as a "hierarchical binomial model"),
  its likelihood, priors, Stan settings and quality checks, what is driven by a prior, how the 4,000
  draws feed Refit support and Uncertainty, the same for the early fit, and links to the explainer.
- **P2.** A "How the models were fitted (Bayesian)" button under Details on the main page.
- **P3.** A paragraph in `refit_in_brief.md` naming the method.
- **P4.** Tests, walkthrough, FIXES row, development record, version 4.6.1. Deploy on separate approval.
- **L1.** Issue 2: the data cannot be reused. Instead, a feasibility step (screen the PDFs,
  double-extract late counts, check identification), then a staged plan like the early fit.

Premise of P1–P2: "Neither the menu nor the Details toggle points to any of this."

## 1. Scope

**Checked:**
- `app.R` (menu, main page, dialogs), `about.md`, `refit_in_brief.md`, `early_fit_in_brief.md`,
  `.rscignore`, `tools/deploy.R`, `tests/testthat/test_deploy_preflight.R`, the walkthrough, and
  `fit/stan/di.stan` and `dd.stan`.
- Every PDF in `References_refit_trials/` and `References/` (120 files), screened with pypdf for
  late-effect terms (`scratchpad/late_scan.py`).
- The late and early dose ceilings, computed with the app's own `compute_all()` in R 4.5.1
  (`scratchpad/binding_late.R`).

**Not checked:**
- Whether any trial PDF gives extractable per-arm late counts. The screen counts terms only.
- The live site in a browser. The code is the deployed code, but its rendering was not reviewed.
- Ceilings under the research (trial-count) early model (see L2).

## 2. Findings on Issue 1

### F1 (major, verified). The premise is partly wrong: a visible link already exists

**Finding.**
- The main page's settings line already shows "How was the tumor model fitted?" and "Trials the tumor
  model was fitted to", without opening Details (`app.R:1652–1653`).
- The same link appears on Uncertainty (`app.R:397`) and on Settings and sources (`app.R:401`).
- What the user cannot find is the METHOD. The tumor dialog (`refit_in_brief.md`) never uses the word
  "Bayesian" and mentions no prior, likelihood, sampler or check.
- The Bayesian rationale (explainer section 5) sits in a separate HTML, linked by a grey line at the
  bottom of the main page and by one sentence inside Model and sources (`about.md:17`).
- The menu has no entry for it.

**Fix.**
- Restate the problem as: the link exists, but what it opens does not explain the method, and the
  menu has no entry.
- Keep P1 (menu page) and P2 (Details button, as Igor asked).
- Retarget the three existing "How was the tumor model fitted?" links to the same new page, under one
  label: "How the models were fitted (Bayesian)".
- Retire the short tumor dialog as a separate destination: its text becomes a section of the new page.
- Do not add a second link with a different name beside the existing one.

**Test.** `tests/testthat/test_v461_fitting_page.R::test_every_fitted_link_opens_the_method_page`:
each of `how_fitted_home`, `how_fitted_unc`, `how_fitted_settings`, the new Details button and the
menu entry select `section == "fitted"` (`shiny::testServer`). To see it fail, leave one link on the
old modal.

**Verification.** `Rscript tools/run_tests.R` gives 0 failed. The walkthrough check
"menu opens How the models were fitted" passes.

### F2 (major, verified). "Hierarchical binomial model" misdescribes the tumor fit

**Finding.** `fit/stan/di.stan` and `dd.stan` give the model as follows:
- one baseline per fitting unit, `c ~ normal(0, prior_c_sd)`, independent with a fixed prior, so no
  common (hyper) distribution;
- an arm-level random effect `sigma * z`;
- a binomial likelihood with tumor control `exp(-exp(eta))`, the Poisson tumor-control form
  (`di.stan:56–62`, `dd.stan:59–66`).

The baselines are not hierarchical. The EARLY fit is hierarchical: its trial baselines are drawn
from a common distribution (`early_fit_in_brief.md`). Calling both models "hierarchical" would be a
factual error on the page that exists to explain them.

**Fix.** Tumor wording: "a Bayesian binomial model of the number of patients with tumor control in
each arm, with the Poisson tumor-control form, a separate baseline for each trial (weak normal prior)
and an arm-level random effect for extra between-arm spread". Keep "hierarchical" for the early fit only.

**Test.** `test_v461_fitting_page.R::test_tumor_section_does_not_call_the_fit_hierarchical`. It
fails if the tumor section of `how_fitted.md` contains "hierarchical".

**Verification.** The test passes. A manual read of the tumor section against `di.stan` lines 56–66.

### F3 (major, verified). Page content in `docs/` would not reach the deployed site

**Finding.**
- `.rscignore` excludes `docs/`, so a page that read `docs/refit_explainer/*.md` at runtime would work
  locally and fail on shinyapps.io.
- `tools/deploy.R` and `test_deploy_preflight.R:54–120` list each served markdown file by name
  (`about.md`, `refit_in_brief.md`, `early_fit_in_brief.md`).

**Fix.**
- Put the page text in `how_fitted.md` at the app root, rendered with `markdown::mark()` as the briefs
  already are.
- Add the file to the deploy need-list and to the four preflight `need <-` vectors.

**Test.** `test_deploy_preflight.R::test_bundle_without_how_fitted_md_is_refused`, modeled on the
existing `early_fit_in_brief.md` case at line 120. To see it fail, remove the file from the need-list.

**Verification.** The preflight tests pass. After deploy, the live walkthrough opens the page and
finds its heading.

### F4 (major, verified). A fifth copy of the fitted numbers invites drift

**Finding.**
- The fitted tumor values already appear in four places: `refit_in_brief.md`, `about.md:125–133`,
  the explainer, and `README.md`.
- They appear on two scales: README "λ/α 0.63 Gy per day", and the brief "0.53 Gy EQD2 per day".
  Both are correct, but a reader comparing them sees a conflict.
- A new page restating the numbers adds a fifth copy that no test ties to the draws.

**Fix.**
- The new page holds the METHOD: likelihood, priors, sampler settings, the quality gate, what prior
  dominance means, and how draws become Refit support and Uncertainty.
- For results it embeds the two existing briefs verbatim, as the app already does through `REFIT_BRIEF`
  and `EARLY_BRIEF`, so each number keeps one source.
- The priors table is the one place with numbers. Test it against the Stan files.

**Test.** `test_v461_fitting_page.R::test_prior_table_matches_stan`. It parses the prior lines of
`fit/stan/di.stan`, `dd.stan` and the early-fit Stan model and compares each to the table. To see it
fail, change `normal(28, 10)` to `normal(28, 12)` in the table.

**Verification.** The test passes. `grep -c "0\.53\|0\.63" how_fitted.md` returns 0.

### F5 (minor, inferred). Jargon on a clinician-facing page

**Finding.** The page will need prior, posterior, likelihood, overdispersion, R-hat, ESS, E-BFMI and
MCMC. The reader-first rules (requirements 4 and 7) and the acronym rule require a gloss at first use.

**Fix.**
- Give each term one plain clause at first use, for example: prior "what was assumed about a parameter
  before seeing the trials"; posterior "the range of values compatible with both the prior and the
  trials".
- Spell out ESS and MCMC at first use.
- Put a three-sentence "in short" paragraph at the top, so a clinician can stop there.

**Test.** `python3 ~/.claude/reader-first/reader_first_audit.py --document how_fitted.md --config
reader_first.yaml --strict`, with the method terms added to the inventory. To see it fail, delete
one gloss.

**Verification.** Exit 0 with a non-empty inventory.

### F6 (minor, verified). Version and menu ripple

**Finding.** These items move with the change:
- `APP_VERSION` is pinned to "4.6.0" in `test_v27_colleague_review.R:259` and `test_v31_stage_e.R:67`.
- The walkthrough asserts six menu links (`walkthrough_v4.R:159`).
- `V4_MENU` and `V4_MENU_NOTES` (`app.R:409–416`).
- `DESCRIPTION` and the README "live" line.
- The v4 user guide's menu description.

**Fix.** Ripple items R1–R5, applied in the same change:
- R1: version 4.6.1;
- R2: walkthrough expects 7 links and opens the new page;
- R3: menu note "Bayesian refits: data, priors, checks";
- R4: README;
- R5: user guide.

**Test.** The existing version and menu-count assertions turn red until updated. That is the check.

**Verification.** Full suite at 0 failed, then the live walkthrough at N of N.

**Dismissed.**
- **The prior tables in explainer 5.2 disagree with the Stan code:** no. λ DI lognormal(log 0.1, 1.5),
  Tk normal(28, 10), λ DD lognormal(log 0.5, 1.5), C lognormal(log 10, 1) and α lognormal(log 0.2, 1)
  all match `di.stan` and `dd.stan`.
- **A menu entry crowds the phone layout:** no evidence. The menu is a sidebar list, and the walkthrough
  already has a phone check (`20_phone_main`), so add the new page to it.

## 3. Findings on Issue 2 (late-toxicity refit)

### L1 (major, verified for the project data; screen only for the PDFs). What exists

**Finding.**
- `data/early_tox/outcomes.csv` has 68 rows, all mucositis endpoint classes.
- The tumor data hold local control only.
- No late count has been extracted, which confirms the plan.
- On the claim that the trials report late effects: 96 of 120 downloaded PDFs mention late-effect
  terms at least 3 times.
- Several are dedicated late-effect reports:
  - Fu 1995, RTOG 83-13 late effects (113 hits);
  - Mortensen 2012, DAHANCA (53);
  - HYPNO, Bentzen 2026 (45);
  - Horiot 1997, EORTC 22851 (32);
  - Denham 1999, consequential late reactions in TROG 91.01 (31);
  - Dische 1997, CHART (23);
  - Bernier 1998 (18);
  - Cox 1991, interfraction interval and late effects (19);
  - Bentzen 1999, repair half-times (15).
- A term count is not an extractable per-arm count.

**Fix.** Feasibility step output: one table per trial giving
- endpoint as written, scale (RTOG/EORTC, LENT-SOMA, CTCAE), and site;
- per-arm numerator and denominator;
- crude or actuarial, time point, and the location in the paper.

Extract twice under the early-fit `extraction/BRIEF.md` rules.

**Test.** Reuse the early-fit builder tests, made specific to late endpoints, for example a test that
only late rows can confirm a count, mirroring `test_only_mucositis_rows_confirm_a_count`.

**Verification.** The feasibility table exists, and every count has two agreeing extractions or is
labeled sensitivity-only.

### L2 (major, derived). A late refit changes the main list only if late α/β falls to about 2 Gy or below

**Finding.** Highest allowed dose per fraction, from the app's `compute_all()`, Fowler early model
(the default), margins −5/+5/+5. The α/β 3 rows reproduce the early-toxicity plan table exactly.

| Fractions | Tumor needs ≥ | Early allows ≤ | Late ≤ at α/β 3 | at α/β 2 | at α/β 1.5 |
|---|---|---|---|---|---|
| 10 | 4.40 | 4.25 | 4.70 | 4.45 | 4.30 |
| 12 | 3.90 | 3.85 | 4.20 | 4.00 | 3.85 |
| 15 | 3.30 | 3.30 | 3.65 | 3.50 | 3.40 |
| 20 | 2.65 | 2.80 | 3.00 | 2.90 | 2.85 |
| 25 | 2.30 | 2.45 | 2.60 | 2.50 | 2.50 |
| 30 | 2.05 | 2.25 | 2.25 | 2.25 | 2.20 |

- At the current α/β of 3 Gy, late binds nowhere below 30 fractions.
- It ties with early at 12 and 16 fractions only at α/β 1.5 Gy, and binds at 30 fractions there.
- So a refit moves the main list only if its posterior sits near or below about 2 Gy.
- It would still replace the assumed late prior (95% 2–4.5 Gy, `R/uncertainty.R:269`) on the
  Uncertainty page with a data-based one, which is the stronger reason to do it.
- My script did not apply the research early model's relative calibration (`early_cal`), so its
  research-early rows are invalid and are not reported. That case is NOT computed.

**Fix.** Make this the first step of the feasibility work: a decision-value table from the app's main
page engine (`home()`), for both early models, across late α/β 1–5 Gy. It shows before any extraction
what a late fit could change.

**Test.**
`tests/python/test_late_tox_step0.py::test_alpha_beta_3_reproduces_the_early_plan_table`. It asserts
the α/β 3, Fowler rows above. To see it fail, set α/β to 4.

**Verification.** The table is printed for both early models. The research-early rows come from the
`home()` path, not `compute_all()`.

### L3 (major, judgment, not computed). Incomplete repair biases a complete-repair late fit

**Finding.**
- The app's late model is LQ with complete repair and no time factor (`R/bed_models.R:224–225`).
- Many informative arms give 2–3 fractions a day: EORTC 22791, RTOG 9003, CHART, Tallari and PMH.
- Late tissues repair slowly, and the project holds Cox 1991 and Bentzen 1999 on exactly this.
- **Mechanism (direction is judgment):** unmodeled incomplete repair makes multi-daily small-fraction
  arms more toxic than complete-repair LQ predicts. The fit can explain that only by giving small
  fractions less sparing, which pushes the estimated late α/β UP.
- A too-high late α/β understates the late cost of large fractions, which is the anti-conservative
  direction for hypofractionation.

**Fix.**
- Fit with the same incomplete-repair term Q as the early fit, with a stated half-time prior and a
  sensitivity grid.
- Report α/β and half-time jointly, because they trade off.
- Applying the result to once-daily schedules with complete repair is acceptable, since repair is
  nearly complete over 24 hours. State this.

**Test.** A simulation check: generate counts from a known α/β (for example 2.5) with a 4-hour
half-time, fit with complete repair and with Q, and assert that the Q fit recovers α/β within its
interval while the complete-repair fit is biased upward. This is the test that shows the direction.

**Verification.** The simulation's bias table, committed with the fit record.

### L4 (major, judgment). Late endpoints do not share a denominator or time scale

**Finding.**
- Late effects are counted among survivors. Arms differ in survival, scale, site and time point
  (2–5 years), and in crude versus actuarial estimates.
- A per-arm binomial on the randomized denominator understates the risk in arms with more deaths.

**Fix.** Prespecify:
- a primary endpoint class, for example grade ≥3 late mucosal or any-site RTOG/EORTC;
- a time window;
- the effective-n rule for actuarial rates, as the tumor fit already uses for Kaplan-Meier rates.

Run crude-only and actuarial-only fits as sensitivities.

**Test.** A builder test that refuses an outcome row with no time point or no crude/actuarial label.

**Verification.** The endpoint-mapping table is approved by Igor before any fit.

### L5 (major, judgment). Consequential late effects can masquerade as a fractionation effect

**Finding.**
- Severe acute reactions in accelerated arms can leave late damage.
- The Denham 1999 PDF on hand addresses this directly for TROG.
- A pure late-LQ fit attributes such excess to fraction size.
- Under that model, arms with the same n and d but different T carry no information on late α/β
  (DAHANCA 6&7, CAIR and other 5-versus-6- or 7-a-week trials). That makes them a free falsification
  check: under the model, their late rates should not differ.

**Fix.**
- Use same-dose, different-time contrasts as a model check, not as fitting data for α/β.
- Flag consequential endpoints (mucosal ulceration and necrosis) separately.
- Run a sensitivity fit with early BED as a covariate.

**Test.** The posterior predictive check on the same-dose contrasts reports observed against
predicted log odds ratios.

**Verification.** The check is reported in the fit record, and a material violation blocks default
promotion.

### L6 (major, verified). Prior double-counting

**Finding.**
- The current late prior cites Stuschke and Thames 1999 (`R/uncertainty.R:269`).
- `about.md:194` says the 3.3–4.0 Gy estimate comes from "one hyperfractionation trial".
- Stuschke 1999 analyses the randomized hyperfractionation trials, which a late refit would include.
- Using it as the prior and fitting that trial's late data counts the trial twice. The tumor refit
  avoided exactly this (explainer 5.2).

**Fix.** Use a weakly informative late α/β prior that does not come from these trials, for example
lognormal centered at 3 Gy with 95% 1–10 Gy. Use Stuschke only as an external comparison, or drop its
trial from the fit if its estimate is used as the prior.

**Test.** `test_prior_provenance` in `sim_provenance.yaml` / `evidence_ledger.yaml` records the prior
source and that the source's trials are not in the fit list. To see it fail, add that trial.

**Verification.** `python3 ~/.claude/simulation-audit/sim_provenance.py check --strict` passes.

### L7 (minor, judgment). Identification of late α/β

**Finding.**
- Under LQ without time, late α/β is informed only by within-trial contrasts that differ in dose per
  fraction, roughly 1.1–2.75 Gy.
- The hyperfractionation trials give most of that leverage. HYPNO, Moon 2014, Liu 2023 and
  Kachhwaha 2021 (hypofractionation, if usable) give the rest.
- Chemotherapy and technique era cancel within a trial only if both arms share them.

**Fix.** Run the early fit's step 4 identification check (rank and condition number of the contrast
design) on the late feasibility table before any Stan run. Keep the planned-chemotherapy indicator.

**Test.** Reuse `test_early_tox_step4`-style tests on the late design matrix.

**Verification.** The identification record states whether α/β is identified, and if not, the fit is
labeled prior-driven.

## 4. Revised plan (version 2) and implementation order

**Decisions for Igor (they block the steps named):**
- **D1 (blocks 1).** Approve F1–F6 as the revised Issue 1 plan, one destination under one label.
  Recommended: yes.
- **D2 (blocks 1).** Page label "How the models were fitted (Bayesian)". Recommended: yes.
- **D3 (blocks 6).** Start the late feasibility after 4.6.1. Recommended: yes, on a separate branch
  from `v4`.
- **D4 (blocks 8).** Primary late endpoint class and time window (L4). Recommended: grade ≥3 late
  mucosal or any-site RTOG/EORTC at 2–5 years, with a crude/actuarial sensitivity.

**Sequence:**
1. Claude: write `how_fitted.md`, covering method, prior table and the two embedded briefs (F2, F4,
   F5). Gate: the reader-first audit exits 0. About 1 hour.
2. Claude: wire the menu entry, the Details button and the three retargeted links. Add the file to
   deploy and preflight (F1, F3). Gate: the new tests fail before the wiring and pass after.
3. Claude: ripple R1–R5 (F6). Gate: `Rscript tools/run_tests.R` at 0 failed on one revision.
4. Claude: local walkthrough, including the phone view. Gate: N of N. Then FIXES row and development
   record, then commit and push to `v4`.
5. Igor approves deploy, then `Rscript tools/deploy.R`, then the live walkthrough. Gate: N of N, and
   README records the release.
6. Claude: late step 0, the decision-value table (L2). Gate: the α/β 3 rows match. About 1 hour.
7. Claude: late feasibility. Screen the PDFs, then double extraction (L1). Igor may need to supply
   PDFs for trials not on hand. Gate: the feasibility table.
8. Igor: D4. Then Claude: endpoint mapping (L4), identification check (L7), and priors with provenance
   (L6). Gate: Igor approves a written late-fit plan. No Stan run before this gate.
9. Fit with repair (L3), simulation bias check, same-dose falsification check (L5), and
   leave-one-family-out validation. This follows the early-fit steps 5–6. Research option only, as
   the early fit was.

**Stopping condition.** When steps 1–5 pass, Issue 1 is closed. Reopen only for a failed gate or a
user report. Issue 2 proceeds no further than step 8 without Igor's written approval of the late-fit plan.
