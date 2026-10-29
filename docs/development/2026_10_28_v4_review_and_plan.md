# Review of version 4.0.0 (commit 61ccaed) and fix plan, 2026-10-28

Four independent reviews (code correctness, radiation oncology, statistics and radiobiology, clinician usability) plus executable checks. Every number below was re-run before this record was written. No code has been changed. Items marked DECISION wait for Igor Shuryak.

## Auditor triage

| Auditor | Verdict | Reason |
|---|---|---|
| paper-code / pipeline (code correctness) | RUN | Subagent review of the v4 diff, plus suite and stale-constant audit |
| radiation-oncology-reviewer | RUN | Clinician-facing radiotherapy tool |
| statistician-review (with radiobiology) | RUN | Selection rule, margins, uncertainty, extrapolation |
| reader-first (as a usability review) | RUN | The deliverable is a clinician page, not a manuscript |
| stale-constant-audit | RUN | `--strict`: 0 errors on 61ccaed |
| biomedical-data (patient-level checks) | RUN | No patient-level data. Data are one row per published trial arm. Screenshots show app tables only. |
| methods-audit | RUN, folded into the statistician review | Plan-level questions (criteria, margins) are reviewed there |
| surprising-result-skeptic | RUN, as re-verification | Each load-bearing number was recomputed (table below) |
| manuscript-audit, evidence-audit, journal compliance, tough-reviewer | N/A | No manuscript or grant in this deliverable |
| simulation-audit | N/A | v4 adds no simulated data (the refit draws are audited in v3 records) |
| ml-practitioner-review | N/A | No learned model |
| radiation-epidemiology-review | N/A | No exposed cohort. Radiobiology is covered by the statistician review. |
| biomedical-data-audit (variable-level) | N/A | No patient dataset |
| subagent-scan | RUN | Four reviewers dispatched in parallel with fixed briefs |

Executable checks on 61ccaed: R suite 5597 expectations, 0 failed (marker written). `pytest tests/python`: 114 passed. Strict stale-constant audit: 0 errors. `.gitignore` lacks the recommended per-patient figure patterns (inherited from v3, low risk because the repository holds no patient data).

## Verified numbers (defaults, Scenario D18)

| Quantity | Value |
|---|---|
| Early cap of the main list (standard + 5) | 58.11 Gy₁₀ |
| HYPNO 55 Gy/20 early BED (v3 default early limit) | 55.87 Gy₁₀ |
| Fowler lower edge with chemotherapy | 54 Gy₁₀ |
| 20 × 2.65 Gy: ΔDI, ΔDD | −4.47, −11.06 Gy₁₀ |
| 20 × 2.75 Gy (HYPNO): ΔDI, ΔDD | −1.39, −8.67 Gy₁₀ |
| 16 × 3.20 Gy: ΔDI, ΔDD, early BED | −3.93, −10.20, 56.49 Gy₁₀ |
| 21 × 2.55 Gy: ΔDI, ΔDD | −4.98, −11.21 Gy₁₀ |
| Share of refit draws meeting the tumor condition (DI / DD / either) | 20 × 2.75: 0.95 / 0.80 / 0.95. 16 × 3.2: 0.77 / 0.62 / 0.77. 21 × 2.55: 0.36 / 0.21 / 0.37 |
| Missed by the 0.05 Gy grid | 15 × 3.330–3.335 Gy qualifies (15 fractions absent from the list) |

## Findings and proposed fixes

### A. Clinical honesty of the main page (highest priority)

A1. **"About as effective" overstates the criteria.** A 5 Gy₁₀ tumor BED loss is about 4.2 Gy EQD2 (roughly two 2-Gy fractions or a week of repopulation). Fix: lead text "Shorter schedules of 5 fractions a week. Their modeled tumor BED is at most 5 Gy₁₀ (about 4 Gy in 2-Gy fractions) below 35 × 2 Gy. Their modeled early and late side-effect BEDs are at most 5 Gy above it. These are candidates for a trial, not equivalent schedules." Units on the three criteria (Gy₁₀, Gy₁₀, Gy₃).

A2. **The DD disagreement is hidden.** Rows 16–21 and 23 pass on DI alone, while DD gives losses of 7–11 Gy₁₀. "More robust choices" overstates agreement, because DD is calibrated to DI and shares its parameters. Fix: add a column "Tumor BED change, DI / DD (Gy₁₀)" showing the range, for example "−1.4 / −8.7" at HYPNO. Replace the footnote with: "The two tumor models share their main parameters and agree by construction at 35 × 2 Gy, so agreement is not independent confirmation. Where only one model meets the condition, the other estimates a larger loss."

A3. **No extrapolation warning on the main page.** The tumor data cover 1.15–2.55 Gy per fraction. Rows 16–21 lie entirely above it. Fix: mark doses above 2.55 Gy with † "outside the fraction sizes of the trials the tumor model was fitted to".

A4. **Key caveats are below the table.** Fix: move a short caution box directly under the criteria: chemotherapy (parameters mostly from radiotherapy-alone trials, and acceleration did not help with cisplatin in RTOG 0129), HPV, late BED at the prescription dose only (not hot spots or organs), treatment gaps not modeled.

A5. **Early margin is looser than the app's own acute limits** (58.1 vs HYPNO 55.9 Gy₁₀). The top doses of 16–22 fractions exceed HYPNO's acute BED, which gave about 51% grade ≥3 mucositis with weekly cisplatin (figure from the radiation-oncology review, to be checked against the HYPNO paper before it is quoted on the page). DECISION 1.

A6. **Trial-arm column.** "HYPNO test arm" next to a 35 × 2 Gy comparison suggests HYPNO tested against it. Fix: "HYPNO test arm (trial control: 66–70 Gy in 33–35 fractions, 6 a week, with cisplatin)". HYHOPE never appears (list starts at 16 fractions), so drop it from the trial list or keep it only in the footnote.

### B. Uncertainty and parameters

B1. **List membership is fragile.** Low-dose edges meet the tumor condition in only about 40–50% of refit draws (21 × 2.55: 37%). DECISION 2.

B2. **Default is the unfitted D18 scenario**, disclosed only as "tumor scenario D18". Under the refit the list is 13–34 fractions and both models agree on most rows. DECISION 3.

### C. Consistency with the menu pages

C1. **Contradictory verdicts.** 18 × 2.9 Gy is listed on the main page, while the Report card (comparison 35 × 2, tumor threshold +1 Gy, early limit HYPNO) says "threshold +1: not met". C2. A row click does not carry the regimen to the Report card. Fix: an "Open in Report card" button in the detail card that sets the regimen (chosen dose) and the main-page rules (comparison 35 × 2 Gy, tumor threshold −5, early and late limits = comparison + 5) and states that it did so. DECISION 4.

### D. Correctness and tests

D1. **0.05 Gy grid** drops fraction numbers that qualify only in a narrower window (15 fractions at D18) and puts range edges up to 0.05 Gy inside the true edges. Fix: say "on a 0.05 Gy dose grid" on the page and in About. (Root-finding the true edges is possible but would show doses like 3.333 Gy, which are not prescribable.)
D2. **DD unavailable** (for example DD rate 0) is not announced on the main page, so every row reads "DI model only" as if DD had failed. Fix: a notice in the count line.
D3. **Tests weaker than the record says.** The preset test compares row counts only, and the empty-list test never checks the reason. Fix both, and correct the 2026-10-28 record.
D4. `home_settings` compares a character vector (via `clock_e`), so the tolerance is inert. Errs safe. Fix: compare numeric fields.
D5. Empty-list reasons reuse Report-card wording ("the comparison regimen"). Fix: main-page wording that names the standard.
D6. `about.md` "every page (page)" artifact. Report card intro "on the left" is wrong on phones.

### E. Usability

E1. Detail card renders about 1,500 px below the click. Fix: scroll to it on selection, and put "Click a row for dose-by-dose numbers" above the table.
E2. Rows saving 1–4 days (31–34 fractions at 2.05–2.20 Gy) push the useful rows down. DECISION 5.
E3. Phone: freeze the Fractions column (DT FixedColumns), shorten the model column text.
E4. Close the menu after a page is chosen on narrow screens. Darken link color for contrast (teal on white is about 2.5:1). Sticky banner covers table headers on menu pages (make it non-sticky or add scroll padding).
E5. About page opens with "What is new in version 3"-style detail after the v4 section, acceptable. Menu labels: add one-line descriptions under Search, Both-model screen and Sensitivity.

## Decisions for Igor (recommendations first)

1. Early limit: **(recommended) keep your +5 Gy rule and flag doses above HYPNO's acute BED** / cap at HYPNO's 55.9 Gy₁₀ / leave as is.
2. Uncertainty on the main page: **(recommended) add a column with the share of refit draws meeting the tumor condition, and mark doses below 80%** / footnote only / nothing.
3. Default parameters: **(recommended) keep D18 but say on the page that it is not fitted to outcome data and that the trial refit is one click away** / switch the default to the refit (also changes every v3 page's defaults).
4. "Open in Report card": **(recommended) yes, and it applies the main-page rules on the menu pages** / carry the regimen only.
5. Near-conventional rows: **(recommended) list only schedules at least 7 days shorter (drops 32–34 fractions, which save 1–3 days)** / keep all.

## Decisions taken (Igor Shuryak, 2026-10-28)

1. Keep the +5 Gy early rule and flag doses above HYPNO's acute BED.
2. Add a refit-support column (share of refit draws meeting the tumor condition) and mark doses below 80%.
3. Keep Scenario D18 as the default and disclose on the page that it is not fitted.
4. "Open in Report card" with the main-page rules (recommended default, applied with the other plain fixes).
5. Keep all rows (no minimum days saved beyond 1).

## Plan once decisions are made

Apply A1–A4, A6, C, D1–D6, E1 and E3–E5 together with the decided options. Each behavior change gets a test that fails first (A2 column, A3 marker, D2 notice, D3 test fixes, C button, decisions 1, 2 and 5). Then re-run the full suite, the browser walkthrough (extended to the new button and markers) and the strict stale-constant audit. Finally, write the FIXES row and bump the version to 4.1.0.

## Applied in version 4.1.0 (2026-10-28)

- A1, A4: new lead text with units and "candidates for a trial, not equivalent schedules". The caution box (chemotherapy, HPV, volume, gaps, organ doses, late effects at the prescription dose) sits under the criteria, above the table. RTOG 0129 is not cited because no project source supports it.
- A2: column "Tumor BED change, DI / DD (Gy₁₀)" with the range over the row's doses. Footnote says the two models share their main parameters and agree by construction at 35 × 2 Gy. The column "Tumor condition holds in" reads "both models", "DI only" or "both at X Gy, else DI only".
- A3, A5 (decision 1), A6, B1 (decision 2): a Notes column marks "outside trial data" (above 2.55 Gy, read from the refit data file), "acute > HYPNO from X Gy" (acute BED above HYPNO 55 Gy/20 under the same parameters), "support < 80% at X Gy", and the HYPNO test arm. The footnote states HYPNO's acute toxicity (50.9% grade ≥3 mucositis in 395 patients, 75% with weekly cisplatin, as in About) and that HYPNO's control arm was 66 Gy in 33 fractions at 6 a week. HYHOPE removed from the trial list (never inside the list). The 3.1 Gy asterisk is replaced by the extrapolation note.
- B1 (decision 2): column "Refit support", the share of the 4000 primary refit draws (DI and calibrated DD paired) meeting the tumor condition under either model, lowest to highest dose of the row, also per dose in the row card.
- B2 (decision 3): the settings line says D18 was chosen by the developer and is not fitted to outcome data, and that a fit to randomized trials is available as a preset.
- C1, C2 (decision 4): the row card has "Open on the Report card", which sets the regimen and the main-page rules (`v4_report_inputs()`) on every menu page and says so in a notification. With these rules the Report card gives "threshold −5: met" and no "exceeds" for a listed dose, and "exceeds" for the next dose step (tested).
- D1: the page and About say doses are checked in 0.05 Gy steps, and About gives the 15-fraction example. D2: a notice when the DD model is not available. D3: tests corrected. D4: numeric comparison of settings. D5: empty-list reasons name the standard. D6: About and Report card wording.
- E1: hint above the table, and the row card scrolls into view. Every page change opens at the top of the page (found during this pass: the Report card opened half-way down). E3: first column frozen on narrow screens (DT FixedColumns). E4: menu closes after a choice below 992 px, link color #0f6e5c (contrast 6.2:1 on white and 5.9:1 on the light grey, against 2.4:1 for the previous teal), banner no longer sticky. E5: one-line descriptions under seven menu entries.
- Decision 5: all rows kept.

Not changed: v3 page contents and defaults, equations, draws, selection rules.
