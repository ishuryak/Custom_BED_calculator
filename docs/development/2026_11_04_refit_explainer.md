# 2026-11-04: refit explainer, its audit, and version 4.3.1

## What and why

Igor Shuryak asked for an explainer of the tumor-model refit behind version 4.3.0. It was to answer five questions. How was the refit done? Which trial arms differ from the 2018 paper? How does each trial get its own result? How do DI and DD compare by AICc? How and why is the fit Bayesian? He then asked for an audit with fixes, and finally for the explainer to be committed and linked from the site so a reader can open it.

- Source: `docs/refit_explainer/HNSCC_BED_v4_refit_explainer_2026_11_04.md`.
- Served copies, built with pandoc 3.8.3: `www/refit_explainer.html` (self-contained, equations as MathML, no external scripts) and `www/refit_explainer.docx`. Shiny serves `www/` at the site root.
- Links: one line on the main page under "read the sources", and a section at the top of About.
- Version 4.3.1. No equation, parameter, draws file or selection rule changed.

## The AICc numbers are new

The refit compared models only by leave-one-trial-out prediction. To answer the AICc question, `docs/refit_explainer/aicc_check.py` repeats the 2018-style maximum-likelihood comparison with one intercept per fitting unit. It uses a binomial likelihood, α/β 10 Gy, and *g* fixed at 0.0154 per day or estimated. Its outputs are `main_fixed.json`, `main_free.json` and `t2018_fixed.json`.

| Data | ΔAICc (DI − DD) | QAICc |
|---|---|---|
| 2018 Table 1 | +3.9 | +1.1 |
| Current main set | −9.0 | −4.4 |
| Current main set, *g* estimated | −8.4 | −4.0 |

The 2018 paper reported +40. The script reproduces the earlier `fit/mle_check.py` log-likelihoods on the 2018 data (4713.10 and 4711.14). In the *g*-free DD fit, *g* sits at its 0.021 bound.

## Audit (three independent read-only verifiers plus text checks)

About 650 values were compared with `data/refit_trials_2026_10_25c_main.csv`, `data/shuryak2018_table1.csv`, the `fit/output_refit_2026_10_25c_*` summaries, the Stan files, `fit/gate.py`, `R/v4_home.R` and the 2018 paper's text. Every equation (lines for TCP, DI, DD, the hinge, BED~LQ~ and effective BED) matches the Stan code. The full text diff is `2026_11_04_refit_explainer/explainer_audit_2026_11_04.patch`.

**Found and fixed:**

1. **Stale calibrated-DD values.** The DD column of the sensitivity table held values from before the 2026-10-26 effective-sample-size repair. Seven of nine rows were affected, by at most 0.1 Gy₁₀ (for example prior sd 10: −3.3 (−6.1, −0.7) → −3.4 (−6.1 to −0.8)). Now pinned by a test.
2. **Missing changes in the 2018 comparison table.** The table of values that differ from 2018 omitted two trials: CAIR vs concomitant boost (65.9/66.3 → 63/65) and Jackson (43.9/48.8 → 44.3/49.1). "Pooled" for CAIR is now marked as an inference.
3. **Leave-one-out description.** Two corrections:
   - the intercept is fixed per fitting unit, from the first arm's observed rate (a plug-in score, not a full predictive density)
   - the DD there is the free DD, with 175 of its 415 divergences in one fold.
   Three of DI's 13 wins are ties below 0.05.
4. **Drop-one influence.** The run used the 2026-10-22 set (47 arms, ARTSCAN included, all-unit baseline −2.4), not the current data. This is now stated.
5. **Rounding of observed rates.** Rounded counts move observed rates by up to 1.4 points, not "a fraction of a point".
6. **Omissions now stated:**
   - the prior on the DD onset constant *C*
   - the gate's tree-depth check
   - the calibrated-DD sampler settings
   - "Refit support" being DI or DD, within 5 Gy₁₀.
7. **AICc wording and explanations.**
   - The QAICc wording is corrected: it uses the larger of the two models' χ²/df.
   - A second explanation for the published 40-unit DD advantage was added. The 2018 DI fit hit its clonogen lower limit of 10, according to that paper's text.
8. **DI vs DD under acceleration.** The difference is 3–4 Gy₁₀, not 2–3.
9. **Writing.**
   - Undefined abbreviations spelled out: QAICc, RTOG 9003 arm codes, PGIMER, LQ, DOI and SE.
   - The Summary glosses calibrated DD, Gy₁₀, effective tumor BED and overdispersion.
   - The double negative "loses −2.0" was removed, long sentences were split, and dated internal run labels were removed.

**Checked and not changed:**

- all 45 posterior predictive rows
- all 19 leave-one-out rows
- the AICc table
- the posterior table
- the schedule table (all 30 cells, 91.9% order agreement over 55 pairs)
- the counts of papers, trials, arms and patients
- all 19 PMIDs
- the data SHA-256.

**Noted, not acted on:**

- `fit/output_refit_2026_10_25c_main/loso_*.json` carry no identity record. The current `fit_refit.py --loto` would refuse to rerun without one, so they come from an earlier code version. Their structure matches the 25c data.
- The free DD posterior's onset BED *C*/α is far above any delivered schedule for many draws.

## Known limitation of the tracked Word file

The tracked-changes Word version (kept outside the repository with the clean version) shows text changes. Changed tables pass through without markup. The table changes are those in items 1, 2 and 6 above.
