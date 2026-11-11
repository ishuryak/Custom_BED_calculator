# Early-toxicity plan, step 1: decisions and target definition

Written 2026-11-10. Applies to plan version 2 (`2026_11_10_early_toxicity_plan.md`). Docs only. No
app code, data or published numbers change, and no refit has been run.

## 1. Decisions taken

Igor Shuryak answered the open decisions on 2026-11-10 with "all as recommended".

| ID | Decision | Consequence |
|---|---|---|
| D1 | Primary endpoint group: maximum grade ≥ 3 acute mucositis, with the scale and its version recorded per arm. Confluent mucositis on the Dische scale is mapped into this group in step 3, with the mapping shown | Step 3 writes the scale-by-scale mapping, and Igor approves it |
| D2 | Primary clinical context: broad-field definitive radiotherapy, with or without concurrent chemotherapy, for locally advanced HNSCC | Section 2 below |
| D3 | Priors for the clinical benchmark are the review's proposals (plan section 6, step 5) | Fixed before fitting. Prior predictive checks run in step 4 |
| D4 | The "+5" margin stays a numerical BED increment at the sampled α/β. An EQD2 version is a sensitivity analysis | Step 7 |
| D5 | No clinical tolerance δ by default. The risk difference is shown without pass/fail labels | Step 7 |
| D6 | Denham 1995 was not obtainable. Its α/β of 9.3 Gy (5.8–17.9, 95% CI, patchy mucositis) is cited as reported by Kaanders 1999, which is labeled a secondary source. Fowler 2003 (10 Gy) stays the center. Whether the hyperfractionation trials inform α/β is checked in step 4 | The α/β sensitivity prior stays lognormal(log 10, 0.354), labeled analyst-chosen |
| D7 | Route: follow the plan's gates (steps 1–4 before the first refit). No exploratory pilot fit | Step 5 starts only after steps 1–4 pass |

Notes on D6:
- Kaanders 1999 reports "an α/β ratio of 9.3 Gy (5.8–17.9, 95% CI) with patchy mucositis as the
  endpoint" (Section 3, read in the local copy). It is a quotation of Denham, not an independent
  estimate.
- Fowler 2003 quotes the range from Kaanders as "5.8–16.4 Gy". That does not match Kaanders' text,
  so Fowler is not cited for the interval.
- Rezvani 1991 and Chougule & Supe 1993 (both in the Joiner table) were not obtained. They are not
  needed under this decision.

## 2. Clinical contexts

Contexts are kept separate. The model is fitted on the primary context only. The other contexts
are used as labeled contextual checks or not at all.

| Context | Definition | Role |
|---|---|---|
| C1 (primary) | Definitive radiotherapy to the primary and neck (broad field), with or without concurrent chemotherapy, for locally advanced squamous cell carcinoma of the oral cavity, oropharynx, hypopharynx or larynx | Fitting, validation and the app's main list |
| C2 | Early glottic cancer with small fields (for example KROG 0201, Yamazaki, DAHANCA 6, RTOG 9512) | Excluded from fitting. Field size changes mucosal exposure |
| C3 | Postoperative radiotherapy (for example Ang 2001, Awwad, Sanguineti, p-CAIR) | Excluded from the primary fit. A sensitivity analysis is possible if endpoints map |
| C4 | Nasopharyngeal carcinoma | Excluded. Different fields and systemic therapy |
| C5 | Reduced-volume, favorable-risk HPV-positive oropharynx (HYHOPE) | Contextual check only |
| C6 | Palliative, split-course or frail elderly schedules (for example ELAN-RT, Ortholan 2026) | Excluded from the primary fit |

**Systemic therapy within C1.** The planned systemic therapy is part of each arm's definition
(plan R7). Most randomized fractionation trials were radiotherapy alone, while the app's reference
setting and HYPNO are chemoradiotherapy. Using a dose-to-risk slope that was learned mainly from
radiotherapy-alone contrasts for chemoradiotherapy is a transportability assumption. It is stated
in the app text and checked against the chemoradiotherapy trials in C1 (for example GORTEC 99-02,
RTOG 0129 and HYPNO). No chemotherapy-by-schedule interaction is fitted unless the contrasts can
inform it (plan step 5).

## 3. Endpoint

- **Primary:** maximum grade ≥ 3 acute mucositis per patient (D1), assessed weekly or more often
  during radiotherapy and for at least 4 weeks after its end.
- **Window rule.** Short courses can peak after the end of treatment. An arm whose follow-up of
  acute reactions stops at the end of radiotherapy, or whose window is not reported, enters only a
  limited-ascertainment sensitivity analysis (plan step 3). The 4-week minimum is Claude's proposal
  and is part of the mapping that Igor approves in step 3.
- **Not pooled with the primary:** prolonged intolerability, duration of confluent mucositis,
  feeding-tube use, unplanned breaks, dose-limiting toxicity (plan step 3).
- **Unit:** patients with the event among evaluable patients, per randomized arm. A percentage is
  converted to a count only when the evaluable denominator is stated and the count is unique
  (plan step 2 rules).

## 4. Reference regimen and exposure

- **Reference regimen:** 70 Gy in 35 fractions over 47 inclusive days (46 elapsed), as in the app.
- **Exposure:** early BED from the single repair-aware evaluator (plan R8). The legacy evaluator
  reproduces version 4.5.0.
- **Dose:** the planned (prescribed) dose and time define the arm. Delivered dose and time are
  recorded but not adjusted for (plan R7).

## 5. Planned predictions and their comparators

Gate for step 1: every planned prediction has a stated context and comparator.

| ID | Prediction | Context | Comparator |
|---|---|---|---|
| P1 | Main-list support: share of joint draws in which a candidate's early BED is within the margin of the reference | C1 | 70 Gy / 35 fractions / 47 days, same draw (plan R9) |
| P2 | Risk difference for a candidate, shown without labels (D5) | C1, with a declared baseline (plan R5) | 70/35 reference arm, same draw |
| P3 | HYPNO benchmark: predicted against observed grade ≥ 3 mucositis | C1, chemoradiotherapy | HYPNO's randomized comparator, 66 Gy / 33 fractions at six a week. 70/35 is refused (plan test 12) |
| P4 | Held-out trial families: new-trial and conditional-contrast prediction | C1 | Each held-out trial's own control arm (plan step 6) |
| P5 | Contextual checks (not scored): Jacinto 2018, Madhava 2006, Sanghera 2007, Jegannathen 2011 (all 55 Gy / 20 fractions with chemotherapy, single arm), Miller 2026 (multicenter hypofractionation cohort, schedule to be confirmed) and HYHOPE. Duplicate patients across the Christie reports (Jegannathen 2010, 2011 and the phase I study) are checked first | C1 single arms, and C5 | None. Reported as observed against predicted risk at a stated baseline |

Single-arm studies (P5) do not enter the randomized-contrast likelihood, because they have no
within-study comparator.

## 6. Gate check

- Every prediction P1–P5 has a context and a comparator: **pass**.
- The plan, audit, review and both exploratory folders are unchanged by this step: **pass**
  (this file and the records listed in section 8 are the only additions).
- The Fowler 2003 evaluation stays a frozen fixture. The legacy-preset requirement (69 doses in
  18 rows) is unchanged: **pass**.

## 7. Starting index for step 2 (unverified)

Found on 2026-11-10 by a text scan of the trial PDFs in `References_refit_trials/eric_2026_10_22/Part 3/`
for statements about mucositis with numbers. **Nothing here is extracted data.** Each value must be
read from the full text in step 2, with evaluable denominators, the scale version and the window.
The context column is provisional.

| Trial (report) | Contrast | Provisional context | Reported mucositis measure (from scan) |
|---|---|---|---|
| DAHANCA 6&7 (Overgaard 2003) | 6 vs 5 fractions a week | C1 (C2 glottic patients included) | Confluent mucositis, with a time-to-event figure |
| GORTEC 99-02 (Bourhis 2012) | Three strategies, with and without chemotherapy | C1, chemoradiotherapy | RTOG grade 3–4 acute mucosal toxicity, counts by arm |
| Hliniak 2002 | Accelerated vs conventional | C1 | Confluent mucositis, per arm (table) |
| Ghoshal 2008 | Concomitant boost vs conventional | C1 | Grade 3 mucositis, percentages |
| Skladowski 2006 / 2000 / 2013 (CAIR) | 7 vs 5 days a week, and CAIR vs concomitant boost | C1 | Confluent mucositis, Dische scale |
| Maciejewski 1996 | 7-day vs 5-day | C1 | Confluent mucositis lasting > 3 weeks (Dische). Duration endpoint, not primary |
| Cummings 2007 | Hyperfractionated vs conventional | C1 (confirm sites) | Grade 3–4 mucositis, percentages |
| Bernier 1998 (EORTC) | Hyperfractionated vs conventional | C1 (oropharynx) | Grade 3 at end of treatment, percentages |
| Daoud 2007 | Twice vs once daily | C1 (confirm) | Severe mucositis, percentages |
| Poddar 2017 | Accelerated vs conventional | C1 | Grade 3, percentages |
| Dobrowsky 2000 (Vienna) | V-CHART ± mitomycin vs conventional | C1 | Grade 3–4, percentages |
| Shivananjappa 2023, Tandon 2018, Mani 2022, Sikdar 2026 | IMRT schedules | C1 (confirm) | Grade 3, percentages |
| Alam 2016 | Two schedules with chemotherapy | C1 (confirm) | Confluent mucositis, percentages |
| RTOG 9003 (Fu 2000), EORTC 22791 (Horiot 1992), EORTC 22851 (Horiot 1997), RTOG 0129 (Nguyen-Tan 2014), CHART (Dische 1997), TROG 91.01 (Poulsen 2001), Jackson 1997 | Various | C1 | No numeric statement found in the text scan. Check tables and figures in step 2 |
| Sanguineti 2005, Suwinski 2005 (p-CAIR) | Accelerated vs conventional | C3 | Confluent mucositis, percentages |
| Li 2006 | Dose escalation | C4 | Grade 3, percentages |

Also to be screened in step 2: HYPNO (Bentzen 2026), the Pardo-Montero 2021 schedule labels (the
historical classifier, plan step 5.3), and Moon 2014 (KROG 0201, C2).

## 8. Records changed in this step

- This file (new).
- `2026_11_10_early_toxicity_plan.md`, section 9: D1–D6 marked decided, with D6 replaced as above.
- `2026_11_10_early_toxicity_handoff.md`: status and next step.
- `References/README.md`: entries for Kaanders 1999, Sanghera 2007, Jegannathen 2011 and Lester
  2024, and the Denham 1995 status.

## 9. Next: step 2

Build `data/early_tox/` (plan step 2), starting with the C1 trials in section 7. Claude extracts.
Igor spot-checks. Every value is read from a downloaded full text, and missing PDFs are reported by
name. Gate: the R1 selection fixtures reproduce, every denominator is explicit, the primary analysis
has no unresolved duplicate evidence, and no evidence family crosses a validation split.
