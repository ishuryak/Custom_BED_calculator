# Early-toxicity plan, step 2: data and provenance tables

Written 2026-11-10. Applies to plan version 2 and the step 1 record
(`2026_11_10_early_toxicity_step1_target_definition.md`). Data tables, builders and tests only. No
app code and no published number changes, and no model has been fitted.

## 1. What was built

All tables are in `data/early_tox/` and are rebuilt by two scripts.

| Table | Rows | Built by | Content |
|---|---:|---|---|
| `labels.csv` | 116 | `tools/early_tox/build_historical.py` | Pardo-Montero 2021 schedule labels, with course scope and the primary and sensitivity selections |
| `schedules.csv`, `fractions.csv` | 97, 4061 | same | Historical schedules, identity from the full dose and time vectors |
| `reports.csv` | 1 | same | Pardo-Montero 2021 (DOI 10.1002/mp.14834, PMID 33704792) |
| `arms.csv` | 43 | `tools/early_tox/build_arms.py` | Randomized arms: role, context, planned systemic therapy, schedule pattern |
| `outcomes.csv` | 45 | same | Arm-level counts with endpoint wording, scale, window, quote and location |
| `arm_schedules.csv`, `arm_fractions.csv` | 23, 911 | same | Planned fraction doses and times for each arm (two variants for CAIR) |
| `evidence_families.csv` | 30 | same | One family per trial, listing every report of it |

Sources of the arm tables, all committed:
- `data/early_tox/source/arm_curation.csv`, `outcome_curation.csv` and `screening_log.csv`
  (Claude's curation and adjudication).
- `docs/development/2026_11_10_early_toxicity/extraction/` (the brief and the four independent
  extractions).

Plan file names are kept, with one addition: trial-arm schedules are in `arm_schedules.csv` and
`arm_fractions.csv`, so the historical tables can be rebuilt without the trial sources.

## 2. How the counts were extracted

- 36 trial reports were split into two batches. Each batch was extracted twice, by two agents
  working blind to each other and reading the PDFs, with a fixed brief (`extraction/BRIEF.md`).
- The builder writes an outcome only if both extractions of its report contain the same
  denominator and event count. Quotes, locations, scales and windows are copied from the
  extraction, not retyped.
- **Agreement.** All 45 counts used are in both extractions of their batch. The disagreements
  found were secondary endpoints recorded by only one extractor, plus two percentage conversions
  (section 4).

## 3. Result

| Role | Families | Arms | Use |
|---|---:|---:|---|
| candidate (context C1, counts by arm) | 15 | 35 | Fitting candidates for step 5, subject to the step 3 endpoint mapping |
| sensitivity | 3 | 6 | IAEA-ACC (ascertainment suspect), Pinto 1991 (own scale, counts in a figure), Poddar 2017 (grade at end of treatment, stage I-II) |
| benchmark | 1 | 2 | HYPNO (plan R6), never fitted |

The candidate families are GORTEC 99-02, RTOG 9003, RTOG 0129, EORTC 22791, EORTC 22851, the
Princess Margaret hyperfractionation trial, Ghoshal 2008, CAIR (Gliwice), Vienna V-CHART, ORO 93-01,
Bartelink 2002, Tandon 2018, Sikdar 2026, Shivananjappa 2023 and Mani 2022. Their 37 candidate
outcome rows cover 4,599 patients, counting the CAIR arms once.

**Endpoints are not yet pooled.** `endpoint_class` keeps what each trial measured: grade 3-4 on
RTOG/EORTC or NCI CTC, grade 3 on CTCAE v4, top-grade diffuse mucositis on the EORTC objective
scale, confluent mucositis, and two CAIR thresholds on a Dische-based scale. Step 3 maps them to
the primary group (D1), and Igor approves the mapping.

## 4. Adjudications

| Item | Decision | Reason |
|---|---|---|
| Skladowski 2006 (94% vs 53%) and CAIR-2 2013 (89% vs 86%) | Not derivable. Extraction A2 derived counts, A1 refused, and A1 is followed | No acute-toxicity denominator is stated. 94% fits both 47/50 and 48/51, and CAIR-2 notes about 25% missing scoring |
| Ghoshal conventional arm | 27/142 (Table 4) | The Discussion says 23/142. The conflict is kept in the note |
| Olmi arm C | 24/50, grade 3+4 (Table 7) | The text's 44% is grade 3 alone |
| CAIR arms | Two outcomes per arm (approximate EORTC grade III-IV, and severe confluent) | Different thresholds. Step 3 chooses |
| CAIR schedule | Two variants per arm, 2.0 Gy era and 1.8 Gy era, weighted by patient numbers | The fraction size changed on 1995-01-01 in both arms. Era counts come from Maciejewski 1996 and the arm totals |
| Duplicate reports | Bernier 1998, Skladowski 2006, Maciejewski 1996, Dobrowsky 1998 and 2000 are refused by the builder | Same patients as an included report |
| Percentage-only reports | DAHANCA 6&7, Hliniak, CHART, TROG 91.01 and CAIR-2 are not derivable | No unique count. The screening log records each one |
| Excluded | Jackson (skin or mucosa), RTOG 8313 (any organ), Alam (weekly snapshot), Daoud (nasopharynx), EORTC 22811, ARTSCAN 2011 and Saunders 2010 (no usable acute data) | Screening log |

**Schedules checked against the papers' own day counts:**
- RTOG 0129: 47 days standard and 40 days accelerated, matching its footnote.
- HYPNO: 26 against 38 days, matching its stated 12-day difference.
- Vienna V-CHART: 17 days.

**Two schedule items are open and flagged in `arms.csv`:**
- **EORTC 22851 accelerated arm.** The second course is "27 fractions and 17 days" at three a day,
  which cannot be 9 consecutive treatment days. The nominal pattern uses 9 weekdays from day 21.
  A sensitivity variant is set up in step 4.
- **Princess Margaret trial.** The paper gives 28 days. The weekday pattern gives 26 inclusive
  days.

**Prescription and mucosal dose.** The IMRT trials report the dose to the high-dose target. Their
mucosal dose is not reported (plan section 8, exposure layer). Mani 2022 has no contrast in the
prescription schedule.

## 5. Gate for step 2

| Gate item | Status |
|---|---|
| R1 selection fixtures reproduce (116 / 115 / 111 / 107 / 103) | Pass (`test_early_tox_historical.py`) |
| Every denominator explicit | Pass. Each outcome has its denominator and its type |
| No unresolved duplicate evidence among arms | Pass. The builder refuses screened duplicates (`test_early_tox_arms.py`) |
| No evidence family crosses a validation split | **Open for the historical labels.** The Pardo-Montero records (K, F and H references) are not yet mapped to trials, so a label from, say, CHART cannot yet be held out with CHART's family. This is needed before step 6, not before steps 3-4, because the two likelihoods are fitted separately (plan R5) |
| Igor spot-check | **Pass, 2026-11-11, done by Claude against the PDFs at Igor's request** (step 5 record, section 1). All 5 items match. Ghoshal: 27/142 confirmed (table grades sum to 142 and 27/142 = 19%), the Discussion's 23/142 is a typographical error |

## 6. Tests

`tests/python/test_early_tox_historical.py` (5 tests) and `tests/python/test_early_tox_arms.py`
(19 tests). Each was turned red by a targeted mutation:
- a partial course treated as complete
- a coarse schedule key
- disputed labels in the primary selection
- an unknown label accepted
- one extraction accepted as enough
- Saturday counted as a W5 day
- the gap off by one
- the same-day interval ignored
- screening ignored
- an unknown pattern token accepted

Python suite: 232 passed. Strict stale-constant audit: 0 errors.

## 7. Next

- Step 3 (endpoint mapping, Igor approves).
- Then step 4: identification, including the Tk profile and whether the hyperfractionation
  contrasts (RTOG 9003, EORTC 22791, PMH, Pinto) inform α/β (D6).
- The historical-label crosswalk to trial families before step 6.
