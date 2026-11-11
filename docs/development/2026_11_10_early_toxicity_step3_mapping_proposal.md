# Early-toxicity plan, step 3: endpoint mapping (PROPOSAL, awaiting approval)

Written 2026-11-10. Plan step 3 is "Claude, 1 day, Igor approves the mapping". This file is the
proposal. Nothing in it has been applied to the data, and no model has been fitted.

## 1. The question

Decision D1 fixed the primary endpoint group as maximum grade ≥ 3 acute mucositis, with confluent
mucositis on the Dische scale mapped into it. The 15 candidate families measured this in five
different ways, over different assessment windows. The mapping decides which arm counts are
pooled in the primary likelihood.

## 2. What each family measured

| Family (arms) | Scale and threshold as reported | Assessment window | Proposed group |
|---|---|---|---|
| GORTEC 99-02 (3) | RTOG acute, grade 3-4 | to 3 months after the end of treatment | G1 |
| RTOG 9003 (4) | RTOG, worst mucous membrane grade 3+4 | within 90 days of the start | G1 |
| RTOG 0129 (2) | NCI CTC v2, grade 3-5 | within 90 days of the start | G1 |
| EORTC 22851 (2) | EORTC/RTOG objective mucosal reaction, grade 3-4 | to 6 weeks after RT | G1 |
| Princess Margaret (2) | RTOG/EORTC, grade 3-4 | about 12 weeks | G1 |
| CAIR (2) | Dische score approximated to EORTC grade III-IV | first 2 months | G1 |
| ORO 93-01 (3) | RTOG acute, grade 3+4 (Table 7) | until major reactions settled | G1 |
| Vienna V-CHART (3) | "confluent extensive mucosal reaction (grade 3)", scale not named | until normalization | G1 |
| EORTC 22791 (2) | EORTC objective scale, grade 3 = diffuse (top grade) | during treatment | G1, window W2 |
| Ghoshal (2) | RTOG/EORTC, grade 3+4 | not stated (weekly during treatment) | G1, window W2 |
| Bartelink (2) | EORTC/RTOG objective "diffuse mucositis" | not stated | G1, window W2 |
| Mani (2) | RTOG, grade 3 | not stated (weekly during RT) | G1, window W2 |
| Tandon (2) | CTCAE v4.0, grade 3 | to 12 weeks after CRT | G2 |
| Sikdar (2) | CTCAE v4.0, grade 3-4 | to 2 weeks after RT | G2, window W2 |
| Shivananjappa (2) | CTCAE v4.0, grade 3 | not stated (weekly during treatment) | G2, window W2 |

Groups:
- **G1, appearance-based.** RTOG/EORTC grade 3 (confluent), NCI CTC v2 grade 3, the EORTC
  objective top grade (diffuse) and Dische-based confluent mucositis. These grade what the
  examiner sees.
- **G2, CTCAE v4 grade 3.** This version grades mainly symptoms and function (pain interfering
  with oral intake) rather than appearance. That difference is recalled from the scale's design
  and must be confirmed against the CTCAE v4.0 document before step 5. HYPNO, the benchmark, also
  uses CTCAE v4.0.

Windows:
- **W1:** the window is stated and reaches at least 4 weeks after the end of RT (or 90 days from
  the start for 7-week courses), or follows patients until healing.
- **W2:** the window is not stated, or ends during or within 2 weeks of the end of RT. The late peak
  that matters is after short courses. Every W2 family has arms of 5 weeks or longer, where the peak
  normally falls within treatment.

Not in the primary group (step 2 roles):
- severe confluent mucositis in CAIR (Dische score above 15), which is a stricter threshold, kept
  as secondary
- IAEA-ACC, Pinto and Poddar, which are sensitivity only (limited or suspect ascertainment)

## 3. Decisions (recommendation first)

- **M1. CTCAE v4 trials (G2).** Pool them with G1, using a fixed scale-family offset in the
  baseline, so a symptom-based grade 3 is not treated as identical to a confluent one. **Recommended.**
  The offset also lets the HYPNO benchmark (CTCAE v4) be predicted on its own scale. Alternative:
  exclude G2 from the primary fit (12 families remain, and HYPNO is then predicted on a scale no
  fitted trial used).
- **M2. Window.** Primary = W1 and W2. Sensitivity analysis = W1 only (9 families: GORTEC, RTOG 9003,
  RTOG 0129, EORTC 22851, Princess Margaret, CAIR, ORO 93-01, Vienna and Tandon). **Recommended.**
  Alternative: W1 only as primary.
- **M3. CAIR threshold.** Use the approximate EORTC III-IV count (48/50 against 35/49) in the
  primary, and severe confluent (31/50 against 13/49) as secondary. **Recommended.**
- **M4. Sensitivity-only trials.** IAEA-ACC, Pinto 1991 and Poddar 2017 stay out of the primary and
  are run as one sensitivity analysis. **Recommended.**

## 4. Spot-check request (step 2 gate)

Please check these against the PDFs (paths in `data/early_tox/outcomes.csv`, column `report_file`):
1. Ghoshal 2008, Table 4 against the Discussion: 27/142 or 23/142 for the conventional arm.
2. Skladowski 2000: the two CAIR thresholds (48/50 and 35/49, 31/50 and 13/49).
3. Olmi 2003, Table 7: grade 3 + grade 4 = 8/54, 21/52 and 24/50.
4. Fu 2000 (RTOG 9003): 67/268, 110/263, 112/274 and 123/268.
5. Tandon 2018: 11/30 and 16/30.
