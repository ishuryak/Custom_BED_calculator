# OAR constraint library: source notes (2026-10-21)

Companion to `oar_constraints.csv` (118 rows). Every value was read from the files in
`HNSCC_BED_calculator/References_OAR/`. Nothing was entered from memory.

## Conventions used in the CSV

- `n_fractions = 0` means the source states the limit for conventional fractionation
  (1.8–2 Gy/fx). In that case `dose_per_fraction_basis` gives what the source says ("2", "1.8-2",
  "<=2", or "conventional" when the source gives no per-fraction number).
- `limit_type` is "upper" on every row, as specified. **For HYHOPE `major_deviation` rows the number
  is a "Max ≥ X" trigger.** It equals the minor-deviation upper bound, or the per-protocol bound when
  the minor level is "None". HYHOPE minor-deviation rows stated as "None" have an empty `limit_gy`.
- HYHOPE "Recommended" goals are **ranges** (for example "mean < 16-21 Gy"). `limit_gy` holds the
  range string exactly as printed. No midpoint or single value was chosen.
- Volume constraints (V50, V44, V58, V35, V70) carry "% volume, not Gy" in `limit_gy`. The metric
  column names the dose level.
- HYHOPE quotes are table cells (organ + column + cell), not running prose. All 78 were compared
  with the page images and with the extracted text. Every non-HYHOPE quote was checked by script and
  is a verbatim substring of its source (whitespace normalized). The three Kirkpatrick S47 quotes
  differ only by a de-hyphenated line break ("myelop-athy").
- Full-text `.txt` files from PMC have no page numbers, so `page` gives the section instead.

## 1. HYHOPE protocol (Moon, STU 2020-1079 v4.1, NCT04580446): 78 rows

**Pages read:** PDF pp. 15–19 (protocol pp. 10–14), sections 4.1.2–4.1.4.3. I checked the page images of
PDF pp. 18–19, and I searched the rest of the document by text.

- **OAR tables:** section 4.1.4.2 "Organ at risk (OAR) constraints". 12 fx and 15 fx mandatory are on PDF
  p. 18 (protocol p. 13). The 15 fx recommended goals and all of 20 fx are on PDF p. 19 (protocol p. 14).
  - **Mandatory organs (5):** spinal cord, spinal cord + 5 mm, brainstem, brainstem + 5 mm, mandible
    (all Max). Each has three levels: per protocol, minor deviation, major deviation.
  - **Recommended organs (10 lines, 11 metrics):** contralateral parotid (mean), contralateral
    submandibular gland (mean), cochlea (maximum and mean), oral cavity excluding PTV (mean),
    constrictors (mean), post-arytenoid and cricoid space (mean), larynx (mean), cervical esophagus
    (mean), brachial plexus (max), brain (max).
  - 5×3×3 + 11×3 = 78 rows.
- **Key values (per protocol / minor / major):**
  - **15 fx spinal cord:** Max < 37 Gy / None / Max ≥ 37 Gy (PDF p. 18).
  - **20 fx (PDF p. 19):**
    - cord: < 40 / None / ≥ 40;
    - cord + 5 mm: < 46 / < 49 / ≥ 49;
    - brainstem: < 46 / < 50 / ≥ 50;
    - brainstem + 5 mm: < 52 / None / ≥ 52;
    - mandible: < 54.6 / < 57.2 / ≥ 57.2;
    - contralateral parotid: mean < 19-23;
    - constrictors: mean < 33-40.
- **The 20 fx submandibular line carries a superscript "+"** ("Contralateral submandibular gland⁺").
  No footnote defining it was found in the extracted text of the protocol. It is unexplained.
- **How the limits were derived: NOT STATED.** The protocol gives no source constraint, no α/β, and no
  conversion method for the OAR tables. The only α/β statement is in 4.1.2 (PDF p. 15), and it
  concerns choosing the prescription:
  - "α/β = 3 for late complications, 10 for tumor/mucosa", α = 0.35 Gy⁻¹,
  - Tk = 21 d (tumor) / 7 d (mucosa), Tp = 3 d / 2.5 d,
  - comparison regimen 60 Gy/30 fx (40 days).

  The BED formula itself is an image and did not extract. The 4.1.2 tables give BED3 for the
  prescriptions: 46.5/15 = 95, 52/20 = 97, 44.4/12 = 99, versus 60/30 = 100.
- **Inference, not a protocol statement (flag for the validation fixture):** I back-calculated EQD2
  from the mandatory point limits with plain LQ. The results are not one EQD2 across 12/15/20 fx:
  - cord at α/β 2: 41.1 / 41.3 / 40.0 Gy;
  - mandible at α/β 3: 64.2 / 61.2 / 62.6 Gy;
  - brainstem + 5 mm at α/β 3: 58.7 / 59.5 / 58.2 Gy.

  A single-α/β LQ conversion from one conventional constraint therefore will **not** reproduce every
  HYHOPE cell exactly. Rounding, a different source per fractionation, or another model may explain
  this, but the protocol does not say which. The fixture needs a stated tolerance, or it should be
  limited to the cells that do agree.

## 2. HYPNO (Bentzen et al. 2026, IJROBP): 0 rows

**Pages read:** all 12 PDF pages by text, with the page image of PDF p. 4.

**Page correction:** the method appears on journal **pp. 390–391** (PDF pp. 4–5), not 389–390.
PDF p. 4 is headed "390 Bentzen et al.".

**Conversion method (verbatim):**

> "Planning dose constraints in the 33-fraction NFX arm were those employed by the Danish Head and
> Neck Cancer Group (DAHANCA), see protocol in Appendix E1. These were converted to constraints in
> the 20F HFX schedule that would produce the same equieffective dose in 2-Gy fractions, EQD2, for a
> specific value of the α/β-ratio. For all endpoints, except the spinal cord and the brain stem, α/β =
> 3 Gy was assumed. For the spinal cord and brain stem, α/β = 2 Gy was assumed. It is important to
> note that these constraints are estimates that presume the validity of the LQ model. For mean doses
> to a structure (ie, when the constraint does not refer to a point (or limited volume) dose), the LQ
> fraction-size adjustment should ideally be applied to each bin in the dose−volume histogram before
> calculating the mean dose. In most cases, however, this will be a relatively small correction with
> the schedules compared in HYPNO."

**What the file does not contain:**

- **Appendix E1 is NOT in this file.** No DAHANCA constraint values appear, and neither do the
  converted 20-fraction values. "cord", "parotid" and "mandib" occur only in the passage above or not
  at all.
- HYPNO therefore contributes the method only.
- The HYPNO schedule is 55 Gy/20 fx versus 66 Gy/33 fx at 6 fx/week, not HYHOPE's 52/20.

## 3. QUANTEC and other sources: 40 rows

| Source (file) | Read | Rows | Notes |
|---|---|---|---|
| Kirkpatrick 2010 spinal cord (PDF) | all 8 pp; Recommended Dose–Volume Limits on PDF p. 6 = S47; abstract S42 | 7 | Full cord at 2 Gy/day: 50 Gy → 0.2%, 60 → 6%, 69 → 50% myelopathy. Abstract: <1% at 54 Gy, <10% at 61 Gy (1.8–2 Gy/fx). Single-fraction 13 Gy and 3-fraction 20 Gy. **States α/β = 0.87 Gy** (95% CI 0.54–1.19). Only the 50 Gy / 0.2% row is tagged `recommended`; 60 and 69 Gy are dose-response points (`NA`). |
| Deasy 2010 salivary (txt + PDF) | section 8; PMC manuscript PDF p. 5 | 3 | One parotid ≈20 Gy, both ≈25 Gy mean, submandibular "<35 Gy to see any effect". No α/β for the constraint: it mentions ≈20 Gy for early effects in rats and a possibly low late α/β, without a value. |
| Mayo 2010 brainstem (txt) | section 8 | 4 | Entire brainstem 54 Gy; 1–10 cc to 59 Gy (≤2 Gy/fx); risk rises markedly above 64 Gy; SRS 12.5 Gy (<5%). No single α/β is stated: the isoeffect curves use 3.3, 2.5 and 2.1 Gy. |
| Rancati 2010 larynx/pharynx (txt) | A.7-8, B.5-8 | 6 | Larynx (edema): V50 ≤27%, mean ≤44 Gy. Non-involved larynx (voice): mean 40–45 Gy, max <63–66 Gy. Proximal esophagus: <65 Gy (first 2 cm) and <60 Gy (first 5 cm), which is Laurell's recommendation as cited. |
| Bhandare 2010 hearing (txt) | section 8 | 2 | Cochlea mean ≤45 Gy, or ≤35 Gy as the conservative option. α/β 3 appears only for converting the Kwong data in a figure legend, so `ab_gy_source` is left blank. |
| Marks 2010 NTCP (txt) | whole file | 0 | **The summary-table body is missing from the txt.** Only the title and footnotes extracted. No values are available. |
| Brodin 2018 (txt) | Table 1, Results, highlights | 15 | Table 1 is a QUANTEC compendium, so it is a **secondary** source: 11 rows. The post-QUANTEC narrative gives 4 rows. The study-level rows of Tables 2–5 were not entered. |
| MDACC 2017 mandible ORN (txt) | abstract and conclusion | 2 | V44 <42% and V58 <25%. Dose per fraction is not stated. The cohort mean was 32.8–33.3 fractions at 68.4–69.2 Gy, so `n_fractions` = NA. |
| Murdoch-Kinch 2008 SMG (txt) | abstract, methods, results | 1 | Mean 39 Gy threshold. The same threshold holds for nominal dose and for BED2 at **α/β 3**. |

## 4. Expected but NOT in these files

| Constraint | Status |
|---|---|
| **HYPNO Appendix E1** (DAHANCA 33-fraction limits and converted 20F limits) | **Absent.** |
| **Derivation of the HYHOPE OAR tables** (source constraint, α/β) | **Not stated** in the protocol. |
| **Pharyngeal constrictor numeric limit from QUANTEC primary (Rancati)** | **Absent.** The recommendation is qualitative: minimize V60 and reduce V50. The "mean <50 Gy, <20%" figure appears only in Brodin's Table 1 (secondary). Rancati gives only study-level points: PC mean 50 Gy → 20% dysphagia, and aspiration above a PC mean of 60 Gy. |
| **Larynx mean <50 Gy (aspiration <30%)** | Only in Brodin Table 1. Not found in Rancati's recommendation text. |
| **Brainstem** | Present (Mayo). |
| **Larynx** | Present (Rancati). |
| **Cochlea** | Present (Bhandare). |
| **Optic nerve / chiasm** | **Absent from the QUANTEC primary.** The Mayo optic paper is abstract-only in this folder. Only Brodin Table 1 (secondary) gives max <55 Gy, <3%. |
| **Brachial plexus** | **No conventional-fractionation constraint** in any of these files. Only HYHOPE's hypofractionated goals exist. |
| **Oral cavity, post-arytenoid/cricoid space, cervical esophagus, brain** (conventional source values) | **Absent.** HYHOPE gives hypofractionated goals only. Brodin gives oral mucosa V30 <73% only in its highlights. |
| **Mandible Dmax** (conventional) | **No point-dose constraint** in these files. MDACC is volume-based (V44/V58). Brodin quotes D2%/mean >60 Gy in ORN cases (De Felice), which is a description of cases, not a limit. |
| **Spinal cord + 5 mm and brainstem + 5 mm PRV** (conventional) | **Absent.** |
| **Marks 2010 QUANTEC summary table** | Not in the local text. A PDF is needed. |
