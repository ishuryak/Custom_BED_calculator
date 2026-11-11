# Extraction brief: acute mucositis by randomized arm (HNSCC fractionation trials)

## Why
We are building a provenance-checked data table for a Bayesian model of acute oral/pharyngeal
mucositis risk as a function of radiotherapy schedule (dose per fraction, number of fractions,
overall time). The model compares arms WITHIN each randomized trial. Each number you extract may
end up in a published calculator, so a wrong count is much worse than a missing one. A second
person is extracting the same trials independently, and the two extractions will be reconciled.
Do not look for or read any other extraction file.

## Primary endpoint (decided)
Maximum grade >= 3 ACUTE mucositis per patient, per arm. Record the exact scale and version
(RTOG/EORTC acute, NCI CTC v1/v2/v3, CTCAE v3/v4/v5, WHO, Dische, or the paper's own) and the
exact threshold as written (e.g. "grade 3-4", "confluent mucositis", "grade >= 3", "severe").
Confluent mucositis is recorded as its own endpoint; it is mapped to the primary group later, not
by you. Also record, when reported, these SECONDARY endpoints in the same structure: duration of
confluent / grade 3 mucositis, feeding tube use during treatment, and unplanned treatment breaks
due to acute toxicity. Do not record late effects.

## Rules (hard)
1. Every number needs a verbatim quote (or exact table cell description) and a location
   (page number of the PDF and table/figure/section). No quote, no number.
2. Counts: record exactly what is reported. If only a percentage is given, record the percentage
   and the denominator the paper states for that analysis. Fill `events` ONLY if the count is
   reported directly OR it is the unique integer consistent with the percentage and the stated
   denominator at the paper's reporting precision; in that case set
   `count_derivation = "from_pct_unique"` and show the arithmetic. Otherwise leave `events` empty
   and set `count_derivation = "not_derivable"` with the reason.
3. Denominator: the number EVALUATED for acute toxicity in that arm if stated; otherwise the number
   randomized/treated, flagged as such (`denominator_type`). Never assume.
4. A missing value is empty, never zero. Do not fill gaps from memory, other papers, abstracts or
   reviews. If the same trial is reported in more than one file you were given, extract each
   report separately and note the overlap.
5. Read tables and figures from the PDF itself (the Read tool renders PDF pages; use `pages`
   ranges of at most 20). The `.txt` versions are machine text extraction and garble tables; use
   them only to locate passages. If a number is only readable from a figure (e.g. a curve),
   record it with `source_kind = "figure_read"` and say how it was read; prefer not to.
6. Schedule per arm as PLANNED in the protocol: total dose, dose per fraction (all components for
   concomitant boost / mixed schedules), fractions per day, interfraction interval (h), treatment
   days per week, planned overall time (days, and say whether it counts calendar days from first
   to last fraction or weeks), any planned gap/split. Also the delivered median dose/time if
   reported (separately).
7. Context per trial: definitive vs postoperative; tumor sites; stage; field size (early glottic /
   small fields?); concurrent systemic therapy PLANNED per arm (drug, dose, schedule) and delivered
   if reported; technique (2D, 3D, IMRT); years of accrual.
8. Assessment of acute toxicity: who scored, how often (weekly?), from when to when (during RT
   only? until N weeks after?). If not stated, say "not stated".
9. Machine load: other heavy jobs are running on this computer. Do not start OCR, installs,
   or any CPU-heavy process. Reading files is fine.

## Output
Write ONE JSON file (path given in your task) containing a list of objects, one per
(report, arm, endpoint). Fields:

```
trial_id            short name (e.g. "DAHANCA 6&7", "RTOG 9003")
report_file         file name you read
pmid                as printed on the paper if visible, else ""
randomized          true/false (is this arm part of a randomized comparison in this report?)
arm_label           as in the paper
context             "definitive" | "postoperative" | "mixed" | other text
sites               text
early_glottic_or_small_field   true/false/unknown
planned_systemic_therapy       text ("none" if none)
delivered_systemic_therapy     text or ""
technique           text or ""
schedule_planned    {total_dose_gy, dose_per_fraction_gy (text if several), n_fractions,
                     fractions_per_day, interfraction_interval_h, days_per_week,
                     planned_overall_time_days, overall_time_definition, gaps, quote, location}
schedule_delivered  {median_total_dose_gy, median_overall_time_days, quote, location} or {}
endpoint            "max_grade_ge3_mucositis" | "confluent_mucositis" | "duration_mucositis" |
                    "feeding_tube" | "treatment_break" | "other:<text>"
endpoint_as_written exact wording
scale               exact scale and version or "not stated"
assessment          {frequency, start, end, quote, location}
n_randomized        integer or ""
denominator         integer or ""
denominator_type    "evaluated_for_acute_toxicity" | "randomized" | "treated" | "other:<text>"
events              integer or ""
percent_reported    number or ""
count_derivation    "direct" | "from_pct_unique" | "not_derivable"
derivation_note     text (arithmetic, or why not derivable)
value_quote         verbatim
value_location      "p.N, Table M" etc.
source_kind         "text" | "table" | "figure_read"
notes               anything a reconciler must know (overlap with another report, ambiguity)
```

Also include at the end one object per report with `"summary": true` listing: what the report
contains on acute mucositis, what is missing, and any reason the trial should not be in the
primary analysis. Then reply with a short summary (<= 250 words): trials done, arms with a usable
grade>=3 or confluent count, and anything surprising.
