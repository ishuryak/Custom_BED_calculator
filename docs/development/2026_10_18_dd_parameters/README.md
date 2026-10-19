# DD model parameter analysis, 2026-10-18

Record of the analyses behind decisions D13, D14 and D18 in `../2026_10_18_fix_plan.md`. Nothing here is part of the app or its deployment (`docs/` is in `.rscignore`).

## Reports

- `DI_DD_clinical_regimens_2026_10_18.md` / `.docx`: DI and DD repopulation onset, doubling time and daily BED loss for 13 regimens (standard, UK 65/30, HYPNO arms, HYHOPE, eight randomized-trial arms), at λ_DD matched to DI (0.666) and at the current default (0.5), with predicted versus observed trial differences.
- `DD_parameter_grid_search_2026_10_18.md` / `.docx`: grid search over tumor α and λ_DD with the current calibration. Result: no combination makes onset, doubling time and daily loss plausible for all 13 regimens; with rate and loss given precedence (Igor Shuryak's decision), a narrow band α 0.265–0.405, λ_DD 0.35–0.41 works.

## Data and figures (`dd_grid_2026_10_18/`)

- `grid_variantA.csv.gz`: all (α, λ_DD) combinations, number of regimens plausible on all three quantities (loose and strict limits).
- `grid_rate_first.csv.gz`: the same grid judged on doubling time and daily loss only, with onset counts kept for ranking.
- `feasible_region_by_alpha.csv`, `candidate_points_regimens.csv`, `best_combo*.csv`: summaries quoted in the reports.
- `fig1`–`fig3`: maps of the number of plausible regimens over (α, λ_DD).

All values are model outputs for published average schedules (no patient-level data).

## Scripts (`scripts/`)

Working scripts as run on 2026-10-18 with Windows R 4.5.1 from `C:\Users\ishur\AppData\Local\Temp`. They source `R/bed_models.R` by absolute path and read `dd_table.csv` (the regimen list, written by `dd_table.R`). To re-run, adjust the paths. `read.csv` reads the `.csv.gz` files directly.

- `dd_table.R`, `dd_table_05.R`: regimen tables at λ_DD matched and 0.5.
- `dd_grid2.R`, `dd_grid3.R`, `dd_grid4.R`, `dd_figs.R`, `dd_cand.R`: grid searches, candidate points, figures.
- `dd_analysis.R`, `dd_floor.R`, `dd_presets.R`, `dd_hypno.R`, `dd_03.R`, `dd_pt.R`, `dd_pt2.R`, `dd_sugg.R`, `dd_fowler.R`, `dd_final.R`: the single-point checks quoted in the plan (λ_DD 0.3, the 7-day onset floor, presets, HYPNO, the final default proposal).
