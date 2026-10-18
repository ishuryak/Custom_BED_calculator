# Fixes

## 2026-10-18: version 2 (branch `v2-toxicity-hypofractionation`)

**Changes published numbers (version 1 website):** yes, for any tumor α/β other than 10 Gy. At α/β = 10 Gy (the default) nothing moves.

| Item | Record |
|---|---|
| What was wrong | Version 1 calibrates the DD model so that the DD and DI tumor BEDs are equal for 35 × 2 Gy in 47 days. The calibration hard-codes that regimen's value for α/β = 10 Gy (the constants 84 and 3948). At any other α/β, the DD tumor BED was wrong at every regimen. |
| Proof it was real | With the v1 formula at α/β = 5 Gy (other defaults unchanged), 35 × 2 Gy gives BED_DD = 73.80 Gy, while BED_DI = 79.00 Gy. The model defines them as equal there. At α/β = 15 Gy: 62.14 vs 60.33. |
| Test that pins it | `tests/testthat/test_bed_models.R::DD is calibrated to DI at the reference for every tumor alpha/beta (v1 is not)`. Confirmed failing when the hard-coded value 84 is restored (deliberate-break check, 2026-10-18). |
| Verification | `Rscript -e "testthat::test_dir('tests/testthat')"` gave 1131 expectations, 0 failed, on 2026-10-18 at this branch's first commit. Baseline before version 2: no tests existed. |
| What it moves | BED_DD for every regimen when tumor α/β ≠ 10 Gy. Example: 35 × 2 Gy at α/β = 5 Gy goes from 73.80 to 79.00 Gy. Version 1 values at α/β = 10 Gy are reproduced to a relative difference below 1e-8. |

Other changes in version 2 (new features, review fixes, and decisions) are recorded in `docs/development/2026_10_18_code_review.md`.
