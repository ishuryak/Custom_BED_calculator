# Custom BED calculator, version 2

R Shiny app for head and neck squamous cell carcinoma (HNSCC). It computes three biologically effective doses (BEDs):

- tumor BED under the dose-independent (DI) and dose-dependent (DD) repopulation models
- early (acute mucosal) toxicity BED
- late toxicity BED.

It also lists shorter, hypofractionated regimens whose tumor BED gain over a reference regimen (default 35 × 2 Gy in 47 days) meets a chosen minimum under the DI or DD model, while early and late toxicity BED stay within chosen limits. Published regimens (the HYPNO trial arms) are shown for comparison.

Version 1 is unchanged: https://ishuryak.shinyapps.io/custom_bed_calculator/ (repo `ishuryak/Custom_BED_calculator`).

## Files

- `R/bed_models.R`: all formulas, as plain functions.
- `app.R`: the user interface.
- `about.md`: the "Model and sources" page.
- `tests/testthat/`: tests of the formulas against an independent implementation, and of the server logic.
- `run_app.bat`: double-click launcher for Windows.
- `.rscignore`: keeps `References/`, `tests/` and `docs/` out of a shinyapps.io deployment.
- `docs/development/`: review records.
- `References/`: sources, with `References/README.md` giving each one's license and the quoted parameter values. Paywalled or non-redistributable PDFs are git-ignored.

## Run locally

```r
shiny::runApp("path/to/Custom_BED_calculator_v2")
testthat::test_dir("tests/testthat")
```

## Change from version 1

Version 1 calibrates the DD model so that the DD and DI tumor BEDs are equal for 35 × 2 Gy in 47 days, but hard-codes that regimen's value at a tumor α/β of 10 Gy. At any other α/β the calibration is wrong. For example, at α/β = 5 Gy the DD tumor BED is 73.8 Gy instead of 79.0 Gy. Version 2 recalculates the calibration for any α/β and any reference regimen. At α/β = 10 Gy the two versions agree to a relative difference below 1e-8. Version 2 also marks the DD model as unavailable when no calibration with a repopulation onset after the start of treatment exists, instead of showing a meaningless number.

Igor Shuryak, MD, PhD, is144@cumc.columbia.edu
