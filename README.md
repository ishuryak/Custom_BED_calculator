# Custom_BED_calculator

**Use the current calculator (version 4):** https://ishuryak.shinyapps.io/custom_bed_calculator_v4/ · **Its source:** the [`v4` branch](https://github.com/ishuryak/Custom_BED_calculator/tree/v4)

A web calculator for designing trials of shorter, hypofractionated radiotherapy for head and neck squamous cell carcinoma (HNSCC). It lists schedules whose modeled tumor biologically effective dose (BED) stays close to 70 Gy in 35 fractions while the modeled early and late side-effect doses rise by no more than chosen limits, under dose-independent (DI) and dose-dependent (DD) tumor repopulation models fitted to randomized trials. It is a research tool for hypothesis generation, not for treatment decisions.

## This branch is the archived version 1

The code on `main` is version 1 (2024), kept as a record. Its DD calibration assumes tumor α/β = 10 Gy and is incorrect at any other α/β (at α/β 5 Gy, 35 × 2 Gy gives a DD BED of 73.8 Gy instead of 79.0). Do not use it for new work.

The websites of versions 1, 2 and 3 were retired on 2026-10-30 and no longer respond. Their code remains in this repository: version 1 on `main`, version 2 on [`v2-toxicity-hypofractionation`](https://github.com/ishuryak/Custom_BED_calculator/tree/v2-toxicity-hypofractionation), version 3 on [`v3`](https://github.com/ishuryak/Custom_BED_calculator/tree/v3).

Igor Shuryak, MD, PhD, is144@cumc.columbia.edu
