---
title: "HNSCC BED calculator, version 4.4: user guide"
author: "Igor Shuryak, MD, PhD – is144@cumc.columbia.edu"
---

# What the calculator is for

The calculator lists shorter radiotherapy schedules for head and neck squamous cell carcinoma (HNSCC) whose modeled doses stay close to the standard 70 Gy in 35 fractions of 2 Gy over 47 days. It is a tool for designing trials and generating hypotheses. A listed schedule is a candidate for a trial, not an equivalent schedule and not a treatment recommendation. Site: https://ishuryak.shinyapps.io/custom_bed_calculator_v4/ (the page footer and the Model and sources page give the version and the deployed revision). The guide for version 3.2.1 (`HNSCC_BED_v321_user_guide.md`) is kept as a historical record and does not describe this layout.

# The main page

## What is listed

Every schedule of 5 fractions a week, started on a Monday, with 5 to 30 fractions (at least 5 fewer than the standard) and a dose per fraction above 2 Gy, checked in 0.05 Gy steps up to 8 Gy. A schedule is listed when, compared with 35 × 2 Gy in 47 days, all three conditions hold:

1. Tumor: biologically effective dose (BED) at most 5 below the standard's (the default limit), under at least one of the two tumor models.
2. Early side effects (acute mucositis): BED at most 5 above the standard's (default).
3. Late side effects: BED at most 5 above the standard's (default), at the prescription dose.

BED is in units of the tissue's α/β: Gy₁₀ for tumor and mucosa and Gy₃ for late effects at the defaults. The page rewrites its units when an α/β is changed. A 5 Gy₁₀ tumor loss is about 4.2 Gy in 2-Gy fractions (EQD2). At another tumor α/β it is 5/(1 + 2/(α/β)) Gy, for example 2.8 Gy at α/β 2.5 Gy. The three margins are chosen exploratory allowances. They are not established non-inferiority or toxicity margins. The boxes under the conditions ("Change the limits") set each margin separately, from 0 to 20 Gy of BED. Smaller values make a condition stricter and larger values looser, and "Reset all to 5" restores the defaults. The list, the plot (one dashed line per limit, in the color of its bars), the refit support, the ZIP export and the Report card follow the values set. While a box is empty or out of range, for example in the middle of retyping a number, the page keeps the last valid limits and a red note beside the boxes says so.

## The two tumor models

The dose-independent (DI) model starts accelerated regrowth on a fixed day. The dose-dependent (DD) model starts it once a set amount of cell kill is reached, at a rate that grows with the daily cell kill. In the default calibrated mode, DD shares α and the background growth rate with DI, and its onset constant is solved so that both models give the same tumor BED at an anchor regimen (by default the standard). Their agreement is therefore not an independent confirmation. In direct mode (Model settings), DD has its own entered parameters and need not agree with DI at any schedule. Either form is a variant of the published DD model, not its independent fit.

A model is used only when it is available and both the standard's and the schedule's tumor BED are positive. When the standard's BED under one model is not positive, the page says so and lists schedules from the other model alone.

## Reading a row

Each row is one number of fractions. The dose column gives the lowest and highest qualifying dose per fraction. Each BED change column gives the lowest and highest change over the qualifying doses. These ranges show the spread within the row. They are not uncertainty intervals. Early, late and DI tumor BEDs rise with dose, so their extremes sit at the row's lowest and highest dose. The DD tumor BED can fall with dose when a higher dose speeds regrowth more than it adds cell kill, so its extremes can sit anywhere in the row. With Details ticked, a worked example for the top row names the doses of each extreme and why the next dose below and above fails.

Worked example at the defaults (scenario D18, calibrated DD): with 20 fractions the qualifying doses are 2.65 to 2.80 Gy in 26 days. At 2.70 Gy the DI change is −2.93 Gy₁₀, which meets the tumor condition, and the DD change is −9.87 Gy₁₀, which does not. The dose is listed because one model is enough, and the Details column marks it "DI only". The standard's DI and DD BEDs are both 71.51 Gy₁₀, its early BED is 53.11 Gy₁₀ and its late BED is 116.67 Gy₃.

## Refit support and Details

The "Refit support" column is always shown. It is the share of the parameter samples of the Bayesian refit to randomized trials in which the dose meets the tumor condition under either model. It uses the refit's own parameters, not the Model settings values, and it is not a probability of tumor control, safety or non-inferiority.

Ticking Details adds two columns. "Tumor condition holds in" says whether both models or only one credit the row. Notes name doses with refit support below 80%, and also mark doses above 2.55 Gy (the largest fraction size in the trials the refit used, so an extrapolation), doses whose acute BED exceeds that of HYPNO's 55 Gy in 20 fractions, and the HYPNO test arm.

Details also shows the button "Trials the tumor models were refitted to". It opens a list of the 45 arms or strata of the 19 randomized trials the Bayesian refit used: each trial's reference with PubMed and DOI links, and each arm's total dose, dose per fraction, number of fractions, overall time and whether that time is planned or delivered, all as used in the fit. A fraction number marked ≈ is an average over the arm. Lettered notes give what the columns cannot show: drugs (nimorazole in DAHANCA 6&7, and previous chemotherapy allowed in RTOG 83-13, and no arm included chemotherapy with radiotherapy) and schedule details, for example where a paper reports only a range of doses or only whole-arm times. A star marks the 8 trials that were not among the trials of the 2018 and 2019 fits by Shuryak, Hall and Brenner.

HYPNO context: 395 patients were randomized to 55 Gy in 20 fractions. Weekly cisplatin was permitted and given to about 75%. Grade 3 or worse acute mucositis was observed in 190 of the 373 evaluable for early effects (50.9%), against 208 of 380 (54.7%) after 66 Gy in 33 fractions at 6 a week. HYPNO did not compare with 35 × 2 Gy. A modeled BED is not an observed toxicity rate.

## Opening a row on the Report card

Clicking a row shows every qualifying dose. "Open on the Report card" copies the chosen dose and the comparison 35 × 2 Gy in 47 days (5 a week from a Monday) to the Report card, which Results also uses. It also copies the tumor threshold and the early and late limits in use (−5 and standard + 5 by default). Only a dose of the row as currently listed is accepted. If the settings changed after the row was selected, the page asks for the row to be selected again.

Other pages keep their own criteria:

| Page | Uses the main-page rules after "Open"? | Own settings that stay |
|---|---|---|
| Report card, Results | Yes | none |
| Search | No | tumor criterion (default both models), ranking, grid, days saved, organ limits |
| Both-model screen | No | its own comparison and margin when overridden |
| Uncertainty | No | its own threshold and toxicity rule |
| Organs at risk | No | the organ limits entered |

A schedule can therefore be listed on the main page and absent from Search: that is a stricter analysis, not a contradiction.

## Saving the list

"Download this list (ZIP)" writes every qualifying dose at full precision (`points.csv`), the rows as displayed (`summary.csv`), the refit support with the SHA-256 of the draws it came from, and every setting (`scenario.json`). `README.txt` gives the R commands that recompute the same tables. The limits in use are recorded, and replaying the file uses them. A file from another version, with other draws, with changed fixed rules (standard, grid or calendar) or edited after export is refused. When the refit draws could not be read, the file records support as not available and never substitutes a value. The Candidate download on the Search page records the Search analysis, not the main list.

# What the calculator does not include

The tumor parameters come mainly from trials of radiotherapy without chemotherapy. With concurrent chemotherapy, shortening treatment may gain less, and no chemotherapy effect is fitted. HPV and smoking status, tumor volume, target volumes and the placement of treatment gaps are not modeled. Late effects are judged at the prescription dose only. The Organs at risk page converts maximum-dose limits with an organ α/β and compares mean-dose limits on physical dose. It cannot certify spinal cord, brainstem, parotid, laryngeal or swallowing safety, which needs the actual organ doses of a plan.

# Settings

The Model settings page holds the tumor, acute and late parameters, the DD mode and the DD anchor. The main list follows them. Its comparison (35 × 2 Gy), margins, calendar and grid are fixed. The default tumor values (scenario D18) were chosen by the developer and are not fitted to outcome data. A fit to randomized trials is available as a preset. The equations and sources are on the Model and sources page.
