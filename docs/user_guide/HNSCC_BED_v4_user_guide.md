---
title: "HNSCC BED calculator, version 4.6: user guide"
author: "Igor Shuryak, MD, PhD – is144@cumc.columbia.edu"
---

# What the calculator is for

The calculator lists shorter radiotherapy schedules for head and neck squamous cell carcinoma (HNSCC) whose modeled doses stay close to the standard 70 Gy in 35 fractions of 2 Gy over 47 days. It is a tool for designing trials and generating hypotheses. A listed schedule is a candidate for a trial, not an equivalent schedule and not a treatment recommendation. Site: https://ishuryak.shinyapps.io/custom_bed_calculator_v4/ (the page footer and the Model and sources page give the version and the deployed revision). The main page is "Suggested regimens". The menu opens Check a regimen, Custom search (advanced), Uncertainty, Organs at risk, and Settings and sources (with the tabs Model settings and Model and sources). The guide for version 3.2.1 (`HNSCC_BED_v321_user_guide.md`) is kept as a historical record and does not describe this layout.

# The main page

## What is listed

Every schedule of 5 fractions a week, started on a Monday, with 5 to 30 fractions (at least 5 fewer than the standard) and a dose per fraction above 2 Gy, checked in 0.05 Gy steps up to 8 Gy. A schedule is listed when, compared with 35 × 2 Gy in 47 days, all three conditions hold:

1. Tumor: biologically effective dose (BED) at most 5 below the standard's (the default limit), under at least one of the two tumor models.
2. Early side effects (acute mucositis): BED at most 5 above the standard's (default), under the early-toxicity model selected on the page (Fowler 2003 by default, see below).
3. Late side effects: BED at most 5 above the standard's (default), at the prescription dose.

BED is in units of the tissue's α/β: Gy₁₀ for tumor and mucosa and Gy₃ for late effects at the defaults. The page rewrites its units when an α/β is changed. A 5 Gy₁₀ tumor loss is about 4.2 Gy in 2-Gy fractions (EQD2). At another tumor α/β it is 5/(1 + 2/(α/β)) Gy, for example 2.8 Gy at α/β 2.5 Gy. The three margins are chosen exploratory allowances. They are not established non-inferiority or toxicity margins. The boxes under the conditions ("Change the limits") set each margin separately, from 0 to 20 Gy of BED. Smaller values make a condition stricter and larger values looser, and "Reset all to 5" restores the defaults. The list, the plot (one dashed line per limit, in the color of its bars), the refit support, the ZIP export and Check a regimen follow the values set. While a box is empty or out of range, for example in the middle of retyping a number, the page keeps the last valid limits and a red note beside the boxes says so.

## The two tumor models

The dose-independent (DI) model starts accelerated regrowth on a fixed day. The dose-dependent (DD) model starts it once a set amount of cell kill is reached, at a rate that grows with the daily cell kill. In the default calibrated mode, DD shares α and the background growth rate with DI, and its onset constant is solved so that both models give the same tumor BED at an anchor regimen (by default the standard). Their agreement is therefore not an independent confirmation. In direct mode (Model settings), DD has its own entered parameters and need not agree with DI at any schedule. Either form is a variant of the published DD model, not its independent fit.

A model is used only when it is available and both the standard's and the schedule's tumor BED are positive. When the standard's BED under one model is not positive, the page says so and lists schedules from the other model alone.

## Reading a row

Each row is one number of fractions. The dose column gives the lowest and highest qualifying dose per fraction. Each BED change column gives the lowest and highest change over the qualifying doses. These ranges show the spread within the row. They are not uncertainty intervals. Early, late and DI tumor BEDs rise with dose, so their extremes sit at the row's lowest and highest dose. The DD tumor BED can fall with dose when a higher dose speeds regrowth more than it adds cell kill, so its extremes can sit anywhere in the row. With Details ticked, a worked example for the top row names the doses of each extreme and why the next dose below and above fails.

At the defaults (the trial-fit tumor values, calibrated DD and the Fowler 2003 early model) the list has 69 doses in 18 rows, from 13 to 30 fractions. Worked example: with 20 fractions the qualifying doses are 2.65 to 2.80 Gy in 26 days. At 2.65 Gy the DI change is −4.98 Gy₁₀, which meets the tumor condition, and the DD change is −5.05 Gy₁₀, which does not. The dose is listed because one model is enough, and the Details column marks it "DI only". The standard's DI and DD BEDs are both 64.80 Gy₁₀, its early BED is 53.11 Gy₁₀ and its late BED is 116.67 Gy₃.

## The early-toxicity model: Fowler 2003 or the trial-count fit

A box on the main page selects the early-toxicity model. Fowler 2003 is the default. It uses fixed mucosal values (α/β 10 Gy, doubling time 2.5 days, onset on day 7), which can be changed on the Model settings tab.

The trial-count fit is a research option. It was fitted, as a Bayesian hierarchical model, to the numbers of patients with grade 3 or worse (or confluent) acute mucositis in 51 arms of 22 randomized trials of altered fractionation. Its posterior medians are a mucosal α/β of 12.87 Gy and a time-loss rate of 1.55 Gy of BED per day from day 7, with incomplete repair between fractions (half-time 3.2 hours). It loses BED to treatment time about twice as fast as Fowler's values, so compressed courses look more mucosally intense. With it, the main list shrinks to 9 doses in 26 to 30 fractions. The links "How was the early-toxicity model fitted?" and "Trials used for the early-toxicity fit" beside the box give the method and each trial with its reference, arms, schedules, concurrent therapy, endpoint and counts.

Read research-mode results with three limits in mind:

- The +5 early limit is a numerical margin on this model's BED scale. Its clinical risk meaning is not established, and the model gives no patient risk.
- The fit predicted held-out trials better on average than Fowler's values, mostly in accelerated trials. For shorter hypofractionated courses the evidence is mixed: HYPNO depends on how its chemotherapy is treated, and ELAN-RT (frail older patients, split course) favored Fowler's values.
- The trials pool different grading scales, observation windows and evaluable populations.

While the fit is selected, every page names it. The shared banner above each page gives the model and its values, with a link back to the box. The Model settings tab shows the fitted values read-only and hides the Fowler fields, which are not used. The fit needs each regimen's actual treatment dates, because repair between fractions depends on their timing. Each regimen and comparison uses its own calendar, including skipped days. When an overall time is typed that its calendar cannot give, the early BED is reported as not available, with the reason, and no gaps are invented. Enter the missed treatment days as skipped days instead. The Fowler zone limit is refused with the fit, which is tied to Fowler's values.

## Refit support and Details

The "Refit support" column is always shown. It is the share of the parameter samples of the Bayesian refit to randomized trials in which the dose meets the tumor condition under either model. It uses the refit's own parameters, not the Model settings values, and it is not a probability of tumor control, safety or non-inferiority.

Ticking Details adds two columns. "Tumor condition holds in" says whether both models or only one credit the row. Notes name doses with refit support below 80%, and also mark doses above 2.55 Gy (the largest fraction size in the trials the refit used, so an extrapolation), doses whose acute BED exceeds that of HYPNO's 55 Gy in 20 fractions, and the HYPNO test arm.

Details also shows the button "Trials the tumor models were refitted to". It opens a list of the 45 arms or strata of the 19 randomized trials the Bayesian refit used: each trial's reference with PubMed and DOI links, and each arm's total dose, dose per fraction, number of fractions, overall time and whether that time is planned or delivered, all as used in the fit. A fraction number marked ≈ is an average over the arm. Lettered notes give what the columns cannot show: drugs (nimorazole in DAHANCA 6&7, and previous chemotherapy allowed in RTOG 83-13, and no arm included chemotherapy with radiotherapy) and schedule details, for example where a paper reports only a range of doses or only whole-arm times. A star marks the 8 trials that were not among the trials of the 2018 and 2019 fits by Shuryak, Hall and Brenner.

HYPNO context: 395 patients were randomized to 55 Gy in 20 fractions. Weekly cisplatin was permitted and given to about 75%. Grade 3 or worse acute mucositis was observed in 190 of the 373 evaluable for early effects (50.9%), against 208 of 380 (54.7%) after 66 Gy in 33 fractions at 6 a week. HYPNO did not compare with 35 × 2 Gy. A modeled BED is not an observed toxicity rate.

## Opening a row on Check a regimen

Clicking a row shows every qualifying dose. "Open on Check a regimen" copies the chosen dose and the comparison 35 × 2 Gy in 47 days (5 a week from a Monday) to Check a regimen. It also copies the tumor threshold and the early and late limits in use (−5 and standard + 5 by default). Only a dose of the row as currently listed is accepted. If the settings changed after the row was selected, the page asks for the row to be selected again.

Other pages keep their own criteria:

| Page | Uses the main-page rules after "Open"? | Own settings that stay |
|---|---|---|
| Check a regimen | Yes | none |
| Custom search (advanced) | No | tumor criterion (default both models), ranking, grid, days saved, organ limits |
| Both-model screen (in Custom search) | No | its own comparison and margin when overridden |
| Uncertainty | No | its own threshold and toxicity rule |
| Organs at risk | No | the organ limits entered |

A schedule can therefore be listed on the main page and absent from Custom search: that is a stricter analysis, not a contradiction. The early-toxicity model selected on the main page applies on every page.

## Saving the list

"Download this list (ZIP)" writes every qualifying dose at full precision (`points.csv`), the rows as displayed (`summary.csv`), the refit support with the SHA-256 of the draws it came from, and every setting (`scenario.json`). `README.txt` gives the R commands that recompute the same tables. The limits in use are recorded, and replaying the file uses them. A file from another version, with other draws, with changed fixed rules (standard, grid or calendar) or edited after export is refused. When the refit draws could not be read, the file records support as not available and never substitutes a value. The Candidate download on the Custom search page records that analysis, not the main list.

A file saved with the trial-count fit also records the fitted model (its identifier, point values, the SHA-256 of its posterior draws and package, its evaluator) and each regimen's calendar. It replays only with the same installed model, and it is refused if any of these differ. Files saved with Fowler 2003 do not record this block and replay as before. On the Uncertainty page with the fit, the early values are the fit's 4,000 joint posterior draws ("Uncertain") or its posterior medians ("Fixed"), and late and organ α/β keep their own sources. The page and the saved file say which prior entries were used and which were not.

# What the calculator does not include

The tumor parameters come mainly from trials of radiotherapy without chemotherapy. The trial-count early fit includes trials with concurrent chemotherapy through one coarse adjustment and does not predict chemotherapy effects. With concurrent chemotherapy, shortening treatment may gain less, and no chemotherapy effect is fitted. HPV and smoking status, tumor volume, target volumes and the placement of treatment gaps are not modeled. Late effects are judged at the prescription dose only. The Organs at risk page converts maximum-dose limits with an organ α/β and compares mean-dose limits on physical dose. It cannot certify spinal cord, brainstem, parotid, laryngeal or swallowing safety, which needs the actual organ doses of a plan.

# Settings

The Model settings tab of Settings and sources holds the tumor, acute and late parameters, the DD mode and the DD anchor. The main list follows them. Its comparison (35 × 2 Gy), calendar and grid are fixed, and its margins are set on the main page. The default tumor values are the posterior medians of the Bayesian refit to 45 arms of 19 randomized trials ("trial-fit point scenario"). Scenario D18, chosen by the developer and not fitted to outcome data, and other published values are available as presets. The acute fields apply to Fowler 2003 only. The equations and sources are on the Model and sources tab.
