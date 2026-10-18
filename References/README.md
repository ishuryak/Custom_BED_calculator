# References for the HNSCC BED calculator v2 (acute mucosal and late toxicity BED)

Compiled 2026-10-18. Every entry below was confirmed to exist via PubMed E-utilities
(esearch/esummary/efetch) and, for DOIs, Crossref and/or Unpaywall. Open-access status
was taken from Unpaywall (queried 2026-10-18). Downloaded files were opened and their
title page checked against the citation.

Verbatim quotes are marked with quotation marks and their location. "Abstract (PubMed)"
means the value was checked against the PubMed abstract only, not the full text.

**Git note.** Only files whose licence permits redistribution (CC BY, CC BY-NC, public protocol)
should ever be committed. Paywalled items are recorded here and were NOT downloaded.

---

## A. Acute oral/pharyngeal mucosal BED

### A1. Fowler JF, Harari PM, Leborgne F, Leborgne JH (2003) - PRIMARY SOURCE for acute parameters
- Citation: Fowler JF, Harari PM, Leborgne F, Leborgne JH. Acute radiation reactions in oral and pharyngeal mucosa: tolerable levels in altered fractionation schedules. Radiother Oncol. 2003;69(2):161-168.
- DOI: 10.1016/s0167-8140(03)00231-7 · PMID: 14643953 · https://doi.org/10.1016/s0167-8140(03)00231-7
- Status: PAYWALLED (Unpaywall: closed). Not downloaded. Full text not seen.
- Supports: acute mucosal BED form and parameters; acute tolerance zone.
- Verbatim, Abstract (PubMed), Methods: "The standard published form of BED with its linear overall time factor was used: BED=nd(1 + d/(alpha/beta) - Ln2(T - T(k))/alpha T(p) (see text for parameters)."
- Verbatim, Abstract (PubMed), Results: "practical schedules used for treatment caused stimulated compensatory proliferation to start at about 7 days ... we suggest choosing a doubling time of 2.5 days as a datum. With these parameters a 'tolerance zone of uncertainty' could be identified which predicted acute-reaction acceptability or not of a schedule within a range of about 2-10 Gy in total BED. If concurrent chemoradiotherapy is used, our provisional suggestion is that this zone should be reduced by up to roughly 3-5 Gy10 in BED"
- Verbatim, Abstract (PubMed), Conclusions: "Parameters of alpha=0.35 Gy-(1), alpha/beta=10 Gy, Tk=7 days and Tp=2.5 days are suggested. The 'acute/ tolerance zone' then turns out to be 59-61 Gy10 for radiation-only treatments. Further information about the decrement caused by concurrent head-and-neck cancer chemoradiotherapy, possibly 3-5 Gy10, is required."
- Note: the abstract's bracket placement is as printed in PubMed; the intended form is nd(1 + d/(α/β)) − ln2·(T − Tk)/(α·Tp).

### A2. Fowler JF (2010) "21 years of biologically effective dose"
- Citation: Fowler JF. 21 years of biologically effective dose. Br J Radiol. 2010;83(991):554-568.
- DOI: 10.1259/bjr/31372149 · PMID: 20603408 · PMCID: PMC3473681 · https://pmc.ncbi.nlm.nih.gov/articles/PMC3473681/
- Status: FREE TO READ in PMC but NOT downloaded: PMC, Europe PMC and the publisher (academic.oup.com) all returned bot-challenge pages to scripted download, and the article is not in the PMC open-access (reuse-licensed) subset. Download manually in a browser if a local copy is wanted.
- Full text read via the PMC HTML page (automated page reader, 2026-10-18); quotes below should be re-checked against the PDF before publication.
- Supports: acute parameters (same as A1), acute "grey zone" upper limit 63 Gy10, late BED3 constraint; cited by the HYHOPE protocol (A4) as its ref 23.
- Verbatim, Table 4 item 3: "Acute BED10 = nd(1+d/10)−0.693 (T−7d)/0.35×2.5d"
- Verbatim, Table 6 caption: "Acute mucosal constraints are 51 Gy EQD with α/β = 10 Gy, α = 0.35 Gy−1, Tk = 7d, Tp = 2.5d."
- Verbatim, text: "proposed a "Grey Zone" of 59–63 Gy10 (corresponding to EQD10/2 = 49–52.5 Gy10/2)"; "the middle of the acute grey zone: (61 Gy10 = EQD of 51 Gy10/2)"
- Verbatim, Table 4 item 1 (late): "late BED3≤117 Gy3≤EQD3/2 70 Gy unless volume reduced" and "α/β = 3 Gy, or 2 Gy for central nervous system and kidney"

### A3. Hartley A, et al. (2011) - the Warwick WRAP record https://wrap.warwick.ac.uk/41625
- Citation: Hartley A, Sanghera P, Kazi W, Mehanna H, McConkey C, Glaholm J, Fowler J. Correlation of currently used radiobiological parameters with local control and acute and late mucosal toxicity in randomised studies of altered fractionation for locally advanced head and neck cancer. Clin Oncol (R Coll Radiol). 2011;23(1):29-33.
- DOI: 10.1016/j.clon.2010.08.007 · PMID: 20829006
- Status: PAYWALLED. The WRAP record (eprint 41625) is METADATA ONLY ("Restricted or Subscription Access"; full text "not available from this repository"). Not downloaded.
- Supports: acute mucosal BED with α = 0.3; late mucosal α/β = 3; empirical check against randomized trials.
- Verbatim, Abstract (PubMed), Methods: "The tumour biologically effective dose (tBED; α/β=10Gy, t(k) [onset time of accelerated repopulation]=22 days, t(p) [average doubling time during accelerated repopulation]=3 days, α=0.3Gy(-1)), acute mucosal biologically effective dose (amBED; α/β=10Gy, t(k)=7 days, t(p)=2.5 days, α=0.3Gy(-1)) and late mucosal biologically effective dose (lmBED; α/β=3Gy) were calculated"
- Verbatim, Abstract (PubMed), Results: "A trend towards a correlation was seen between the percentage difference in amBED and the percentage difference in acute grade 3 mucositis (P=0.06). A significant correlation was observed between the percentage difference in lmBED and the percentage difference in grade 3 late mucosal toxicity (P=0.02). However, a 15% decrease in lmBED between control and experimental arms of the study was necessary for any sparing of late mucosal toxicity to be observed."

### A4. HYHOPE protocol (UTSW STU 2020-1079; NCT04580446)
- Citation: Moon D (PI). HYHOPE: Phase I study of de-intensified hypofractionated radiation therapy for human papilloma virus-associated oropharynx cancer. Protocol STU 2020-1079, version 4.1 (document footer "Mod_17, 04-02-24"). UT Southwestern Simmons Comprehensive Cancer Center. ClinicalTrials.gov NCT04580446.
- URL: https://cdn.clinicaltrials.gov/large-docs/46/NCT04580446/Prot_SAP_000.pdf
- Status: DOWNLOADED -> `Moon_2024_HYHOPE_protocol_NCT04580446.pdf` (47 pp). Publicly posted on ClinicalTrials.gov. Note: page 1 carries a "CONFIDENTIAL ... property of the UTSW Simmons Comprehensive Cancer Center" template header despite public posting; do not commit to a public repository without checking.
- Supports: an applied use of the Fowler acute parameters and of the "BED3 ≤ standard, BED10 ≥ standard" design rule (the same rule as calculator v2).
- Verbatim, section 4.1.2 "Dose and fractionation", p. 10 (PDF p. 15): "Tk = time delay of repopulation (21 days for tumor, 7 days for mucosa) / Tp = cell doubling time (3 days for tumor, 2.5 days for mucosa) / α/β = 3 for late complications, 10 for tumor/mucosa / α = 0.35 Gy-1"; "Conventionally fractionated RT: 40 days (30 fractions) / Hypofractionated RT: 19 days (12-15 fractions), 26 days (20 fractions)"
- Verbatim, same section: "Dose and fractionation were selected to satisfy two parameters (see Table): ≤ BED3 vs conventional fractionation (60 Gy), i.e. late toxicity dose to normal tissue; ≥ BED10 vs conventional fractionation (60 Gy), i.e. dose to tumor"
- Table (same page), BED tumor / acute mucosa / late (Gy): 44.4 Gy/12 fx: 61 / 52 / 99; 46.5 Gy/15 fx: 61 / 51 / 95; 52 Gy/20 fx: 62 / 50 / 97; comparison 60 Gy/30 fx: 59 / 46 / 100.
- Verbatim, same section: "the dose ranges are well below the limits deemed tolerable (<59-63 Gy10) based on analysis of toxicities associated with various head and neck cancer fractionation regimens published to date (23)." Ref 23 = Fowler 2010 (A2).
- Spot check (this compilation): comparison arm acute BED = 60×1.2 − 0.792×(40−7) = 45.9, matching "46".
- Section 4.1.2 also records that dose level 1 (12 fx) was eliminated after "an apparent increase in acute grade 2 and grade 3 toxicities related to oral/pharyngeal mucositis and pain" in a parallel trial, despite its modeled acute BED of 52 Gy10.

### A5. Fenwick JD, et al. (2008) - alternative acute tolerance boundary
- Citation: Fenwick JD, Lawrence GP, Malik Z, Nahum AE, Mayles WP. Early mucosal reactions during and after head-and-neck radiotherapy: dependence of treatment tolerance on radiation dose and schedule duration. Int J Radiat Oncol Biol Phys. 2008;71(2):625-634.
- DOI: 10.1016/j.ijrobp.2008.02.008 · PMID: 18474317
- Status: PAYWALLED. Not downloaded.
- Supports: a non-linear, time-dependent tolerance boundary on cell-kill BED (no Tk/Tp time factor).
- Verbatim, Abstract (PubMed), Conclusions: "All schedules analyzed (regular, gap, and accelerating) with BED(CK) values below BED(CKboundary)(T)=69.5x(T/32.2)/sin((T/32.2)((radians)))-3.5Gy(10)(forT< or =50 days) are tolerable, and many lying above the boundary are intolerable."
- Verbatim, Methods: "Total cell-kill biologically effective doses (BED(CK)) have been calculated for 84 schedules, including incomplete repair effects, but making no other corrections for the effect of schedule duration T."

### A6. Pardo-Montero J, Parga-Pazos M, Fenwick JD (2021) - re-fit of Fowler, Fenwick and Strigari models on 108 schedules
- Citation: Pardo-Montero J, Parga-Pazos M, Fenwick JD. Classification of tolerable/intolerable mucosal toxicity of head-and-neck radiotherapy schedules with a biomathematical model of cell dynamics. Med Phys. 2021;48(7):4075-4084.
- DOI: 10.1002/mp.14834 · PMID: 33704792 · PMCID: PMC8362027
- Status: DOWNLOADED -> `PardoMontero_2021_mucosal_tolerance_model.pdf` (from the PMC open-access S3 bucket). Licence CC BY-NC.
- Also DOWNLOADED -> `PardoMontero_2021_SupportingInfo_108_schedules.xls` (supporting information: per-schedule dose per fraction, fraction times in hours, and Tol/Intol classification). Useful as a validation set for the calculator.
- Verbatim, Appendix A.1 (Eq. A1): "BEDboundary(T) = BED0 + log(2)(T − Tk)/(αTp) ... where Tk is the accelerated proliferation kick-off time, and Tp is the doubling time of stem cells in the mucosa." Appendix A.2: "T is the treatment time (in days), including the first day of treatment".
- Verbatim, Methods p. 4: "We have fixed the values of α/β = 10 Gy (as typically done for mucositis12), and the incomplete repair parameter, T1/2 = 3.2 h ... The BED model of Fowler has three free-fitting parameters, αTd, Tk, and BED0"
- Table I (p. 5), refitted "Fowler Ref. [6]" column (read from the rendered table): BED0 = 67.38 [Gy10], Tk = 19.44 d, αTd = 0.677 d Gy^-1, α/β = 10 Gy (fixed), T1/2 = 3.2 h (fixed); sensitivity 1.000, specificity 0.785. Fenwick column: BED0 = 64.25, Tc = 45.60 d. Strigari column: TD50(0) = 48.67, Δ = 84.44 d.
- Verbatim, Abstract: "The predictive value of our model is very similar to that of the phenomenological models of Fowler (1.785), Fenwick (1.806), and Strigari (1.774). When performing a k = 2 cross-validation ... For Fowler, the worsening is higher, down to 1.49."

### A7. Bentzen SM, Saunders MI, Dische S, Bond SJ (2001) - CHART mucosal Dprolif
- Citation: Bentzen SM, Saunders MI, Dische S, Bond SJ. Radiotherapy-related early morbidity in head and neck cancer: quantitative clinical radiobiology as deduced from the CHART trial. Radiother Oncol. 2001;60(2):123-135.
- DOI: 10.1016/s0167-8140(01)00358-9 · PMID: 11439207
- Status: PAYWALLED. Not downloaded.
- Supports: an independent clinical estimate of the mucosal time factor.
- Verbatim, Abstract (PubMed): "The dose recovered per 1-day protraction of overall treatment time, D(prolif), was estimated at 0.80 Gy with 95% confidence limits 0.7 and 1.1 Gy/day for human mucosa." Also: "Confluent mucositis developed earlier after the start of treatment (2.9 vs. 4.9 weeks)" and "a significant dose-volume or dose-area effect exists".

### A8. Denham JW, et al. (1995) - mucosal α/β
- Citation: Denham JW, Hamilton CS, Simpson SA, O'Brien MY, Ostwald PM, Kron T, Dear KB. Acute reaction parameters for human oropharyngeal mucosa. Radiother Oncol. 1995;35(2):129-137.
- DOI: 10.1016/0167-8140(95)01545-r · PMID: 7569021
- Status: PAYWALLED (not checked in Unpaywall separately; Elsevier 1995, assumed closed). Not downloaded.
- Verbatim, Abstract (PubMed): "Analysis of duration by direct and indirect methods suggest alpha/beta ratios in the range 7-10 Gy and half-time (t1/2) values in the range 0.27-0.5 h".
- Joiner 2023 slides (B5) list Denham et al. (1995) mucositis α/β = 9.3 (95% CL 5.8; 17.9).

### A9. Fowler JF (2007) and (2008), optimum overall time papers
- Fowler JF. Is there an optimum overall time for head and neck radiotherapy? A review, with new modelling. Clin Oncol (R Coll Radiol). 2007;19(1):8-22. DOI 10.1016/j.clon.2006.09.008 · PMID 17305251. PAYWALLED. Erratum: Clin Oncol 2008;20(2):124-6 (listed in PubMed).
- Fowler JF. Optimum overall times II: Extended modelling for head and neck radiotherapy. Clin Oncol (R Coll Radiol). 2008;20(2):113-126. DOI 10.1016/j.clon.2007.11.003 · PMID 18155893. PAYWALLED.
- Verbatim, 2008 Abstract (PubMed): "Standard linear quadratic modelling was used, assuming at first alpha/beta=10 Gy, alpha=0.35 ln/Gy, Tk=21 days, Tp=3 days for tumours, but with Tk=7 days, Tp=2.5 Gy for acute mucosal reactions, as before. A late complications constraint of 70 Gy was accepted, and an acute constraint of 51 Gy (both at 2 Gy fractions)." ("Tp=2.5 Gy" is a typo for days in the abstract.)

### A10. Other acute-tolerance models (cited, not downloaded, all PAYWALLED)
- Strigari L, et al. A new model for predicting acute mucosal toxicity in head-and-neck cancer patients undergoing radiotherapy with altered schedules. Int J Radiat Oncol Biol Phys. 2012;83(5):e697-702. DOI 10.1016/j.ijrobp.2012.02.004 · PMID 22578541. Abstract: AUC 0.84 (95% CI 0.75-0.92) for tolerable vs intolerable.
- Meade S, Sanghera P, Glaholm J, Hartley A. Models of acute mucosal tolerance to radiotherapy alone applied to synchronous chemoradiation schedules in head and neck cancer. Tumour Biol. 2014;35(3):2017-23. DOI 10.1007/s13277-013-1267-0 · PMID 24104502. Verbatim abstract: "a global value for the contribution of chemotherapy of 5.1 Gy10".
- Fenwick JD, Pardo-Montero J, Nahum AE, Malik ZI. Impact of schedule duration on head and neck radiotherapy: accelerated tumor repopulation versus compensatory mucosal proliferation. Int J Radiat Oncol Biol Phys. 2012;82(2):1021-30. DOI 10.1016/j.ijrobp.2010.12.048 · PMID 21377280.
- Fowler JF. Sensitivity analysis of parameters in linear-quadratic radiobiologic modeling. Int J Radiat Oncol Biol Phys. 2009;73(5):1532-7. DOI 10.1016/j.ijrobp.2008.11.039 · PMID 19306749.

---

## B. Late normal-tissue α/β and HN tumor α/β

### B1. Stuschke M, Thames HD (1999)
- Citation: Stuschke M, Thames HD. Fractionation sensitivities and dose-control relations of head and neck carcinomas: analysis of the randomized hyperfractionation trials. Radiother Oncol. 1999;51(2):113-121.
- DOI: 10.1016/s0167-8140(99)00042-0 · PMID: 10435801
- Note: published in Radiotherapy and Oncology, NOT IJROBP. The 1997 IJROBP paper by the same authors (Stuschke M, Thames HD. Hyperfractionated radiotherapy of human tumors: overview of the randomized clinical trials. Int J Radiat Oncol Biol Phys. 1997;37(2):259-67; DOI 10.1016/s0360-3016(96)00511-1; PMID 9069295) is a different paper.
- Status: PAYWALLED. Not downloaded.
- Verbatim, Abstract (PubMed): "The tumours showed small fractionation sensitivities characterized by an overall alpha/beta of 10.5 (6.5-29) Gy. One trial allowed quantitative estimation of the alpha/beta values for late normal tissue damage: The alpha/beta estimate for late effects of grade 2 + was 4.0 (3.3-5.0) Gy, assuming a fixed gamma50 of 5 and was even smaller for smaller gamma50 values." Also: "An overall gamma50 of 3.1 (1.5-4.7)".

### B2. van Leeuwen CM, et al. (2018)
- Citation: van Leeuwen CM, Oei AL, Crezee J, Bel A, Franken NAP, Stalpers LJA, Kok HP. The alfa and beta of tumours: a review of parameters of the linear-quadratic model, derived from clinical radiotherapy studies. Radiat Oncol. 2018;13:96.
- DOI: 10.1186/s13014-018-1040-z · PMID: 29769103 · PMCID: PMC5956964
- Status: DOWNLOADED -> `vanLeeuwen_2018_alpha_beta_tumours.pdf` and `vanLeeuwen_2018_alpha_beta_tumours_AdditionalFile1.pdf`. Licence CC BY 4.0.
- Supports: HN tumor α/β (tumours only; no normal-tissue values).
- Verbatim, Results p. 3: "Estimates of α/β for head & neck, cervix, bladder and liver tumours generally indicated low fractionation sensitivity (mostly, α/β = 10 Gy or α/β = −10 Gy)"; "I2 estimates of α/β for tumours of the head & neck, prostate and skin were 87, 94 and 97%, respectively."
- Fig. 1 (p. 4, read from the rendered figure): 39 HN estimates ranging from −83.6 to 30.0 Gy; no pooled HN value is given. Entry "[252] Stuschke 1999 (US) 10.6 [6.5,28.6]".

### B3. Rezvani M, et al. (1991) - BIR 3 vs 5 fractions/week trial
- Citation: Rezvani M, Alcock CJ, Fowler JF, Haybittle JL, Hopewell JW, Wiernik G. Normal tissue reactions in the British Institute of Radiology Study of 3 fractions per week versus 5 fractions per week in the treatment of carcinoma of the laryngo-pharynx by radiotherapy. Br J Radiol. 1991;64(768):1122-1133.
- DOI: 10.1259/0007-1285-64-768-1122 · PMID: 1773272
- Status: PAYWALLED. Not downloaded.
- Abstract (PubMed) states only: "For a number of end-points, alpha/beta ratios and N- and T-exponents of a modified nominal standard dose (NSD) formula have been calculated." No values in the abstract.
- Values available only from the secondary source B5: HN "Various late effects 3.5 (1.1; 5.9)" and oral mucosa "Mucositis 15 (–15; 45)", both attributed to "Rezvani et al. (1991)". NOT verified against the primary paper.

### B4. Maciejewski B, Taylor JM, Withers HR (1986) - supraglottic larynx late α/β
- Citation: Maciejewski B, Taylor JM, Withers HR. Alpha/beta value and the importance of size of dose per fraction for late complications in the supraglottic larynx. Radiother Oncol. 1986;7(4):323-6.
- DOI: 10.1016/s0167-8140(86)80061-5 · PMID: 3809591 · PAYWALLED.
- Verbatim, Abstract (PubMed): "Using two different methods, an alpha/beta value of 3.8 Gy was calculated. A strong correlation between the risk of late effect in the larynx and size of dose per fraction was found; there was no influence of overall treatment time."

### B5. University of Toronto file https://radonc.utoronto.ca/media/2690/download (identified)
- Identity: lecture slides, "Ch.9 The Linear-Quadratic approach to fractionation", Michael Joiner, Toronto 2023, Basic Clinical Radiobiology course (28 slides). Slide 16 reproduces "Table 10.1: α/β values for human normal tissues and tumors" (from Joiner & van der Kogel, Basic Clinical Radiobiology, 5th ed.; slide 11 cites "Table 9.1, Basic Clinical Radiobiology 5th Ed" for the animal table).
- Status: DOWNLOADED -> `Joiner_2023_LQ_fractionation_UofT_lecture.pdf`. Publicly downloadable teaching material; licence not stated. SECONDARY source.
- A summary box ("Mean Late 2.9 / Mean Early 10.6 / H&N, Lung tumors high, Breast, Prostate tumors low") overlays the middle of the table on screen; the values below were read from the PDF text layer underneath.
- Table 10.1 rows relevant to HN (α/β Gy; 95% CL; source), verbatim:
  - Oral mucosa, Mucositis: 9.3; 5.8; 17.9; Denham et al. (1995)
  - Oral mucosa, Mucositis: 15; –15; 45; Rezvani et al. (1991)
  - Oral mucosa, Mucositis: ~8; N/A; Chogule and Supe (1993)
  - Subcutis, Fibrosis: 1.7; 0.6; 2.6; Bentzen and Overgaard (1991)
  - Skin/vasculature, Telangiectasia: 2.8; 1.7; 3.8; Turesson and Thames (1989)
  - Spinal cord, Myelopathy: <3.3; N/A; Dische et al. (1981)
  - Nerve, Brachial plexopathy: <3.5; N/A; Olsen et al. (1990) and ~2; N/A; Powell et al. (1990)
  - Head and neck, Various late effects: 3.5; 1.1; 5.9; Rezvani et al. (1991)
  - Head and neck, Various late effects: 4.0; 3.3; 5.0; Stuschke and Thames (1999)
  - Supraglottic larynx, Various late effects: 3.8; 0.8; 14; Maciejewski et al. (1986)
  - Oral cavity + oropharynx, Various late effects: 0.8; –0.6; 2.5; Maciejewski et al. (1990)
  - Tumours, Head and neck, Various: 10.5; 6.5; 29; Stuschke and Thames (1999)
  - Larynx 14.5 (4.9; 24) Rezvani 1993; Vocal cord ~13 ('Wide') Robertson 1993; Buccal mucosa 6.6 (2.9; ∞) and Tonsil 7.2 (3.6; ∞) Maciejewski 1989; Nasopharynx 16 (–11; 43) Lee 1995.

### B6. Kirkpatrick JP, van der Kogel AJ, Schultheiss TE (2010) - QUANTEC spinal cord
- Citation: Kirkpatrick JP, van der Kogel AJ, Schultheiss TE. Radiation dose-volume effects in the spinal cord. Int J Radiat Oncol Biol Phys. 2010;76(3 Suppl):S42-49.
- DOI: 10.1016/j.ijrobp.2009.04.095 · PMID: 20171517 · PAYWALLED (Unpaywall: closed).
- Verbatim, Abstract (PubMed): "the estimated risk of myelopathy is <1% and <10% at 54 Gy and 61 Gy, respectively, with a calculated strong dependence on dose/fraction (alpha/beta = 0.87 Gy.)"

### B7. Bentzen SM, Christensen JJ, Overgaard J, Overgaard M (1988)
- Citation: Some methodological problems in estimating radiobiological parameters from clinical data. Alpha/beta ratios and electron RBE for cutaneous reactions in patients treated with postmastectomy radiotherapy. Acta Oncol. 1988;27(2):105-16. DOI 10.3109/02841868809090330 · PMID 3390341.
- Status: not downloaded (abstract gives no α/β value; breast, not HN). Listed because "Bentzen & Overgaard" was a requested candidate; the fibrosis α/β 1.7 in B5 is attributed to "Bentzen and Overgaard (1991)", which was NOT located/verified here.

---

## C. Hypofractionated HNSCC comparison regimens

### C1. HYPNO trial - full peer-reviewed publication EXISTS
- Citation: Bentzen SM, Gupta T, Jacinto AA, Rosenblatt E, Bhasker S, Napoles M, Binia S, Chansilpa Y, Aguiar S, Nagarajan M, Laskar SG, De Marchi P, Goloubeva O, Hopkins K, Zubizarreta E, Agarwal JP, Abdel-Wahab M. Hypofractionated Versus Normofractionated Accelerated Radiation Therapy With or Without Cisplatin for Locally Advanced Head and Neck Squamous Cell Carcinoma (HYPNO): A Randomized, Open-Label, Phase 3, Noninferiority Trial. Int J Radiat Oncol Biol Phys. 2026;125(2):387-398. Epub 2026 Feb 19.
- DOI: 10.1016/j.ijrobp.2025.12.056 · PMID: 41711624 · NCT02765503 (as given by the requester; not separately checked here)
- Status: OPEN ACCESS (hybrid, CC BY-NC-ND 4.0 per Unpaywall and the Elsevier API) but NOT downloaded: redjournal.org and sciencedirect.com returned bot-challenge / HTTP 403 pages to scripted requests, and the Elsevier API returned metadata only. Download manually from https://www.redjournal.org/article/S0360-3016(26)00020-9/pdf.
- Verbatim, Abstract (PubMed): "a 20-fraction hypofractionated (HFX) schedule delivering 55 Gy in 20 fractions, 5 fractions per week, over 4 weeks would be noninferior to a 33-fraction accelerated, normofractionated (NFX) 2 Gy per fraction schedule, delivering 66 Gy in 2-Gy fractions, 6 fractions per week over 5.5 weeks"; "792 patients were centrally randomized: 395 to HFX and 397 to NFX"; "The HYPNO test arm passed the separate noninferiority tests for both loco-regional tumor control (P = .04) and grade 3+ late adverse events (P = .004). At 3 years, the absolute difference in outcome between the 2 arms was ≤1.4 percentage points for overall survival, progression-free survival, loco-regional control, and grade 3+ late adverse events." "Concurrent chemoradiation therapy with 35 mg/m2 cisplatin weekly was permitted."
- Acute mucositis rates by arm: not in the abstract; full text not read.

### C2. Gupta T, et al. (2020) - HYPNO rationale / 55 Gy in 20 fractions review
- Citation: Gupta T, Ghosh-Laskar S, Agarwal JP. Resource-sparing curative-intent hypofractionated-accelerated radiotherapy in head and neck cancer: More relevant than ever before in the COVID era. Oral Oncol. 2020;111:105045.
- DOI: 10.1016/j.oraloncology.2020.105045 · PMID: 33091846 · PMCID: PMC7572312
- Status: DOWNLOADED -> `Gupta_2020_hypofractionated_HYPNO_perspective.pdf` (PMC copy). Licence: Elsevier COVID-19 resource centre free-access permission (cover page), NOT a CC licence; do not commit.
- Verbatim, p. 2: "Fractionation sensitivity of HNSCC is quantified by the parameter α/β (in Gy) and the estimate derived from this systematic overview is somewhat lesser than the often quoted 'textbook' value of α/β = 10 Gy [16]." Ref 16 = Stuschke & Thames 1999. Later: "If the actual α/β for HNSCC is 6.4 Gy, i.e. the best estimate from the systematic overview and meta-analysis".
- Verbatim: "assuming 6-hour inter-fraction interval and half-life for repair (T1/2) of 4.4 h. Acceleration is assumed to kick in before 26 days after the first fraction of RT."
- Table 1 lists UK/Christie/Birmingham 55 Gy/20 fx series (Sanghera 2007, Jegannathen 2010/2011, Tobias 2010, Chan 2011, Benghiat 2014, Jacinto 2018).

### C3. Jacinto AA, et al. (2018) - 55 Gy/20 fx with cisplatin (open access)
- Citation: Jacinto AA, Batalha Filho ES, Viana LS, De Marchi P, et al. Feasibility of concomitant cisplatin with hypofractionated radiotherapy for locally advanced head and neck squamous cell carcinoma. BMC Cancer. 2018;18:1026.
- DOI: 10.1186/s12885-018-4893-5 · PMID: 30352576 · PMCID: PMC6199702
- Status: DOWNLOADED -> `Jacinto_2018_hypofractionated_55Gy20fx_cisplatin.pdf`. Licence CC BY 4.0.
- Verbatim, Abstract: "The median overall treatment time was 29 days (27–34)." "The rates of grade 3 dermatitis and mucositis were 30% and 40%, respectively, with spontaneous resolution." (n = 20)
- Verbatim, Methods: "Assuming the α/β ratio to be 10Gy for head and neck carcinoma, the estimated BED for HYP-RT is about 66Gy10, which is almost equivalent to CFRT [15]."

### C4. UK / Christie 55 Gy in 20 fractions series - status
- No open-access full paper from the Christie series was found. Identified:
  - Sanghera P, et al. Hypofractionated accelerated radiotherapy with concurrent chemotherapy for locally advanced squamous cell carcinoma of the head and neck. Int J Radiat Oncol Biol Phys. 2007;67(5):1342-51. DOI 10.1016/j.ijrobp.2006.11.015 · PMID 17241752. PAYWALLED (Birmingham, not Christie).
  - Jegannathen A, et al. 2010 (Clin Oncol 22:185-91) and 2011 (Clin Oncol 23:149-58), Christie series, cited in Gupta 2020 Table 1; PubMed IDs not retrieved here.
  - Miller D, et al. A Multicentre Retrospective Cohort Study of Hypofractionated Radiotherapy in Head and Neck Cancer. Clin Oncol (R Coll Radiol). 2026;59:104332. DOI 10.1016/j.clon.2026.104332 · PMID 42759386. Unpaywall: closed. Abstract not read.

---

## D. Tumor repopulation

### D1. Withers HR, Taylor JM, Maciejewski B (1988)
- Citation: Withers HR, Taylor JMG, Maciejewski B. The hazard of accelerated tumor clonogen repopulation during radiotherapy. Acta Oncol. 1988;27(2):131-146.
- DOI: 10.3109/02841868809090333 · PMID: 3390344
- Status: DOWNLOADED -> `Withers_1988_accelerated_repopulation.pdf` (16 pp) from Medical Journals Sweden (https://medicaljournalssweden.se/actaoncologica/article/view/33350), whose page links a CC BY 4.0 licence.
- Verbatim, Abstract p. 131: "clonogen repopulation in squamous cell carcinomas of the head and neck accelerates only after a lag period of the order of 4±1 weeks after initiation of radiotherapy and that a dose increment of about 0.6 Gy per day is required to compensate for this repopulation. Such a dose increment is consistent with a 4-day clonogen doubling rate, compared with a median of about 60 days in published reports of unperturbed tumor growth rates." (The PDF text layer renders "±" as "+".)
- Verbatim, p. 140: "If it is assumed that it requires about 2.4 Gy in a 2-Gy/fraction regimen to compensate for each doubling in clonogen number, then a daily increment of 0.6 Gy in TCD50 is consistent with a clonogen doubling rate of 2.4/0.6=4 days."
- Verbatim, p. 140: "the onset of mucosal repopulation at about 2 weeks and the onset of tumor clonogen repopulation at about 3 to 5 weeks."

---

## Files in this folder

| File | Licence | Source |
|---|---|---|
| Withers_1988_accelerated_repopulation.pdf | CC BY 4.0 (MJS page) | medicaljournalssweden.se |
| vanLeeuwen_2018_alpha_beta_tumours.pdf | CC BY 4.0 | PMC OA S3 (PMC5956964) |
| vanLeeuwen_2018_alpha_beta_tumours_AdditionalFile1.pdf | CC BY 4.0 | PMC OA S3 |
| PardoMontero_2021_mucosal_tolerance_model.pdf | CC BY-NC | PMC OA S3 (PMC8362027) |
| PardoMontero_2021_SupportingInfo_108_schedules.xls | CC BY-NC | PMC OA S3 |
| Jacinto_2018_hypofractionated_55Gy20fx_cisplatin.pdf | CC BY 4.0 | PMC OA S3 (PMC6199702) |
| Gupta_2020_hypofractionated_HYPNO_perspective.pdf | Elsevier COVID free access (not CC) | PMC OA S3 (PMC7572312) |
| Moon_2024_HYHOPE_protocol_NCT04580446.pdf | public ClinicalTrials.gov document; "CONFIDENTIAL" template header | cdn.clinicaltrials.gov |
| Joiner_2023_LQ_fractionation_UofT_lecture.pdf | not stated (public teaching file) | radonc.utoronto.ca |

Not downloaded, paywalled: Fowler 2003; Hartley 2011; Fenwick 2008, 2012; Bentzen 2001; Denham 1995; Fowler 2007, 2008, 2009; Strigari 2012; Meade 2014; Stuschke & Thames 1999; Rezvani 1991; Maciejewski 1986; Kirkpatrick 2010; Sanghera 2007; Miller 2026.
Not downloaded although free to read (bot-blocked): Fowler 2010 BJR (PMC3473681); HYPNO 2026 IJROBP (CC BY-NC-ND).

### E1. Shuryak I, Hall EJ, Brenner DJ (2018) - DI and DD repopulation models
- Citation: Shuryak I, Hall EJ, Brenner DJ. Dose dependence of accelerated repopulation in head and neck cancer: supporting evidence and clinical implications. Radiother Oncol. 2018;127(1):20-26.
- DOI: 10.1016/j.radonc.2018.02.015
- Status: Elsevier, not redistributable. Local copy outside this folder: `Office/Publications/Articles/accelerated_repopulation.pdf` (the author's own paper). Not copied here.
- Supports: the DI model (Eq. 1) and the DD model (Eqs. 2-5) used by the calculator; the background growth rate g = 0.0154 per day, which the calculator omits.
- Checked 2026-10-18 against the local copy: the DD onset C = (αD + βD²/n)·TkDD/T and the rate γDD = γ·Pkill with Pkill = 1 − exp(−(αd + βd²)n/T) match `R/bed_models.R`.
