"""Build data/fit_trial_list.csv: reference and drug information for every arm of the
tumor-model fitting data, for the main page's list of fitted trials (version 4.4.0).

One row per arm of the fitting file, keyed by ref_no and arm. Doses, fraction numbers
and times are NOT copied here: the app reads them from the fitting file itself, so the
list cannot disagree with what was fitted. This file adds only what the fitting file
lacks:
- trial: a readable trial name (the names of the refit explainer's source list);
- citation and doi: from the PubMed records saved by the trial search;
- chemotherapy: the `chemo` field of the full-text extraction rows of that arm;
- other_drugs: a radiosensitizer given within the arm, from the population notes that
  fit/build_provenance.py records (read in the publications on 2026-10-26);
- arm_label: a readable arm name where the fitting file uses a bare abbreviation;
- trial_note: one note per trial for the dialog's footnotes: drugs given and what a reader
  should know about the schedule values used in the fit. Each sentence was checked on
  2026-11-06 against the trial paper (fulltext_text/<PMID>.txt) or, where the PDF is not in
  hand, against the full-text extraction quote (all_extract.json);
- shuryak2018_ref: the trial's reference number in Shuryak, Hall and Brenner 2018
  (Radiother Oncol 127:20-26, PMID 29534828, Table 1 and references [23]-[55], read in the
  full text on 2026-11-06), or empty when the trial was not in that paper. CHART was cited
  there through Dische 1997 and RTOG 9003 through its first report (Fu 2000); the refit
  uses Saunders 2010 for CHART. Shuryak, Hall and Brenner 2019 (Int J Radiat Oncol Biol Phys
  104:272-278, PMID 30776451) used the same 16 publications (its references 5, 6 and 12-25),
  so this column covers both papers.

The builder refuses an arm with no extraction row, extraction rows that disagree on
chemotherapy, or a trial with no PubMed record.

Run from the repository root (the extraction folder and records are outside the repo):
  python3 tools/build_fit_trial_list.py data/refit_trials_2026_10_25c_main.csv data/fit_trial_list.csv \
      --search-dir ../References_refit_trials/search
"""
import argparse, csv, json, os, sys

TRIAL_NAMES = {
    "17011446": "CAIR (Skladowski 2006)",
    "22836063": "CAIR vs concomitant boost (Skladowski 2013)",
    "20394851": "CHART (Saunders 2010)",
    "14511925": "DAHANCA 6&7 (Overgaard 2003)",
    "1480768": "EORTC 22791 (Horiot 1992)",
    "7569030": "EORTC 22811 (Van den Bogaert 1995)",
    "9288839": "EORTC 22851 (Horiot 1997)",
    "18343310": "Ghoshal 2008",
    "11830307": "Hliniak 2002",
    "20382075": "IAEA-ACC (Overgaard 2010)",
    "9165135": "Jackson 1997",
    "24161568": "KROG 0201 (Moon 2014)",
    "32055717": "Princess Margaret (O'Sullivan 2020)",
    "2180866": "RTOG 83-13 (Cox 1990)",
    "7790243": "RTOG 88-09 (Fu 1995)",
    "10924966": "RTOG 9003 (Fu 2000)",
    "25035199": "RTOG 9512 (Trotti 2014)",
    "11439206": "TROG 91.01 (Poulsen 2001)",
    "16169681": "Yamazaki 2006",
}

# Reference number of the trial in Shuryak, Hall and Brenner 2018 (Table 1 sources).
SHURYAK_2018_REF = {
    "17011446": "52", "22836063": "53", "20394851": "44", "14511925": "49", "1480768": "23",
    "9288839": "24", "11830307": "47", "9165135": "48", "10924966": "45", "11439206": "51",
    "16169681": "55",
}

# Radiosensitizers given within a fitted arm (fit/build_provenance.py, POPULATION; Overgaard
# 2003: "All patients, except those with glottic cancers, also received the hypoxic
# radiosensitiser nimorazole").
OTHER_DRUGS = {
    "DAHANCA67": "nimorazole (hypoxic radiosensitizer) for all patients except glottic cancer, in both arms",
}

# Readable labels for arms the fitting file names only by an abbreviation, keyed by (ref_no, arm).
ARM_LABELS = {
    ("CAIR_vs_CB", "CAIR"): "CAIR: continuous accelerated, 7 days a week",
    ("CAIR_vs_CB", "CB"): "CB: concomitant boost, 7 fractions in 5 days a week",
    ("KROG0201", "CONV"): "Conventional, 2.0 Gy (CONV)",
    ("KROG0201", "HYPO"): "Hypofractionated, 2.25 Gy (HYPO)",
    ("PMH_HF", "HF"): "Hyperfractionated, 1.45 Gy twice daily (HF)",
    ("PMH_HF", "CF"): "2.55 Gy once daily (CF)",
    ("EORTC22791", "CF (arm 1)"): "Conventional fractionation (CF, arm 1)",
    ("EORTC22791", "HF (arm 2)"): "Hyperfractionation (HF, arm 2)",
}
YAMAZAKI = {"Arm A-1": "Arm A (2.0 Gy), tumor up to 2/3 of the vocal cord",
            "Arm A-2": "Arm A (2.0 Gy), tumor over 2/3 of the vocal cord",
            "Arm B-1": "Arm B (2.25 Gy), tumor up to 2/3 of the vocal cord",
            "Arm B-2": "Arm B (2.25 Gy), tumor over 2/3 of the vocal cord"}

# Trial-level notes on drugs and on the schedule values used in the fit, keyed by PMID.
TRIAL_NOTES = {
    "2180866": "Previous chemotherapy was allowed if not within 6 weeks of randomization. No arm included "
               "chemotherapy with radiotherapy.",
    "17011446": "Patients were treated first with 2.0 Gy and later with 1.8 Gy fractions, so each arm is shown as "
                "two periods. The paper gives overall times only for whole arms (median 40 days for CAIR and 54 "
                "days for control, with a planned CAIR time of 33 to 36 days). The period times used in the fit "
                "(36 and 51 days, 40 and 54 days) are therefore not period medians, although the fitting data "
                "label them so. Total doses ranged from 66 to 76 Gy.",
    "7569030": "The conventional arm allowed 70 Gy in 35 fractions over 7 weeks or 75 Gy in 44 fractions over "
               "9 weeks, and the split is not reported. The fit uses 70 Gy in 35 fractions. The trial's third arm "
               "(three fractions a day with misonidazole) is not used.",
    "20382075": "Both arms received 66 to 70 Gy in 33 to 35 fractions. The fit uses 68 Gy in 34 fractions.",
    "9165135": "Planned overall times were 45 to 48 days (once daily) and 22 to 25 days (twice daily). The fit uses "
               "47.6 and 24.25 days, labeled as mean delivered times. These match the 2018 paper's Table 1 "
               "(47.6 and 24.3 days) but were not checked against the trial paper.",
    "24161568": "Schedules differed by T stage within each arm (T1: 66 Gy in 33 fractions or 63 Gy in 28, T2: "
                "70 Gy in 35 or 67.5 Gy in 30). The fit uses fraction numbers weighted by each arm's T-stage mix.",
    "32055717": "The overall times (median 29 days) are from the trial's earlier report (2007), not from the "
                "2020 paper cited.",
    "16169681": "Arms A (2.0 Gy) and B (2.25 Gy) each had two total doses by tumor extent, and the fit uses "
                "these strata separately.",
    "20394851": "The 2018 and 2019 papers cited this trial through an earlier report (Dische 1997).",
}


def chemo_text(s):
    s = str(s).strip()
    if s.lower().startswith("none"):
        rest = s[4:].strip(" ()")
        rest = rest.replace("prior chemo allowed earlier", "earlier chemotherapy allowed")
        return "none" + (f" ({rest})" if rest else "")
    return s


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("data"); ap.add_argument("out"); ap.add_argument("--search-dir", required=True)
    a = ap.parse_args(argv)
    ex = json.load(open(os.path.join(a.search_dir, "extraction", "all_extract.json"), encoding="utf-8"))
    pm = {}
    for f in ("pubmed_records.json", "supp_records.json"):
        for r in json.load(open(os.path.join(a.search_dir, f), encoding="utf-8")):
            pm.setdefault(str(r["pmid"]), r)
    rows = list(csv.DictReader(open(a.data, encoding="utf-8")))
    out, problems = [], []
    for r in rows:
        pmid = str(r["pmid"])
        m = [e for e in ex if str(e["pmid"]) == pmid and str(e["arm"]).strip() == r["arm"].strip()]
        if not m:
            problems.append(f"{r['ref_no']} / {r['arm']}: no extraction row of PMID {pmid}")
            continue
        chemo = {chemo_text(e.get("chemo", "")) for e in m}
        if len(chemo) != 1:
            problems.append(f"{r['ref_no']} / {r['arm']}: extraction rows disagree on chemotherapy: {sorted(chemo)}")
            continue
        if pmid not in pm or pmid not in TRIAL_NAMES:
            problems.append(f"{r['ref_no']}: no PubMed record or trial name for PMID {pmid}")
            continue
        p = pm[pmid]
        title = p["title"].rstrip(".")
        note = [f"Other drugs: {OTHER_DRUGS[r['ref_no']]}." if r["ref_no"] in OTHER_DRUGS else "",
                TRIAL_NOTES.get(pmid, "")]
        out.append(dict(ref_no=r["ref_no"], arm=r["arm"], pmid=pmid, trial=TRIAL_NAMES[pmid],
                        arm_label=ARM_LABELS.get((r["ref_no"], r["arm"]), YAMAZAKI.get(r["arm"], "") if pmid == "16169681" else ""),
                        trial_note=" ".join(x for x in note if x),
                        citation=f"{p['first_author']} et al. {title}. {p['journal']} {p['year']}.",
                        doi=p.get("doi", ""), chemotherapy=chemo.pop(),
                        other_drugs=OTHER_DRUGS.get(r["ref_no"], ""),
                        shuryak2018_ref=SHURYAK_2018_REF.get(pmid, "")))
    if problems:
        sys.exit("REFUSED:\n" + "\n".join(problems))
    with open(a.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), lineterminator="\n")
        w.writeheader(); w.writerows(out)
    print(len(out), "arms,", len({o["pmid"] for o in out}), "trials ->", a.out)


if __name__ == "__main__":
    main(sys.argv[1:])
