"""Plan step 5.3: separate historical classifier on the Pardo-Montero tolerable/intolerable labels.

Usage:
    python3 tools/early_tox/step5_historical.py [--out DIR]

Records: complete-course records only (primary 103, 13 intolerable; sensitivity 107, 15). Partial
courses never enter (plan R1, acceptance test 2). Each model is a logistic regression of the label
on one exposure, logit p = a + b (E - 60) / 10, with weak priors a ~ N(0, 2.5), b ~ N(0, 2.5), fitted
by maximum a posteriori. Cohort structure: leave-one-source-out cross-validation over the source_ref
groups (61 in the primary set); no record of a held-out source is used in its fold, including the
choice of Tk for the time-boundary model.

Exposures:
    fowler       Fowler 2003 set, complete repair (ab 10, r 0.79217, Tk 7), as in the app
    boundary     same, with Tk chosen from 0-35 days inside each training fold (alternative time
                 boundary)
    step5        step 5 primary candidate at its posterior medians (ab, r from
                 step5/overdispersed/U02), 3.2 h repair, Tk 7: a consistency check, not a validation
    null         intercept only
    fowler_h3.2  Fowler's parameters with 3.2 h repair     } matched-repair comparisons (review F8):
    step5_complete  step 5 point with complete repair      } compare like with like
A descriptive compatibility check, not an independent validation: Fowler 2003 was tuned on largely
the same historical schedules. Not used to predict patient risk.
"""
import argparse
import csv
import importlib.util
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import minimize

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_spec = importlib.util.spec_from_file_location("evaluator", os.path.join(ROOT, "tools", "early_tox", "evaluator.py"))
ev = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ev)

OUT_DEFAULT = os.path.join(ROOT, "docs", "development", "2026_11_10_early_toxicity", "step5", "historical")
TK_GRID = list(range(0, 36))


def read(name):
    with open(os.path.join(ROOT, "data", "early_tox", name), newline="") as f:
        return list(csv.DictReader(f))


def load_records(selection="primary"):
    fr = {}
    for r in read("fractions.csv"):
        fr.setdefault(r["schedule_id"], []).append((int(r["fraction"]), float(r["dose_gy"]), float(r["time_h"])))
    fr = {k: ([d for _, d, _ in sorted(v)], [t for _, _, t in sorted(v)]) for k, v in fr.items()}
    recs = []
    for r in read("labels.csv"):
        if r[f"include_{selection}"] != "1":
            continue
        if r["course_scope"] != "complete":
            raise ValueError(f"{r['record_id']}: {r['course_scope']} course in the {selection} selection")
        doses, times = fr[r["schedule_id"]]
        recs.append(dict(record_id=r["record_id"], source=r["source_ref"], schedule_id=r["schedule_id"],
                         y=int(r[f"event_{selection}"]), doses=doses, times=times))
    return recs


STEP5_DEFAULT = os.path.join(ROOT, "docs", "development", "2026_11_10_early_toxicity", "step5")


def step5_point(step5_dir=STEP5_DEFAULT):
    import pandas as pd
    s = pd.read_csv(os.path.join(step5_dir, "overdispersed", "U02", "summary.csv"), index_col=0)
    return float(s.loc["ab", "50%"]), float(s.loc["r", "50%"])


def exposure(recs, model, Tk=7.0, step5=None):
    """Exposure per record. Matched-repair variants (review F8): fowler_h3.2 is Fowler's parameters
    with 3.2 h repair; step5_complete is the step 5 point with complete repair."""
    if model == "null":
        return np.zeros(len(recs))
    if model in ("fowler", "boundary"):
        return np.array([ev.early_bed(x["doses"], x["times"], ab=ev.FOWLER_2003["ab"], Tk=Tk) for x in recs])
    if model == "fowler_h3.2":
        return np.array([ev.early_bed(x["doses"], x["times"], ab=ev.FOWLER_2003["ab"], Tk=7.0, halftime_h=3.2)
                         for x in recs])
    ab, r = step5
    h = None if model == "step5_complete" else 3.2
    return np.array([ev.early_bed(x["doses"], x["times"], ab=ab, r=r, Tk=7.0, halftime_h=h) for x in recs])


def fit(E, y, slope=True):
    x = (E - 60) / 10

    def nlp(t):
        lin = t[0] + (t[1] * x if slope else 0)
        return -(np.sum(y * lin - np.logaddexp(0, lin)) - 0.5 * (t[0] / 2.5) ** 2 - 0.5 * (t[1] / 2.5) ** 2)
    return minimize(nlp, np.zeros(2), method="BFGS").x


def predict(t, E, slope=True):
    lin = t[0] + (t[1] * (E - 60) / 10 if slope else 0)
    return 1 / (1 + np.exp(-lin))


def log_score(p, y):
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return y * np.log(p) + (1 - y) * np.log(1 - p)


def folds(recs):
    srcs = sorted({x["source"] for x in recs})
    return [(s, np.array([x["source"] != s for x in recs]), np.array([x["source"] == s for x in recs])) for s in srcs]


def choose_tk(recs):
    """Tk on TK_GRID maximizing the fitted log likelihood of these records (training data only)."""
    y = np.array([x["y"] for x in recs])

    def ll(tk):
        E = exposure(recs, "boundary", tk)
        return np.sum(log_score(predict(fit(E, y), E), y))
    return float(max(TK_GRID, key=ll))


def loso(recs, model, step5=None):
    """Held-out log score per record, Tk for the boundary model chosen on training data only."""
    y = np.array([x["y"] for x in recs])
    out, tks = np.zeros(len(recs)), []
    for _, tr, te in folds(recs):
        Tk = 7.0
        if model == "boundary":
            Tk = choose_tk([recs[i] for i in np.flatnonzero(tr)])
            tks.append(Tk)
        E = exposure(recs, model, Tk, step5)
        t = fit(E[tr], y[tr], slope=model != "null")
        out[te] = log_score(predict(t, E[te], slope=model != "null"), y[te])
    return out, tks


def auc(score, y):
    pos, neg = score[y == 1], score[y == 0]
    return float(np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg]))


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT_DEFAULT)
    ap.add_argument("--step5-dir", default=STEP5_DEFAULT, help="step 5 output folder holding overdispersed/U02")
    a = ap.parse_args(argv[1:])
    os.makedirs(a.out, exist_ok=True)
    s5pt = step5_point(a.step5_dir)
    res = dict(step5_dir=os.path.relpath(a.step5_dir, ROOT), step5_point=dict(ab=s5pt[0], r=s5pt[1]), tk_grid=[TK_GRID[0], TK_GRID[-1]], selections={})
    for sel in ("primary", "sensitivity"):
        recs = load_records(sel)
        y = np.array([x["y"] for x in recs])
        r = dict(records=len(recs), events=int(y.sum()), sources=len(folds(recs)), models={})
        ls = {}
        for m in ("null", "fowler", "boundary", "step5", "fowler_h3.2", "step5_complete"):
            ls[m], tks = loso(recs, m, s5pt)
            Tk_full = 7.0
            if m == "boundary":
                Tk_full = choose_tk(recs)
            E = exposure(recs, m, Tk_full, s5pt)
            t = fit(E, y, slope=m != "null")
            r["models"][m] = dict(
                loso_log_score=round(float(ls[m].sum()), 3),
                full_data_a=round(float(t[0]), 3), full_data_b_per_10=round(float(t[1]), 3) if m != "null" else None,
                exposure_at_p50=round(60 - 10 * t[0] / t[1], 2) if m != "null" and t[1] > 0 else None,
                in_sample_auc=round(auc(E, y), 3) if m != "null" else None,
                Tk_full_data=Tk_full if m == "boundary" else None,
                Tk_folds_range=[min(tks), max(tks)] if tks else None)
        for m in ("fowler", "boundary", "step5", "fowler_h3.2", "step5_complete"):
            d = ls[m] - ls["fowler"] if m != "fowler" else ls[m] - ls["null"]
            per_src = {}
            for i, x in enumerate(recs):
                per_src[x["source"]] = per_src.get(x["source"], 0) + d[i]
            ref = "fowler" if m != "fowler" else "null"
            r["models"][m][f"difference_vs_{ref}"] = round(float(d.sum()), 3)
            r["models"][m][f"sources_better_vs_{ref}"] = int(sum(v > 1e-9 for v in per_src.values()))
            r["models"][m][f"sources_worse_vs_{ref}"] = int(sum(v < -1e-9 for v in per_src.values()))
        r["matched_repair"] = {
            "3.2 h repair (step5 vs fowler_h3.2)": round(float((ls["step5"] - ls["fowler_h3.2"]).sum()), 3),
            "complete repair (step5_complete vs fowler)": round(float((ls["step5_complete"] - ls["fowler"]).sum()), 3)}
        res["selections"][sel] = r
        with open(os.path.join(a.out, f"records_{sel}.csv"), "w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(["record_id", "source", "schedule_id", "y", "E_fowler", "E_step5"] + [f"loso_{m}" for m in ls])
            Ef, Es = exposure(recs, "fowler"), exposure(recs, "step5", step5=s5pt)
            for i, x in enumerate(recs):
                w.writerow([x["record_id"], x["source"], x["schedule_id"], x["y"], round(Ef[i], 4), round(Es[i], 4)]
                           + [round(float(ls[m][i]), 5) for m in ls])
    with open(os.path.join(a.out, "summary.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv)
