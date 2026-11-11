"""Export the fitted early-toxicity model (step 5 primary, review generation) for the app's research mode.

Usage:
    python3 tools/early_tox/export_fit_for_app.py [--unit DIR]

Reads the primary unit (default docs/development/2026_11_10_early_toxicity/step5_review/overdispersed/U02),
checks its manifest (gate passed, draws hash unchanged) and writes:
    draws/early_fit_draws.csv      4,000 joint posterior draws: ab_e (Gy), r_e (BED per day), kappa, gamma
    draws/early_fit_manifest.json  identity, source, sha256 of the CSV, medians, fixed settings, endpoint
    tests/testthat/fixtures/early_fit_evaluator.json
                                   schedules evaluated with the Python evaluator used for fitting, so the
                                   R evaluator can be checked against it (cross-language fixture)
The app's research mode uses the medians as point values and the joint draws under uncertainty.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_spec = importlib.util.spec_from_file_location("evaluator", os.path.join(ROOT, "tools", "early_tox", "evaluator.py"))
ev = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ev)

UNIT_DEFAULT = os.path.join(ROOT, "docs", "development", "2026_11_10_early_toxicity", "step5_review", "overdispersed", "U02")
TK, HALFTIME_H = 7.0, 3.2
COLUMNS = (("ab", "ab_e"), ("r", "r_e"), ("kappa", "kappa"), ("gamma[1]", "gamma"))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def read_unit(unit):
    man = json.load(open(os.path.join(unit, "manifest.json")))
    if not man.get("gate_pass"):
        raise RuntimeError(f"{unit}: gate did not pass; refusing to export")
    if man.get("smoke"):
        raise RuntimeError(f"{unit}: smoke run; refusing to export")
    path = os.path.join(unit, "draws.csv.gz")
    if sha256(path) != man["outputs"]["draws.csv.gz"]:
        raise RuntimeError(f"{path}: changed since the run; refusing to export")
    with gzip.open(path, "rt") as f:
        rows = list(csv.DictReader(f))
    return man, rows


def days_to_hours(days):
    return [ev.HOURS_PER_DAY * (d - days[0]) for d in days]


def fixture_cases():
    """Once-daily schedules on app calendars (calendar day numbers, 1 = first treatment day)."""
    def cal(m, per_week=5, start=1, skip=()):
        out, day = [], 0
        while len(out) < m:
            day += 1
            wd = (start - 1 + day - 1) % 7 + 1
            if wd <= per_week and day not in skip:
                out.append(day)
        return out
    spec = [("35x2 Mon 5/wk", 35, 2.0, 5, 1, ()), ("20x2.75 Mon 5/wk", 20, 2.75, 5, 1, ()),
            ("30x2.2 6/wk", 30, 2.2, 6, 1, ()), ("33x2 6/wk Wed", 33, 2.0, 6, 3, ()),
            ("25x2.4 skip 12,13", 25, 2.4, 5, 1, (12, 13)), ("5x5 Mon", 5, 5.0, 5, 1, ())]
    return [(n_, m, d, cal(m, pw, st, sk), dict(fx_per_week=pw, start_weekday=st, skip_days=list(sk)))
            for n_, m, d, pw, st, sk in spec]


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--unit", default=UNIT_DEFAULT)
    a = ap.parse_args(argv[1:])
    man, rows = read_unit(a.unit)
    out_csv = os.path.join(ROOT, "draws", "early_fit_draws.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["draw"] + [c for _, c in COLUMNS])
        for i, r in enumerate(rows, 1):
            w.writerow([i] + [r[k] for k, _ in COLUMNS])
    vals = {c: np.array([float(r[k]) for r in rows]) for k, c in COLUMNS}
    med = {c: float(np.median(v)) for c, v in vals.items()}
    diag = json.load(open(os.path.join(a.unit, "diagnostics.json")))
    manifest = dict(
        model_id=f"early-fit-{man['identity'][:12]}",
        label="Trial-count fit (research)",
        created_by="tools/early_tox/export_fit_for_app.py",
        source_unit=os.path.relpath(a.unit, ROOT), source_identity=man["identity"],
        source_revision=diag.get("revision"), chain_archive=man.get("chain_archive"),
        endpoint="maximum grade >= 3 acute mucositis (or confluent), arm-level counts, 22 randomized trial families "
                 f"({diag['totals']['arms']} arms, {diag['totals']['evaluable']} evaluable patients)",
        model="logit p = eta_trial + kappa * (E - E_ref) + gamma * chemo + arm effect; "
              "E = D + Q/ab - r * max(0, T_elapsed - Tk), Q with incomplete repair",
        fixed=dict(Tk_e_elapsed_days=TK, repair_halftime_h=HALFTIME_H, hinge="exact", clock="elapsed days"),
        point_values=dict(ab_e=round(med["ab_e"], 4), r_e=round(med["r_e"], 4), rule="posterior medians"),
        n_draws=len(rows), gate=dict(passed=True, rhat=diag["rhat"], ess_bulk=diag["ess_bulk"], ess_tail=diag["ess_tail"],
                                     divergences=diag["divergences"], ebfmi_min=diag["ebfmi_min"]),
        status="research scenario; Fowler 2003 is the app default. Validation: better held-out log score across "
               "22 trial families (+0.28 per family); HYPNO depends on its chemotherapy assumption; ELAN-RT (frail, "
               "split course) favored the Fowler-parameter benchmark. Record: "
               "docs/development/2026_11_11_early_toxicity_review_response.md",
        files={"early_fit_draws.csv": dict(sha256=sha256(out_csv), n_draws=len(rows))})
    with open(os.path.join(ROOT, "draws", "early_fit_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    # cross-language fixture
    cases = []
    for name, m, d, days, calendar in fixture_cases():
        doses, times = [d] * m, days_to_hours(days)
        for ab, r in ((med["ab_e"], med["r_e"]), (10.0, 0.79217), (7.5, 2.4)):
            cases.append(dict(name=name, m=m, d=d, days=days, T=days[-1], ab=ab, r=r, calendar=calendar,
                              E=ev.early_bed(doses, times, ab=ab, r=r, Tk=TK, halftime_h=HALFTIME_H)))
    fx = os.path.join(ROOT, "tests", "testthat", "fixtures", "early_fit_evaluator.json")
    with open(fx, "w") as f:
        json.dump(dict(Tk=TK, halftime_h=HALFTIME_H, cases=cases), f, indent=1)
    print(f"{len(rows)} draws; medians ab {med['ab_e']:.3f}, r {med['r_e']:.4f}; model {manifest['model_id']}")


if __name__ == "__main__":
    main(sys.argv)
