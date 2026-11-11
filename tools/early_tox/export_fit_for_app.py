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


# The one scientific profile this release may export (review R4): the step 5 primary candidate.
PRIMARY_PROFILE = dict(settings={"use_od": 1}, Tk=TK, halftime_h=HALFTIME_H, n_families=22, n_arms=51,
                       totals={"arms": 51, "families": 22, "evaluable": 6042, "events": 2945})
PACKAGE_SCHEMA = "early-fit-package-1"
EVALUATOR_ID = "early-bed-fit-1: E = D + Q/ab - r*max(0, T_elapsed - Tk), exact hinge, repair half-time h on fraction times"


def _manifest_problem(m):
    spec = importlib.util.spec_from_file_location("step5_fit", os.path.join(ROOT, "tools", "early_tox", "step5_fit.py"))
    s5 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(s5)
    return s5.manifest_problem(m)


def read_unit(unit):
    """The primary unit's manifest, diagnostics and draws, after every check (review R4)."""
    man = json.load(open(os.path.join(unit, "manifest.json")))
    problem = _manifest_problem(man)
    if problem:
        raise RuntimeError(f"{unit}: malformed manifest ({problem}); refusing to export")
    if man["gate_pass"] is not True:
        raise RuntimeError(f"{unit}: gate did not pass; refusing to export")
    if man["smoke"] is not False:
        raise RuntimeError(f"{unit}: smoke run; refusing to export")
    for fn, digest in man["outputs"].items():
        path = os.path.join(unit, fn)
        if not os.path.exists(path) or sha256(path) != digest:
            raise RuntimeError(f"{path}: changed since the run; refusing to export")
    diag = json.load(open(os.path.join(unit, "diagnostics.json")))
    if diag.get("identity") != man["identity"]:
        raise RuntimeError(f"{unit}: diagnostics belong to another run; refusing to export")
    if diag.get("gate_pass") is not True or diag.get("divergences") != 0 or not diag.get("rhat", 9) < 1.01 \
            or not diag.get("ess_bulk", 0) > 400 or not diag.get("ess_tail", 0) > 400 or not diag.get("ebfmi_min", 0) >= 0.3:
        raise RuntimeError(f"{unit}: diagnostics do not meet the gate; refusing to export")
    for k, v in PRIMARY_PROFILE.items():
        if diag.get(k) != v:
            raise RuntimeError(f"{unit}: not the approved primary profile ({k} = {diag.get(k)!r}, expected {v!r}); "
                               "refusing to export it under the primary description")
    with gzip.open(os.path.join(unit, "draws.csv.gz"), "rt") as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 4000:
        raise RuntimeError(f"{unit}: {len(rows)} draws, expected 4000")
    for k, _ in COLUMNS:
        v = np.array([float(r[k]) for r in rows])
        if not np.all(np.isfinite(v)) or (k in ("ab", "r", "kappa") and np.any(v <= 0)):
            raise RuntimeError(f"{unit}: posterior column {k} is not finite and positive")
    return man, rows, diag


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


ROLE_LABEL = {"primary": "fitted", "sensitivity": "sensitivity check only", "benchmark": "external benchmark (never fitted)"}


def assessment_text(last, gap="", note=""):
    if last == "":
        head = "follow-up after RT not stated"
    elif last == "resolution":
        head = "followed until reactions settled"
    elif float(last) == 0:
        head = "scored during RT only"
    else:
        head = f"last scored {float(last):g} week{'' if float(last) == 1 else 's'} after RT"
    if gap not in ("", None) and last not in ("", "0"):
        g = float(gap)
        head += f"; largest gap after RT {g:g} week{'' if g == 1 else 's'}"
    return f"{head} ({note})" if note else head


def trial_table():
    """One row per arm of every trial in the early-toxicity evidence (fitted, sensitivity, benchmark),
    with its verified citation (data/early_tox/source/report_citations.csv, PubMed esummary). Every
    modelled schedule variant and weight, the schedule notes, the assessment schedule and population,
    and the structured planned therapy are kept (review R2); a pooled count is never split by variant."""
    spec = importlib.util.spec_from_file_location("app_table", os.path.join(ROOT, "tools", "early_tox", "applicability_table.py"))
    app = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(app)
    feats = {r["arm_id"]: r for r in app.build("/nonexistent")}
    src = os.path.join(ROOT, "data", "early_tox", "source")
    cit = {r["family_id"]: r for r in csv.DictReader(open(os.path.join(src, "report_citations.csv")))}
    asc = {r["family_id"]: r for r in csv.DictReader(open(os.path.join(src, "ascertainment.csv")))}
    pt = {r["arm_id"]: r for r in csv.DictReader(open(os.path.join(src, "planned_therapy.csv")))}
    arms = {r["arm_id"]: r for r in csv.DictReader(open(os.path.join(ROOT, "data", "early_tox", "arms.csv")))}
    outs = {r["outcome_id"]: r for r in csv.DictReader(open(os.path.join(ROOT, "data", "early_tox", "outcomes.csv")))}
    fr = {}
    for r in csv.DictReader(open(os.path.join(ROOT, "data", "early_tox", "arm_fractions.csv"))):
        fr.setdefault(r["schedule_id"], []).append((int(r["fraction"]), float(r["dose_gy"]), float(r["time_h"])))
    rows = []
    for m in csv.DictReader(open(os.path.join(ROOT, "data", "early_tox", "endpoint_mapping.csv"))):
        if m["analysis_set"] not in ROLE_LABEL:
            continue
        f, a, o, c = feats[m["arm_id"]], arms[m["arm_id"]], outs[m["outcome_id"]], cit[m["family_id"]]
        variants = []
        for v in a["schedule_variants"].split():
            sid, w = v.split(":")
            x = sorted(fr[sid])
            vf = app.schedule_features([d for _, d, _ in x], [t for _, _, t in x])
            variants.append((float(w), vf))
        sched = "; ".join(
            (f"weight {w:g}: " if len(variants) > 1 else "") +
            f"{vf['n_fractions']} fractions, {vf['total_dose_gy']:g} Gy, {int(vf['elapsed_days']) + 1} days"
            + (f", {vf['max_fractions_per_day']} a day" if vf["max_fractions_per_day"] > 1 else "")
            for w, vf in variants)
        t = pt.get(m["arm_id"])
        therapy = a["planned_systemic_therapy"]
        if t:
            courses = f", {t['planned_courses']} courses planned" if t["planned_courses"] else ""
            share = (f", intended for about {float(t['proportion_intended']):.0%} of patients"
                     if t["proportion_intended"] not in ("", "1", "0") else "")
            status = {"known": "", "assumed_equal": " (assumed equal to the control arm: planned courses not stated)"}[t["status"]]
            therapy = f"{t['agents']}{courses}{share}{status}"
        ac = asc.get(m["family_id"], {})
        rows.append(dict(
            role=ROLE_LABEL[m["analysis_set"]], trial=a["trial_id"],
            reference=f"{c['first_author']} et al., {c['journal']} {c['year']}", pmid=c["pmid"], doi=c["doi"],
            arm=a["arm_label"], control_arm=int(f["reference_arm"]),
            schedule=sched, schedule_variants=len(variants), schedule_notes=a["pattern_notes"],
            dose_per_fraction_gy=f["mean_dose_per_fraction"], fractions_per_day=f["max_fractions_per_day"],
            schedule_class=f["schedule_class"], concurrent_therapy=therapy,
            therapy_note=(t or {}).get("basis", ""),
            endpoint=o["endpoint_as_written"], endpoint_timepoint=ac.get("endpoint_timepoint", ""),
            grading_group=m["scale_group"], events=o["events"], patients=o["denominator"],
            analysis_population=ac.get("analysis_population", ""),
            assessment=assessment_text(m["last_assessment_weeks_after_rt"], m["max_post_rt_gap_weeks"],
                                       ac.get("schedule_note", ""))))
    order = {"fitted": 0, "sensitivity check only": 1, "external benchmark (never fitted)": 2}
    rows.sort(key=lambda r: (order[r["role"]], r["trial"], -r["control_arm"], r["arm"]))
    return rows


def publish_package(stage, root=ROOT):
    """Check the staged package against its own manifest, then move it in with the manifest last, so
    an interrupted export never leaves draws and manifest from different generations (review R4)."""
    man = json.load(open(os.path.join(stage, "early_fit_manifest.json")))
    for fn, meta in man["files"].items():
        if sha256(os.path.join(stage, fn)) != meta["sha256"]:
            raise RuntimeError(f"staged {fn} does not match its manifest; nothing published")
    dest = {"early_fit_draws.csv": os.path.join(root, "draws"), "early_fit_trials.csv": os.path.join(root, "draws"),
            "early_fit_evaluator.json": os.path.join(root, "tests", "testthat", "fixtures")}
    old = os.path.join(root, "draws", "early_fit_manifest.json")
    if os.path.exists(old):
        os.remove(old)                         # the app refuses the package until the new manifest is in place
    for fn, d in dest.items():
        os.replace(os.path.join(stage, fn), os.path.join(d, fn))
    os.replace(os.path.join(stage, "early_fit_manifest.json"), old)
    os.rmdir(stage)


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--unit", default=UNIT_DEFAULT)
    a = ap.parse_args(argv[1:])
    man, rows, diag = read_unit(a.unit)
    stage = os.path.join(ROOT, "draws", f".early_fit.partial-{os.getpid()}")
    os.makedirs(stage, exist_ok=True)
    out_csv = os.path.join(stage, "early_fit_draws.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["draw"] + [c for _, c in COLUMNS])
        for i, r in enumerate(rows, 1):
            w.writerow([i] + [r[k] for k, _ in COLUMNS])
    vals = {c: np.array([float(r[k]) for r in rows]) for k, c in COLUMNS}
    med = {c: float(np.median(v)) for c, v in vals.items()}
    manifest = dict(
        schema=PACKAGE_SCHEMA, evaluator=EVALUATOR_ID,
        row_order="draws in the order of the source unit's draws.csv.gz (4 chains x 1000), one joint row per draw",
        model_id=f"early-fit-{man['identity'][:12]}",
        label="Trial-count fit (research)",
        created_by="tools/early_tox/export_fit_for_app.py",
        source_unit=os.path.relpath(a.unit, ROOT), source_identity=man["identity"],
        source_revision=diag.get("revision"), chain_archive=man.get("chain_archive"),
        endpoint="severe acute mucositis (maximum grade >= 3, or confluent), arm-level counts pooled across grading "
                 "scales, observation windows and evaluable, as-treated or completer populations; "
                 f"{diag['n_families']} randomized trial families ({diag['totals']['arms']} arms, "
                 f"{diag['totals']['evaluable']} evaluable patients); 2 of 22 families met the adequate post-RT "
                 "ascertainment rule",
        model="logit p = eta_trial + kappa * (E - E_ref) + gamma * chemo + arm effect; "
              "E = D + Q/ab - r * max(0, T_elapsed - Tk), Q with incomplete repair",
        fixed=dict(Tk_e_elapsed_days=TK, repair_halftime_h=HALFTIME_H, hinge="exact", clock="elapsed days"),
        point_values=dict(ab_e=round(med["ab_e"], 4), r_e=round(med["r_e"], 4), rule="posterior medians"),
        n_draws=len(rows), gate=dict(passed=True, rhat=diag["rhat"], ess_bulk=diag["ess_bulk"], ess_tail=diag["ess_tail"],
                                     divergences=diag["divergences"], ebfmi_min=diag["ebfmi_min"]),
        margin_meaning="+5 is a numerical margin on this model's BED scale; its clinical risk meaning is not established",
        chemotherapy_term="one pooled coarse indicator (greater or added planned chemotherapy exposure); not a "
                          "treatment-specific causal effect; not used by the calculator",
        status="research scenario; Fowler 2003 is the app default. Validation: better held-out log score across "
               "22 trial families (+0.28 per family); HYPNO depends on its chemotherapy assumption; ELAN-RT (frail, "
               "split course) favored the Fowler-parameter benchmark. Record: "
               "docs/development/2026_11_11_early_toxicity_review_response.md",
        files={"early_fit_draws.csv": dict(sha256=sha256(out_csv), n_draws=len(rows))})
    trials = trial_table()
    out_trials = os.path.join(stage, "early_fit_trials.csv")
    with open(out_trials, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(trials[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(trials)
    manifest["files"]["early_fit_trials.csv"] = dict(sha256=sha256(out_trials), n_rows=len(trials),
                                                     citations="PubMed esummary, verified 2026-11-11")
    with open(os.path.join(stage, "early_fit_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    # cross-language fixture
    cases = []
    for name, m, d, days, calendar in fixture_cases():
        doses, times = [d] * m, days_to_hours(days)
        for ab, r in ((med["ab_e"], med["r_e"]), (10.0, 0.79217), (7.5, 2.4)):
            cases.append(dict(name=name, m=m, d=d, days=days, T=days[-1], ab=ab, r=r, calendar=calendar,
                              E=ev.early_bed(doses, times, ab=ab, r=r, Tk=TK, halftime_h=HALFTIME_H)))
    with open(os.path.join(stage, "early_fit_evaluator.json"), "w") as f:
        json.dump(dict(Tk=TK, halftime_h=HALFTIME_H, cases=cases), f, indent=1)
    publish_package(stage)
    print(f"{len(rows)} draws; medians ab {med['ab_e']:.3f}, r {med['r_e']:.4f}; model {manifest['model_id']}")


if __name__ == "__main__":
    main(sys.argv)
