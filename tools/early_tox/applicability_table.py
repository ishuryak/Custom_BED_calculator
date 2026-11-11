"""Applicability table for the early-toxicity evidence (review F3), generated from the data.

Usage:
    python3 tools/early_tox/applicability_table.py [--step6 DIR] [--out PATH]

One row per arm: schedule features computed from data/early_tox/arm_fractions.csv (dose per
fraction, fractions per week and per day, elapsed days, longest pause), planned therapy, technique,
scale group, assessment schedule and role, plus the family's step 6 score difference (candidate
minus Fowler-parameter benchmark, conditional task) for primary families and external checks.
Writes a CSV; a schedule class (accelerated / hyperfractionated / hypofractionated / split /
conventional) is assigned by a stated rule, not by hand.
"""
import argparse
import csv
import importlib.util
import os
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_spec = importlib.util.spec_from_file_location("step5_fit", os.path.join(ROOT, "tools", "early_tox", "step5_fit.py"))
s5 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s5)
OUT_DEFAULT = os.path.join(ROOT, "docs", "development", "2026_11_10_early_toxicity", "applicability_table.csv")
STEP6_DEFAULT = os.path.join(ROOT, "docs", "development", "2026_11_10_early_toxicity", "step6_review")


def schedule_features(doses, times):
    days = sorted({int(t // 24) for t in times})
    per_day = defaultdict(int)
    for t in times:
        per_day[int(t // 24)] += 1
    elapsed = (times[-1] - times[0]) / 24
    pauses = [b - a - 1 for a, b in zip(days, days[1:])]
    return dict(n_fractions=len(doses), total_dose_gy=round(sum(doses), 2),
                max_dose_per_fraction=round(max(doses), 2), mean_dose_per_fraction=round(sum(doses) / len(doses), 3),
                max_fractions_per_day=max(per_day.values()), elapsed_days=round(elapsed, 2),
                fractions_per_week=round(len(doses) / ((elapsed + 1) / 7), 2), longest_pause_days=max(pauses or [0]))


def schedule_class(f):
    """Stated rule, in order: split course if a pause > 4 days; hypofractionated if the mean dose
    per fraction > 2.1 Gy; hyperfractionated if > 1 fraction a day with mean dose <= 1.5 Gy;
    accelerated if more than 11 Gy a week (calendar weeks, first to last fraction); else
    conventional (70 Gy / 35 over 47 days is 10.4 Gy a week)."""
    weekly = f["total_dose_gy"] / ((f["elapsed_days"] + 1) / 7)
    if f["longest_pause_days"] > 4:
        return "split course"
    if f["mean_dose_per_fraction"] > 2.1:
        return "hypofractionated"
    if f["max_fractions_per_day"] > 1 and f["mean_dose_per_fraction"] <= 1.5:
        return "hyperfractionated"
    if weekly > 11:
        return "accelerated"
    return "conventional"


def build(step6_dir):
    fr = {}
    for r in s5.read("arm_fractions.csv"):
        fr.setdefault(r["schedule_id"], []).append((int(r["fraction"]), float(r["dose_gy"]), float(r["time_h"])))
    arms = {a["arm_id"]: a for a in s5.read("arms.csv")}
    mapping = {m["arm_id"]: m for m in s5.read("endpoint_mapping.csv") if m["analysis_set"] != "secondary"}
    pt = s5.planned_therapy()
    scores = {}
    pf = os.path.join(step6_dir, "per_family.csv")
    if os.path.exists(pf):
        for r in csv.DictReader(open(pf)):
            scores[r["family"]] = ("held-out primary", r["d_conditional"])
    ex = os.path.join(step6_dir, "external.csv")
    if os.path.exists(ex):
        ext = defaultdict(dict)
        for r in csv.DictReader(open(ex)):
            ext[r["family"]][r["model"]] = float(r["conditional"])
        hyp = []
        for f, v in ext.items():
            if "candidate" not in v or "benchmark" not in v:
                continue
            d = round(v["candidate"] - v["benchmark"], 4)
            if f.startswith("HYPNO"):
                hyp.append(f"{f[7:-1]} {d:+.3f}")                  # 'HYPNO [code_0]' -> 'code_0'
            else:
                scores[f] = ("external", d)
        if hyp:
            scores["HYPNO"] = ("external, chemotherapy scenarios", "; ".join(hyp))
    rows = []
    for aid, m in mapping.items():
        a = arms[aid]
        variants = a["schedule_variants"].split()
        sid = variants[0].split(":")[0]
        v = sorted(fr[sid])
        f = schedule_features([d for _, d, _ in v], [t for _, _, t in v])
        role = {"primary": "primary (held-out validated)", "sensitivity": "external sensitivity",
                "benchmark": "external benchmark (HYPNO)"}.get(m["analysis_set"], m["analysis_set"])
        fam = m["family_id"]
        sc = scores.get(fam, ("", ""))
        is_ref = aid == s5.REF_ARM.get(fam) or aid == "HYPNO_NFX"
        rows.append(dict(family=fam, arm_id=aid, reference_arm=int(is_ref), role=role,
                         schedule_class=schedule_class(f), **f,
                         schedule_variants=len(variants), technique=a["technique"],
                         planned_therapy=(pt.get(aid, {}).get("agents") or a["planned_systemic_therapy"]),
                         therapy_status=pt.get(aid, {}).get("status", "none" if a["planned_systemic_therapy"] == "none" else "text only"),
                         scale_group=m["scale_group"], nominal_window=m["window"],
                         last_assessment_weeks_after_rt=m["last_assessment_weeks_after_rt"],
                         max_post_rt_gap_weeks=m["max_post_rt_gap_weeks"],
                         adequate_ascertainment=m["adequate_ascertainment"],
                         analysis_population=m["analysis_population"],
                         family_score_task=sc[0], family_d_conditional=sc[1]))
    rows.sort(key=lambda r: (r["role"], r["family"], -r["reference_arm"], r["arm_id"]))
    return rows


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--step6", default=STEP6_DEFAULT)
    ap.add_argument("--out", default=OUT_DEFAULT)
    a = ap.parse_args(argv[1:])
    rows = build(a.step6)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} arms written to {os.path.relpath(a.out, ROOT)}")


if __name__ == "__main__":
    main(sys.argv)
