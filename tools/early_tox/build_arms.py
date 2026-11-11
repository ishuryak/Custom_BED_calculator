"""Build the randomized-arm early-toxicity tables (plan step 2).

Usage:
    python3 tools/early_tox/build_arms.py [REPO_ROOT]

Inputs (committed):
    data/early_tox/source/arm_curation.csv       one row per arm: role, context, therapy, schedule pattern
    data/early_tox/source/outcome_curation.csv   one row per arm outcome: expected denominator and events
    docs/development/2026_11_10_early_toxicity/extraction/out_{A1,A2,B1,B2}.json
                                                 two independent extractions per batch of reports
Outputs (data/early_tox/):
    arms.csv, outcomes.csv, arm_schedules.csv, arm_fractions.csv, evidence_families.csv

An outcome is written only if BOTH independent extractions of its report contain the same
denominator and event count. Quotes, locations, scales and assessment windows are copied from the
extraction, never typed here.

Schedule pattern grammar (one variant; variants separated by " | ", optional weight "{w}"):
    W5 | W6 | W7       treatment weekdays: Mon-Fri, Mon-Sat, every day (day 1 is a Monday)
    I<h>               hours between fractions on the same day (default 6)
    <N>d[d1;d2;...]    N treatment days, each with these fraction doses in Gy
    gap<N>             N calendar days without treatment after the last treatment day
    resume@<D>         next treatment day is the first allowed weekday on or after calendar day D
Fraction times are hours from the first fraction.
"""
import csv
import hashlib
import json
import os
import re
import sys
from collections import defaultdict

WEEKDAYS = {"W5": {0, 1, 2, 3, 4}, "W6": {0, 1, 2, 3, 4, 5}, "W7": set(range(7))}
EXTRACTIONS = {"A": ("out_A1.json", "out_A2.json"), "B": ("out_B1.json", "out_B2.json")}


def parse_variant(text):
    """Return (fractions [(dose_gy, time_h)], weight) for one pattern variant."""
    m = re.search(r"\{([0-9.]+)\}\s*$", text)
    weight = float(m.group(1)) if m else 1.0
    text = text[: m.start()] if m else text
    allowed, interval, day, last_day, out = None, 6.0, 0, None, []
    for tok in text.split():
        if tok in WEEKDAYS:
            allowed = WEEKDAYS[tok]
        elif re.fullmatch(r"I[0-9.]+", tok):
            interval = float(tok[1:])
        elif re.fullmatch(r"gap[0-9]+", tok):
            if last_day is None:
                raise ValueError(f"{tok}: gap before any treatment")
            day = last_day + 1 + int(tok[3:])
        elif re.fullmatch(r"resume@[0-9]+", tok):
            day = max(day, int(tok[7:]) - 1)
        elif re.fullmatch(r"[0-9]+d\[[0-9.;]+\]", tok):
            if allowed is None:
                raise ValueError("pattern must start with W5, W6 or W7")
            n, doses = tok.split("d[")
            doses = [float(x) for x in doses.rstrip("]").split(";")]
            for _ in range(int(n)):
                while day % 7 not in allowed:
                    day += 1
                out += [(d, day * 24 + k * interval) for k, d in enumerate(doses)]
                last_day, day = day, day + 1
        else:
            raise ValueError(f"unknown schedule token {tok!r}")
    if not out:
        raise ValueError(f"empty schedule {text!r}")
    t0 = out[0][1]
    return [(d, round(t - t0, 3)) for d, t in out], weight


def schedule_id(fr):
    key = ";".join(f"{d:.3f}@{t:.2f}" for d, t in fr)
    return "S" + hashlib.sha256(key.encode()).hexdigest()[:10]


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def as_int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def load_extractions(xdir):
    ex = {}
    for batch, files in EXTRACTIONS.items():
        ex[batch] = [[r for r in json.load(open(os.path.join(xdir, f))) if not r.get("summary")] for f in files]
    return ex


def matches(rows, prefix, n, e):
    return [r for r in rows if r["report_file"].startswith(prefix)
            and as_int(r.get("denominator")) == n and as_int(r.get("events")) == e]


def build(root):
    src = os.path.join(root, "data", "early_tox", "source")
    xdir = os.path.join(root, "docs", "development", "2026_11_10_early_toxicity", "extraction")
    arms_in = read_csv(os.path.join(src, "arm_curation.csv"))
    outs_in = read_csv(os.path.join(src, "outcome_curation.csv"))
    ex = load_extractions(xdir)

    arms, schedules, fractions, seen = [], [], [], set()
    by_arm = {}
    for a in arms_in:
        if a["arm_id"] in by_arm:
            raise ValueError(f"duplicate arm_id {a['arm_id']}")
        by_arm[a["arm_id"]] = a
        variants = []
        for v in a["pattern"].split("|"):
            fr, w = parse_variant(v.strip())
            sid = schedule_id(fr)
            variants.append(f"{sid}:{w:g}")
            if sid not in seen:
                seen.add(sid)
                last = fr[-1][1]
                schedules.append(dict(
                    schedule_id=sid, evidence_kind="trial_arm", n_fractions=len(fr),
                    total_dose_gy=round(sum(d for d, _ in fr), 3), last_time_h=last,
                    elapsed_days=round(last / 24, 4), inclusive_days=int(last // 24) + 1,
                    dose_time_status="planned (protocol); see arm pattern_notes",
                    clock="hours from the first fraction; day 1 is a Monday"))
                fractions += [dict(schedule_id=sid, fraction=k + 1, dose_gy=d, time_h=t)
                              for k, (d, t) in enumerate(fr)]
        arms.append(dict(arm_id=a["arm_id"], family_id=a["family_id"], trial_id=a["trial_id"],
                         report_prefix=a["report_prefix"], arm_label=a["arm_label"], role=a["role"],
                         context=a["context"], sites=a["sites"],
                         planned_systemic_therapy=a["planned_systemic_therapy"],
                         technique=a["technique"], pattern=a["pattern"],
                         schedule_variants=" ".join(variants), pattern_notes=a["pattern_notes"]))

    screening = {s["report_prefix"]: s["decision"] for s in read_csv(os.path.join(src, "screening_log.csv"))}
    outcomes = []
    for o in outs_in:
        decision = screening.get(o["report_prefix"])
        if decision not in ("included", "included_sensitivity", "benchmark"):
            raise ValueError(f"{o['outcome_id']}: report {o['report_prefix']} is screened as {decision!r}")
        arm = by_arm.get(o["arm_id"])
        if arm is None:
            raise ValueError(f"{o['outcome_id']}: unknown arm {o['arm_id']}")
        if arm["report_prefix"] != o["report_prefix"]:
            raise ValueError(f"{o['outcome_id']}: report differs from the arm's report")
        n, e = as_int(o["denominator"]), as_int(o["events"])
        if n is None or e is None or not 0 <= e <= n:
            raise ValueError(f"{o['outcome_id']}: invalid count {o['events']}/{o['denominator']}")
        hits = None
        for batch, (x1, x2) in ex.items():
            h1, h2 = matches(x1, o["report_prefix"], n, e), matches(x2, o["report_prefix"], n, e)
            if h1 or h2:
                if not (h1 and h2):
                    raise ValueError(f"{o['outcome_id']}: {e}/{n} found in only one extraction of batch {batch}")
                hits = (batch, h1[0], h2[0])
        if hits is None:
            raise ValueError(f"{o['outcome_id']}: {e}/{n} not found in any extraction of {o['report_prefix']}")
        batch, r1, r2 = hits
        outcomes.append(dict(
            outcome_id=o["outcome_id"], arm_id=o["arm_id"], family_id=arm["family_id"], role=arm["role"],
            endpoint_class=o["endpoint_class"], endpoint_as_written=r1.get("endpoint_as_written", ""),
            scale=r1.get("scale", ""), assessment=json.dumps(r1.get("assessment", {}), ensure_ascii=False),
            denominator=n, denominator_type=r1.get("denominator_type", ""), events=e,
            percent_reported=o["percent"], count_derivation=r1.get("count_derivation", ""),
            status=o["status"], value_quote=r1.get("value_quote", ""), value_location=r1.get("value_location", ""),
            report_file=r1["report_file"], extraction_agreement=f"batch {batch}: both extractions",
            note=o["note"]))

    families = defaultdict(lambda: dict(roles=set(), reports=set(), arms=0))
    for a in arms:
        f = families[a["family_id"]]
        f["roles"].add(a["role"])
        f["reports"].add(a["report_prefix"])
        f["arms"] += 1
    for s in read_csv(os.path.join(src, "screening_log.csv")):
        families[s["family_id"]]["reports"].add(s["report_prefix"])
    fam_rows = [dict(family_id=k, roles=";".join(sorted(v["roles"])) or "none",
                     reports=";".join(sorted(v["reports"])), n_arms=v["arms"])
                for k, v in sorted(families.items())]
    return arms, outcomes, schedules, fractions, fam_rows


def main(argv):
    root = argv[1] if len(argv) > 1 else os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    arms, outcomes, schedules, fractions, fams = build(root)
    out = os.path.join(root, "data", "early_tox")
    write_csv(os.path.join(out, "arms.csv"), arms)
    write_csv(os.path.join(out, "outcomes.csv"), outcomes)
    write_csv(os.path.join(out, "arm_schedules.csv"), schedules)
    write_csv(os.path.join(out, "arm_fractions.csv"), fractions)
    write_csv(os.path.join(out, "evidence_families.csv"), fams)
    print(f"{len(arms)} arms, {len(outcomes)} outcomes, {len(schedules)} schedules, "
          f"{len(fractions)} fractions, {len(fams)} evidence families")


if __name__ == "__main__":
    main(sys.argv)
