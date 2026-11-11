"""Build the randomized-arm early-toxicity tables (plan step 2).

Usage:
    python3 tools/early_tox/build_arms.py [REPO_ROOT]

Inputs (committed):
    data/early_tox/source/arm_curation.csv       one row per arm: role, context, therapy, schedule pattern
    data/early_tox/source/outcome_curation.csv   one row per arm outcome: expected denominator and events
    data/early_tox/source/extraction_map.csv     the exact row of each extraction that confirms each outcome
    docs/development/2026_11_10_early_toxicity/extraction/out_{A,B,C,D,E}{1,2}.json
                                                 two independent extractions per batch of reports
Outputs (data/early_tox/):
    arms.csv, outcomes.csv, arm_schedules.csv, arm_fractions.csv, evidence_families.csv

An outcome is written only if BOTH independent extractions confirm it through the rows named in
extraction_map.csv: same report, same denominator and event count, a mucositis endpoint compatible
with the outcome's class and never a grade 2 threshold, and no row confirming two outcomes
(identity matching added 2026-11-11, review finding F7). Quotes, locations, scales and assessment windows are copied from the
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
EXTRACTIONS = {"A": ("out_A1.json", "out_A2.json"), "B": ("out_B1.json", "out_B2.json"),
               "C": ("out_C1.json", "out_C2.json"), "D": ("out_D1.json", "out_D2.json"),
               "E": ("out_E1.json", "out_E2.json")}
# A count found in only one extraction is accepted only for a sensitivity-only arm, with
# status "adjudicated_single" and a note saying how it was confirmed (added 2026-11-11 for the
# ELAN-RT row hidden under a print banner).
SINGLE_OK_STATUS = "adjudicated_single"


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


def is_mucositis_endpoint(endpoint):
    """Only mucositis rows can confirm a count (a pharyngeal or feeding-tube row with the same
    numbers is a coincidence, found 2026-11-11 for ELAN-RT 11/91)."""
    ep = (endpoint or "").lower()
    return ep in ("max_grade_ge3_mucositis", "confluent_mucositis") or (ep.startswith("other:") and "mucos" in ep)


# Accepted endpoint pairings (review R7, 2026-11-11): a curated outcome class is confirmed only by
# its standard extraction code, or by a source-specific adjudicated equivalence recorded in
# data/early_tox/source/endpoint_equivalence.csv. Nothing is inferred from the word "mucositis".
STANDARD_ENDPOINTS = {
    "grade_3_4_mucositis": "max_grade_ge3_mucositis", "grade_3_mucositis": "max_grade_ge3_mucositis",
    "grade_3_plus_mucositis": "max_grade_ge3_mucositis", "grade_3_4_mucositis_from_dische": "max_grade_ge3_mucositis",
    "diffuse_mucositis_top_grade": "max_grade_ge3_mucositis", "grade_3_mucositis_end_of_rt": "max_grade_ge3_mucositis",
    "confluent_mucositis": "confluent_mucositis", "confluent_mucositis_own_scale": "confluent_mucositis",
    "severe_confluent_dische": "confluent_mucositis", "diffuse_mucositis": "confluent_mucositis",
}
EQUIVALENCES = set()


def load_equivalences(src):
    out = set()
    path = os.path.join(src, "endpoint_equivalence.csv")
    if os.path.exists(path):
        for r in read_csv(path):
            if not r["basis"].strip():
                raise ValueError(f"{r['report_prefix']}: endpoint equivalence without a basis")
            out.add((r["report_prefix"], r["extraction_endpoint"], r["outcome_class"]))
    return out


def row_problem(r, o, n, e):
    """Why an extraction row cannot confirm outcome o (n events e), or '' if it can."""
    if not r["report_file"].startswith(o["report_prefix"]):
        return f"is from {r['report_file'][:40]}, not {o['report_prefix']}"
    if as_int(r.get("denominator")) != n or as_int(r.get("events")) != e:
        return f"holds {r.get('events')}/{r.get('denominator')}, not {e}/{n}"
    ep = (r.get("endpoint") or "")
    cls = o["endpoint_class"]
    if STANDARD_ENDPOINTS.get(cls) == ep or (o["report_prefix"], ep, cls) in EQUIVALENCES:
        return ""
    if not is_mucositis_endpoint(ep):
        return f"endpoint {ep!r} is not mucositis"
    return f"endpoint {ep!r} is not an accepted equivalent of class {cls} (no standard pairing or recorded equivalence)"


def build(root):
    src = os.path.join(root, "data", "early_tox", "source")
    xdir = os.path.join(root, "docs", "development", "2026_11_10_early_toxicity", "extraction")
    arms_in = read_csv(os.path.join(src, "arm_curation.csv"))
    outs_in = read_csv(os.path.join(src, "outcome_curation.csv"))
    ex = load_extractions(xdir)
    EQUIVALENCES.clear()
    EQUIVALENCES.update(load_equivalences(src))

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
    xmap = {}
    for m in read_csv(os.path.join(src, "extraction_map.csv")):
        if m["outcome_id"] in xmap:
            raise ValueError(f"{m['outcome_id']}: duplicate row in extraction_map.csv")
        xmap[m["outcome_id"]] = m
    outcomes, seen_outcomes, seen_arm_endpoint, used_rows = [], set(), set(), {}
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
        if o["outcome_id"] in seen_outcomes:
            raise ValueError(f"{o['outcome_id']}: duplicate outcome_id")
        seen_outcomes.add(o["outcome_id"])
        key = (o["arm_id"], o["endpoint_class"])
        if key in seen_arm_endpoint:
            raise ValueError(f"{o['outcome_id']}: second outcome for arm {key[0]} and endpoint {key[1]}")
        seen_arm_endpoint.add(key)
        mp = xmap.get(o["outcome_id"])
        if mp is None:
            raise ValueError(f"{o['outcome_id']}: no row in extraction_map.csv")
        batch = mp["batch"]
        if batch not in ex:
            raise ValueError(f"{o['outcome_id']}: unknown extraction batch {batch!r}")
        picked = []
        for k, col in ((0, "row_extraction_1"), (1, "row_extraction_2")):
            if mp[col] == "":
                picked.append(None)
                continue
            idx = int(mp[col])
            rows_k = ex[batch][k]
            if not 0 <= idx < len(rows_k):
                raise ValueError(f"{o['outcome_id']}: extraction {batch}{k + 1} has no row {idx}")
            r = rows_k[idx]
            tag = (batch, k, idx)
            if tag in used_rows:
                raise ValueError(f"{o['outcome_id']}: extraction row {batch}{k + 1}#{idx} already confirms {used_rows[tag]}")
            used_rows[tag] = o["outcome_id"]
            problem = row_problem(r, o, n, e)
            if problem:
                raise ValueError(f"{o['outcome_id']}: extraction {batch}{k + 1}#{idx} {problem}")
            picked.append((idx, r))
        if picked[0] is None and picked[1] is None:
            raise ValueError(f"{o['outcome_id']}: {e}/{n} not found in any extraction of {o['report_prefix']}")
        if picked[0] is None or picked[1] is None:
            single_ok = (o["status"] == SINGLE_OK_STATUS and arm["role"] == "sensitivity" and len(o["note"]) > 40)
            if not single_ok:
                raise ValueError(f"{o['outcome_id']}: {e}/{n} found in only one extraction of batch {batch}")
        first = picked[0] or picked[1]
        r1, r2 = first[1], (picked[1][1] if picked[0] and picked[1] else None)
        rows_used = ";".join(f"{batch}{k + 1}#{p[0]}" for k, p in enumerate(picked) if p)
        outcomes.append(dict(
            outcome_id=o["outcome_id"], arm_id=o["arm_id"], family_id=arm["family_id"], role=arm["role"],
            endpoint_class=o["endpoint_class"], endpoint_as_written=r1.get("endpoint_as_written", ""),
            scale=r1.get("scale", ""), assessment=json.dumps(r1.get("assessment", {}), ensure_ascii=False),
            denominator=n, denominator_type=r1.get("denominator_type", ""), events=e,
            percent_reported=o["percent"], count_derivation=r1.get("count_derivation", ""),
            status=o["status"], value_quote=r1.get("value_quote", ""), value_location=r1.get("value_location", ""),
            report_file=r1["report_file"],
            extraction_agreement=(f"batch {batch}: both extractions" if r2 is not None
                                  else f"batch {batch}: one extraction, adjudicated (see note)"),
            extraction_rows=rows_used,
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
