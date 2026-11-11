"""Build the historical early-toxicity tables from the Pardo-Montero 2021 supporting file.

Plan: docs/development/2026_11_10_early_toxicity_plan.md, step 2 (R1, R2).

Usage:
    python3 tools/early_tox/build_historical.py PATH/TO/PardoMontero_2021_SupportingInfo_108_schedules.xls [OUT_DIR]

The spreadsheet itself is CC BY-NC and is not committed. Its schedule facts (fraction doses,
fraction times, labels) are written to OUT_DIR (default data/early_tox/):

    schedules.csv  one row per distinct schedule (identity = full dose and relative-time vectors)
    fractions.csv  one row per fraction of each schedule
    labels.csv     one row per spreadsheet record, with original and adjudicated labels,
                   course scope and the primary / sensitivity selections
    reports.csv    the source report

Rules carried from the plan: no label is inferred, unknown labels are refused, partial courses are
kept as rows but excluded from both selections, and schedule identity never uses a coarse key.
"""
import csv
import hashlib
import os
import sys

# Original label -> (adjudicated label, counted as an event in the sensitivity selection)
LABELS = {
    "Tol": ("tolerable", False),
    "Intol": ("intolerable", True),
    "Disp_Intol": ("disputed_leaning_intolerable", True),
    "Disp_Tol": ("disputed_leaning_tolerable", False),
    "Disp": ("disputed_unresolved", None),
}

# Pardo-Montero 2021, section 2.D: the supporting file holds 108 analysed schedules "plus eight
# additional tolerable schedules that correspond to the first part of gap treatments that have
# not been used in the analysis". They are the final eight records of the file (spreadsheet rows
# 541-576, 1-based), all labeled Tol. Position, count and label were checked against the text on
# 2026-11-10 (verify_review_claims.py).
N_ANALYSED = 108
N_PARTIAL = 8
PARTIAL_BASIS = ("Pardo-Montero 2021 section 2.D: eight first-part-of-gap schedules in the "
                 "supporting file, not analysed; they are its final eight records")

REPORT = dict(
    report_id="PardoMontero2021",
    citation=("Pardo-Montero J, Parga-Pazos M, Fenwick JD. Classification of tolerable/intolerable "
              "mucosal toxicity of head-and-neck radiotherapy schedules with a biomathematical model "
              "of cell dynamics. Med Phys. 2021;48(7):4075-4084."),
    doi="10.1002/mp.14834",
    pmid="33704792",
    local_file="Custom_BED_calculator_v2/References/PardoMontero_2021_SupportingInfo_108_schedules.xls (git-ignored)",
    locator="Supporting information spreadsheet; text section 2.D",
    licence="CC BY-NC (supporting file); not committed",
    extraction="schedules, fraction times and labels read by tools/early_tox/build_historical.py",
)

DOSE_DP = 3
TIME_DP = 2


def adjudicate(label):
    """Map an original label to (adjudicated label, sensitivity event). Unknown labels are refused."""
    key = str(label).strip()
    if key not in LABELS:
        raise ValueError(f"unknown schedule label {label!r}: add it to LABELS with a source, do not guess")
    return LABELS[key]


def schedule_id(doses, rel_times):
    """Identity from the full fraction-dose and relative-time vectors (plan R2)."""
    key = ";".join(f"{d:.{DOSE_DP}f}@{t:.{TIME_DP}f}" for d, t in zip(doses, rel_times))
    return "S" + hashlib.sha256(key.encode()).hexdigest()[:10]


def read_records(xls_path):
    import pandas as pd

    x = pd.read_excel(xls_path, header=None)
    recs, i = [], 0
    while i < len(x):
        if str(x.iat[i, 0]).strip() == "Reference":
            doses = pd.to_numeric(x.iloc[i + 1, 1:], errors="coerce").dropna().tolist()
            times = pd.to_numeric(x.iloc[i + 2, 1:], errors="coerce").dropna().tolist()
            if len(doses) != len(times) or not doses:
                raise ValueError(f"row {i + 1}: {len(doses)} doses against {len(times)} times")
            if any(b < a for a, b in zip(times, times[1:])):
                raise ValueError(f"row {i + 1}: fraction times are not ordered")
            recs.append(dict(row=i + 1, ref=str(x.iat[i, 1]).strip(), doses=doses, times=times,
                             label=str(x.iat[i + 3, 1]).strip()))
            i += 4
        else:
            i += 1
    if len(recs) != N_ANALYSED + N_PARTIAL:
        raise ValueError(f"expected {N_ANALYSED + N_PARTIAL} records, found {len(recs)}")
    return recs


def build(recs):
    schedules, fractions, labels = {}, [], []
    for k, r in enumerate(recs, start=1):
        t0 = r["times"][0]
        rel = [round(t - t0, TIME_DP) for t in r["times"]]
        doses = [round(d, DOSE_DP) for d in r["doses"]]
        sid = schedule_id(doses, rel)
        if sid not in schedules:
            elapsed_h = rel[-1]
            schedules[sid] = dict(
                schedule_id=sid, evidence_kind="historical_label", n_fractions=len(doses),
                total_dose_gy=round(sum(doses), DOSE_DP), first_time_h=0.0, last_time_h=elapsed_h,
                elapsed_days=round(elapsed_h / 24, 4), inclusive_days=int(elapsed_h // 24) + 1,
                dose_time_status="as tabulated by Pardo-Montero (median values where reported; 6 h "
                                 "between same-day fractions when not stated)",
                clock="hours from the first fraction")
            fractions += [dict(schedule_id=sid, fraction=j + 1, dose_gy=d, time_h=t)
                          for j, (d, t) in enumerate(zip(doses, rel))]
        adj, event = adjudicate(r["label"])
        partial = k > N_ANALYSED
        if partial and adj != "tolerable":
            raise ValueError(f"record {k}: partial-course record is not labeled Tol")
        primary = (not partial) and adj in ("tolerable", "intolerable")
        sensitivity = (not partial) and event is not None
        labels.append(dict(
            record_id=f"PM{k:03d}", report_id=REPORT["report_id"], spreadsheet_row=r["row"],
            source_ref=r["ref"], schedule_id=sid, label_original=r["label"], label_adjudicated=adj,
            course_scope="partial" if partial else "complete",
            scope_basis=PARTIAL_BASIS if partial else "analysed set (Pardo-Montero 2021 section 2.D)",
            include_primary=int(primary), include_sensitivity=int(sensitivity),
            event_primary="" if not primary else int(adj == "intolerable"),
            event_sensitivity="" if not sensitivity else int(bool(event)),
            exclusion_reason=("partial course (first part of a gap schedule)" if partial else
                              "" if primary else
                              "disputed label: excluded from primary" + ("" if sensitivity else
                              " and sensitivity (no leaning given)"))))
    return list(schedules.values()), fractions, labels


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(argv):
    if len(argv) < 2:
        sys.exit(__doc__)
    out = argv[2] if len(argv) > 2 else os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "data", "early_tox")
    os.makedirs(out, exist_ok=True)
    schedules, fractions, labels = build(read_records(argv[1]))
    write_csv(os.path.join(out, "schedules.csv"), schedules)
    write_csv(os.path.join(out, "fractions.csv"), fractions)
    write_csv(os.path.join(out, "labels.csv"), labels)
    write_csv(os.path.join(out, "reports.csv"), [REPORT])
    print(f"{len(labels)} records, {len(schedules)} schedules, {len(fractions)} fractions -> {out}")


if __name__ == "__main__":
    main(sys.argv)
