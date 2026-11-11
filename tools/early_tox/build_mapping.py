"""Apply the step 3 endpoint mapping (decisions M1-M4, approved 2026-11-10) to data/early_tox/outcomes.csv.

Usage: python3 tools/early_tox/build_mapping.py [REPO_ROOT]
Writes data/early_tox/endpoint_mapping.csv, one row per outcome.

Scale groups:
    G1  grade 3 defined by mucosal appearance as reported (RTOG/EORTC, NCI CTC v2, EORTC objective
        top grade, Dische-based confluent)
    G2  CTCAE v4: grade 3 = "Severe pain; interfering with oral intake" (checked against the NCI
        CTCAE v4.03 quick reference on 2026-11-10)
Windows:
    W1  stated and reaching at least 4 weeks after RT (or 90 days from the start of a 7-week
        course), or follow-up to healing
    W2  not stated, or ending during or within 2 weeks of the end of RT
Analysis sets:
    primary      M1: G1 and G2 pooled, with a scale-family offset. M2: W1 and W2. M3: CAIR uses the
                 approximate EORTC III-IV count
    secondary    CAIR severe confluent (Dische score above 15)
    sensitivity  M4: IAEA-ACC, Pinto 1991, Poddar 2017
    benchmark    HYPNO (plan R6), never fitted
"""
import csv
import os
import sys

# family -> (scale group, window)
FAMILY_RULES = {
    "GORTEC9902": ("G1", "W1"), "RTOG9003": ("G1", "W1"), "RTOG0129": ("G1", "W1"),
    "EORTC22851": ("G1", "W1"), "PMH_HF": ("G1", "W1"), "CAIR1": ("G1", "W1"),
    "ORO9301": ("G1", "W1"), "Vienna": ("G1", "W1"), "Tandon": ("G2", "W1"),
    "EORTC22791": ("G1", "W2"), "Ghoshal": ("G1", "W2"), "Bartelink2002": ("G1", "W2"),
    "Mani": ("G1", "W2"), "Sikdar": ("G2", "W2"), "Shivananjappa": ("G2", "W2"),
    "IAEA_ACC": ("G1", "none"), "Pinto": ("other", "W2"), "Poddar": ("G2", "W2"),
    "HYPNO": ("G2", "W2"),
}
SECONDARY_CLASSES = {"severe_confluent_dische"}


def analysis_set(o):
    if o["role"] == "benchmark":
        return "benchmark"
    if o["role"] == "sensitivity":
        return "sensitivity"
    if o["role"] != "candidate":
        raise ValueError(f"{o['outcome_id']}: unknown role {o['role']!r}")
    return "secondary" if o["endpoint_class"] in SECONDARY_CLASSES else "primary"


def build(root):
    with open(os.path.join(root, "data", "early_tox", "outcomes.csv"), newline="") as f:
        outcomes = list(csv.DictReader(f))
    rows = []
    for o in outcomes:
        if o["family_id"] not in FAMILY_RULES:
            raise ValueError(f"{o['outcome_id']}: family {o['family_id']} has no mapping rule")
        group, window = FAMILY_RULES[o["family_id"]]
        aset = analysis_set(o)
        rows.append(dict(outcome_id=o["outcome_id"], arm_id=o["arm_id"], family_id=o["family_id"],
                         endpoint_class=o["endpoint_class"], scale_group=group, window=window,
                         analysis_set=aset,
                         in_primary=int(aset == "primary"),
                         in_w1_sensitivity=int(aset == "primary" and window == "W1")))
    return rows


def main(argv):
    root = argv[1] if len(argv) > 1 else os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    rows = build(root)
    with open(os.path.join(root, "data", "early_tox", "endpoint_mapping.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} outcomes mapped; primary {sum(r['in_primary'] for r in rows)}, "
          f"W1 sensitivity {sum(r['in_w1_sensitivity'] for r in rows)}")


if __name__ == "__main__":
    main(sys.argv)
