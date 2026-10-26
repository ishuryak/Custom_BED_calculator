"""Provenance table beside a refit data file (review R08, 2026-10-26).

Writes data/<data file stem>_provenance.csv, one row per fitted arm, keyed by row, ref_no and arm,
from the same extraction rows the builder matched. The fitting file itself is not changed, so no
fit is invalidated: the table records what the fit does not use.

Columns: the exact source value and whether it is a failure rate; the source interval as
extracted, with its confidence level when stated and a reason when absent; estimator, endpoint,
time point and source location; the quote-check status; where the patient number came from;
population exceptions within the arm; prescription categories and whether their patient weights
are known; and the rule that chose the representative schedule.

Population and prescription notes were read in the trial publications (2026-10-26):
- RTOG 83-13 (Cox 1990, PMID 2180866), Table 2: the analyzed cases include nasopharynx and sinus sites.
- DAHANCA 6/7 (Overgaard 2003, PMID 14511925): eligibility included the nasopharynx; all patients
  except those with glottic cancer received nimorazole; 66 Gy (tumor and nodes <= 4 cm) or 68 Gy,
  and well-differentiated stage I glottic tumors could receive 62 Gy (31 fractions). The share of
  patients in each category is not given, so the representative schedule is an assumption.

Run from the repository root:  python3 fit/build_provenance.py data/refit_trials_2026_10_25c_main.csv [builder flags]
"""
import csv, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_refit_dataset as B

POPULATION = {
    "RTOG8313": "includes nasopharyngeal and sinus cases (Cox 1990, Table 2)",
    "DAHANCA67": "eligibility included the nasopharynx; nonglottic patients received nimorazole (Overgaard 2003)",
}
PRESCRIPTION = {
    "DAHANCA67": ("66 Gy (tumor and nodes <= 4 cm) or 68 Gy; 62 Gy in 31 fractions allowed for well-differentiated "
                  "stage I glottic tumors (Overgaard 2003)", "unknown"),
}


def ci_fields(ci):
    ci = (ci or "").strip()
    if not ci:
        return "", "", "not extracted (an interval may exist in the source)"
    m = re.search(r"(\d+(?:\.\d+)?)\s*%\s*CI", ci)
    return ci, (m.group(1) + "%" if m else "not stated in the extraction"), ""


def quote_status(r):
    q = r.get("_quote_ok") or {}
    if r.get("value_source") in ("figure_read", "table_pdf_only"):
        return f"no text quote: {r.get('value_source')}"
    if not q:
        return "no quote check recorded"
    bad = [k for k, v in q.items() if not v]
    return "all quotes matched" if not bad else "failed: " + ", ".join(bad)


def provenance(argv):
    matches = []
    out, problems = B.build(argv, matches)
    if problems:
        sys.exit("REFUSED:\n" + "\n".join(problems))
    if "--km-only" in argv:     # the same filter the builder applies, on the same row objects
        keep = {id(o) for o in B.km_only(out)}
        matches = [(o, r) for o, r in matches if id(o) in keep]
    rows = []
    for o, r in matches:
        ci, level, missing = ci_fields(r.get("value_ci"))
        presc, weights = PRESCRIPTION.get(o["first_author"], ("", ""))
        rows.append(dict(row=o["row"], ref_no=o["ref_no"], first_author=o["first_author"], arm=o["arm"], pmid=o["pmid"],
                         source_value_pct=o["source_value_pct"], source_value_is_failure=o["source_value_is_failure"],
                         source_ci=ci, ci_level=level, ci_missing_reason=missing,
                         estimate_type=o["estimate_type"], endpoint=o["endpoint"], timepoint_years=o["timepoint_years"],
                         value_source=o["value_source"], location=o["location"], quote_status=quote_status(r),
                         n_source=o["n_source"], population_exceptions=POPULATION.get(o["first_author"], ""),
                         prescription_categories=presc, prescription_weights=weights,
                         representative_schedule_rule=o["T_source_rule"] + ("; " + o["notes"] if o["notes"] else "")))
    return rows


def main(argv):
    rows = provenance(argv)
    path = re.sub(r"\.csv$", "_provenance.csv", argv[1])
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(len(rows), "arms,", sum(bool(r["source_ci"]) for r in rows), "with a source interval ->", path)


if __name__ == "__main__":
    main(sys.argv)
