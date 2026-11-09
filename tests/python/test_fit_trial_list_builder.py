"""tools/build_fit_trial_list.py (version 4.4.0, sources and notes corrected 2026-11-09): every fitted arm gets its reference and drug
fields from the extraction and PubMed records, and anything unmatched or ambiguous is refused."""
import csv, json, os, subprocess, sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BUILDER = os.path.join(ROOT, "tools", "build_fit_trial_list.py")

DATA = [dict(row=1, ref_no="DAHANCA67", arm="6 fractions per week (accelerated)", pmid="14511925"),
        dict(row=2, ref_no="RTOG8313", arm="67.2 Gy", pmid="2180866")]


def setup(tmp_path, extract, records):
    s = tmp_path / "search"; (s / "extraction").mkdir(parents=True)
    (s / "extraction" / "all_extract.json").write_text(json.dumps(extract))
    (s / "pubmed_records.json").write_text(json.dumps(records))
    (s / "supp_records.json").write_text("[]")
    d = tmp_path / "data.csv"
    with open(d, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(DATA[0])); w.writeheader(); w.writerows(DATA)
    return s, d


def run(tmp_path, s, d):
    out = tmp_path / "out.csv"
    p = subprocess.run([sys.executable, BUILDER, str(d), str(out), "--search-dir", str(s)],
                       capture_output=True, text=True)
    return p, out


REC = [dict(pmid="14511925", first_author="Overgaard", year="2003", journal="Lancet",
            title="Five compared with six fractions.", doi="10.1/x"),
       dict(pmid="2180866", first_author="Cox", year="1990", journal="IJROBP", title="Dose-response.", doi="")]
EX = [dict(pmid="14511925", arm="6 fractions per week (accelerated)", chemo="none"),
      dict(pmid="14511925", arm="6 fractions per week (accelerated)", chemo="none"),
      dict(pmid="2180866", arm="67.2 Gy", chemo="none (no chemotherapy within 6 weeks; prior chemo allowed earlier)")]


def test_builds_one_row_per_arm_with_reference_and_drugs(tmp_path):
    s, d = setup(tmp_path, EX, REC)
    p, out = run(tmp_path, s, d)
    assert p.returncode == 0, p.stderr
    rows = list(csv.DictReader(open(out)))
    assert [r["ref_no"] for r in rows] == ["DAHANCA67", "RTOG8313"]
    assert rows[0]["citation"] == "Overgaard et al. Five compared with six fractions. Lancet 2003."
    assert rows[0]["trial"] == "DAHANCA 6&7 (Overgaard 2003)"
    assert "nimorazole" in rows[0]["other_drugs"] and rows[1]["other_drugs"] == ""
    assert rows[0]["shuryak2018_ref"] == "49" and rows[1]["shuryak2018_ref"] == ""
    assert rows[1]["chemotherapy"] == "none (no chemotherapy within 6 weeks; earlier chemotherapy allowed)"
    # The outcome publication and, where it differs or a table was checked, the schedule source.
    assert rows[0]["outcome_source"] == "Overgaard et al. Lancet 2003 (PMID 14511925)"
    assert rows[0]["schedule_source"] == "" and rows[1]["schedule_source"] == ""


def test_schedule_sources_and_notes_follow_the_full_texts():
    """Review of 2026-11-09, item 2.1 and 2.2: CAIR period medians are reported (Skladowski 2000, Table 1),
    Jackson's times are from its Table 3 and its twice-daily dose was amended, and the Princess Margaret
    times are from Cummings 2007, Table 2."""
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import build_fit_trial_list as b
    cair = b.TRIAL_NOTES["17011446"]
    assert "not period medians" not in cair and "Table 1, p.102" in cair
    for k in ("36 (CAIR) and 51 (control) days at 2.0 Gy", "40 and 54 days at 1.8 Gy", "35 and 35, and 40 and 40",
              "from 60 to 76 Gy", "below 60 to 75.6 Gy", "Table 2, p.711"):
        assert k in cair, k
    jac = b.TRIAL_NOTES["9165135"]
    assert "not checked" not in jac
    for k in ("47.6 and 24.25 days (trial paper, Table 3, p.41)", "last ten of its 41 patients received 62.8 Gy",
              "reduced to 1.8 Gy", "one pooled outcome per arm"):
        assert k in jac, k
    assert "Cummings 2007, Table 2, p.12" in b.TRIAL_NOTES["32055717"]
    assert set(b.SCHEDULE_SOURCES) == {"17011446", "9165135", "32055717"}
    assert "PMID 10799721" in b.SCHEDULE_SOURCES["17011446"]
    assert "PMID 17920715" in b.SCHEDULE_SOURCES["32055717"]


def test_refuses_an_arm_without_an_extraction_row(tmp_path):
    s, d = setup(tmp_path, EX[:2], REC)
    p, _ = run(tmp_path, s, d)
    assert p.returncode != 0 and "RTOG8313 / 67.2 Gy: no extraction row" in p.stderr


def test_refuses_extraction_rows_that_disagree_on_chemotherapy(tmp_path):
    ex = EX[:1] + [dict(EX[1], chemo="cisplatin 100 mg/m2")] + EX[2:]
    s, d = setup(tmp_path, ex, REC)
    p, _ = run(tmp_path, s, d)
    assert p.returncode != 0 and "disagree on chemotherapy" in p.stderr


def test_refuses_a_trial_without_a_pubmed_record(tmp_path):
    s, d = setup(tmp_path, EX, REC[:1])
    p, _ = run(tmp_path, s, d)
    assert p.returncode != 0 and "no PubMed record" in p.stderr


@pytest.mark.skipif(not os.path.isdir(os.path.join(ROOT, "..", "References_refit_trials", "search", "extraction")),
                    reason="extraction folder not present")
def test_committed_list_is_what_the_builder_makes(tmp_path):
    p, out = run(tmp_path, os.path.join(ROOT, "..", "References_refit_trials", "search"),
                 os.path.join(ROOT, "data", "refit_trials_2026_10_25c_main.csv"))
    assert p.returncode == 0, p.stderr
    assert open(out).read() == open(os.path.join(ROOT, "data", "fit_trial_list.csv")).read()
