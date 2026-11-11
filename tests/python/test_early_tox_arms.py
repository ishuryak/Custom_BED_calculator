"""Early-toxicity plan step 2, randomized-arm tables (plan tests 1, 4, 5 and 12).

Runs the builder into a temporary copy, so it checks the committed sources, and compares the
result with the committed tables.
"""
import csv
import importlib.util
import os
import shutil

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def load_builder():
    spec = importlib.util.spec_from_file_location(
        "build_arms", os.path.join(ROOT, "tools", "early_tox", "build_arms.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rows(name):
    with open(os.path.join(ROOT, "data", "early_tox", name), newline="") as f:
        return list(csv.DictReader(f))


@pytest.fixture(scope="module")
def built():
    return load_builder().build(ROOT)


def test_committed_tables_match_a_fresh_build(built):
    arms, outcomes, schedules, fractions, fams = built
    for name, fresh in [("arms.csv", arms), ("outcomes.csv", outcomes), ("arm_schedules.csv", schedules),
                        ("evidence_families.csv", fams)]:
        assert [{k: str(v) for k, v in r.items()} for r in fresh] == rows(name), name
    assert len(rows("arm_fractions.csv")) == len(fractions)


def test_every_count_is_in_both_independent_extractions(built):
    _, outcomes, _, _, _ = built
    assert len(outcomes) == 68
    single = {o["outcome_id"] for o in outcomes if not o["extraction_agreement"].endswith("both extractions")}
    assert single == {"O067", "O068"}            # ELAN-RT, adjudicated, sensitivity only
    assert all(o["role"] == "sensitivity" for o in outcomes if o["outcome_id"] in single)
    assert all(o["value_quote"] and o["value_location"] for o in outcomes)
    assert all(0 <= int(o["events"]) <= int(o["denominator"]) for o in outcomes)


def copy_sources(tmp_path):
    for sub in ("data/early_tox/source", "docs/development/2026_11_10_early_toxicity/extraction"):
        shutil.copytree(os.path.join(ROOT, sub), tmp_path / sub)
    return tmp_path / "data/early_tox/source"


def add_rows(src, arm_row, outcome_row):
    with open(src / "arm_curation.csv", "a") as f:
        f.write(arm_row + "\n")
    with open(src / "outcome_curation.csv", "a") as f:
        f.write(outcome_row + "\n")


def test_a_count_found_in_neither_extraction_is_refused(tmp_path):
    src = copy_sources(tmp_path)
    p = src / "outcome_curation.csv"
    # The Ghoshal conventional-arm conflict: the table says 27, the Discussion says 23.
    p.write_text(p.read_text().replace("O017,Ghoshal_CF,Ghoshal_S_2008,grade_3_4_mucositis,142,27,",
                                       "O017,Ghoshal_CF,Ghoshal_S_2008,grade_3_4_mucositis,142,23,"))
    with pytest.raises(ValueError, match="O017: 23/142 not found"):
        load_builder().build(str(tmp_path))


def test_a_count_found_in_only_one_extraction_is_refused(tmp_path):
    # Horiot 1992: only extraction A1 recorded a functional mucosal G3-G4 row (64/158).
    src = copy_sources(tmp_path)
    add_rows(src, "X_EO,EORTC22791,EORTC 22791,Horiot_JC_1992,CF,candidate,C1,x,none,x,W5 35d[2.0],x",
             "O900,X_EO,Horiot_JC_1992,grade_3_4_mucositis,158,64,,usable,x")
    with pytest.raises(ValueError, match="O900: 64/158 found in only one extraction"):
        load_builder().build(str(tmp_path))


def test_adjudicated_single_extraction_is_refused_for_a_candidate_arm(tmp_path):
    # The ELAN-RT exception holds only for sensitivity arms with an explanatory note.
    src = copy_sources(tmp_path)
    p = src / "arm_curation.csv"
    p.write_text(p.read_text().replace("ELAN_HSC,ELAN,ELAN-RT,Ortholan_C_2026,HSC-RT 55 Gy / 20 split course,sensitivity",
                                       "ELAN_HSC,ELAN,ELAN-RT,Ortholan_C_2026,HSC-RT 55 Gy / 20 split course,candidate"))
    with pytest.raises(ValueError, match="O068: 11/91 found in only one extraction"):
        load_builder().build(str(tmp_path))


def test_only_mucositis_rows_confirm_a_count():
    # E1 has a pharyngeal grade 3 row of 11/91 for ELAN-RT; it must not confirm the mucositis count.
    b = load_builder()
    assert b.is_mucositis_endpoint("confluent_mucositis") and b.is_mucositis_endpoint("other:Mucositis grade 4")
    assert not b.is_mucositis_endpoint("other:pharyngeal_toxicity_grade3")
    assert not b.is_mucositis_endpoint("feeding_tube")
    ex = b.load_extractions(os.path.join(ROOT, "docs", "development", "2026_11_10_early_toxicity", "extraction"))
    e1, e2 = ex["E"]
    assert not b.matches(e1, "Ortholan_C_2026", 91, 11) and b.matches(e2, "Ortholan_C_2026", 91, 11)


def test_duplicate_reports_add_no_counts_even_when_both_extractions_agree(tmp_path):
    # Maciejewski 1996 is the first 44 CAIR patients: 11/23 is in both extractions, but the
    # screening log marks the report a duplicate, so the builder must refuse it.
    src = copy_sources(tmp_path)
    add_rows(src, "X_MAC,CAIR1,CAIR,Maciejewski_B_1996,Arm A,candidate,C1,x,none,x,W7 35d[2.0],x",
             "O901,X_MAC,Maciejewski_B_1996,confluent_mucositis,23,11,,usable,x")
    with pytest.raises(ValueError, match="O901: report Maciejewski_B_1996 is screened as 'duplicate'"):
        load_builder().build(str(tmp_path))


DUPLICATE_REPORTS = ("Maciejewski_B_1996", "Skladowski_K_2006", "Bernier_J_1998", "Dobrowsky_W_1998",
                     "Dobrowsky_W_2000", "Skladowski_K_2013")


def assert_no_duplicate_reports(outcomes):
    for o in outcomes:
        assert not o["report_file"].startswith(DUPLICATE_REPORTS), o["report_file"]


@pytest.mark.parametrize("arm,total,n,days", [
    ("RTOG0129_SFX", 70.0, 35, 47),        # paper: 47 days if starting Monday
    ("RTOG0129_AFXC", 72.0, 42, 40),       # paper: 40 days if starting Monday
    ("HYPNO_HFX", 55.0, 20, 26),           # paper: 12 days shorter than NFX if both start Monday
    ("HYPNO_NFX", 66.0, 33, 38),
    ("Vienna_VCHART", 55.3, 33, 17),       # paper: days 1-17
    ("GORTEC9902_ACCCT", 70.0, 40, 40),    # 70 Gy in 6 weeks
    ("GORTEC9902_VACC", 64.8, 36, 24),     # 64.8 Gy in 3.5 weeks
    ("RTOG9003_AFXS", 67.2, 42, 43),       # 6 weeks including a 2-week rest
    ("Bartelink_MFD", 72.0, 45, 47),       # weeks 1, 4 and 7 only
    ("Ghoshal_CB", 67.5, 40, 33),
])
def test_schedule_patterns_reproduce_the_protocols(built, arm, total, n, days):
    arms, _, schedules, _, _ = built
    sched = {r["schedule_id"]: r for r in schedules}
    a = {r["arm_id"]: r for r in arms}[arm]
    (sid, _), = [v.split(":") for v in a["schedule_variants"].split()]
    s = sched[sid]
    assert (float(s["total_dose_gy"]), int(s["n_fractions"]), int(s["inclusive_days"])) == (
        pytest.approx(total), n, days)


def test_fraction_times_and_same_day_intervals():
    b = load_builder()
    fr, _ = b.parse_variant("W5 I6 18d[1.8] 12d[1.8;1.5]")
    assert fr[18 - 1][1] == 17 * 24 + 3 * 48        # 18th weekday is Wed of week 4
    boost = [(d, t) for d, t in fr if d == 1.5]
    assert len(boost) == 12 and all((t % 24) == 6 for _, t in boost)
    times = [t for _, t in fr]
    assert times == sorted(times)
    assert all(t % 24 not in range(1, 6) for t in times)


def test_unknown_pattern_token_is_refused():
    with pytest.raises(ValueError, match="unknown schedule token"):
        load_builder().parse_variant("W5 35x2.0")


def test_benchmark_is_never_a_fitting_candidate_and_duplicates_add_no_counts(built):
    arms, outcomes, _, _, fams = built
    assert {a["role"] for a in arms if a["family_id"] == "HYPNO"} == {"benchmark"}
    assert all(";" not in f["roles"] for f in fams), "a family mixes roles"
    assert_no_duplicate_reports(outcomes)


def test_hypno_counts_match_the_verified_benchmark(built):
    _, outcomes, _, _, _ = built
    o = {r["arm_id"]: (r["events"], r["denominator"]) for r in outcomes}
    assert o["HYPNO_HFX"] == (190, 373) and o["HYPNO_NFX"] == (208, 380)
