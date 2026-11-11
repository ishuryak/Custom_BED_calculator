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


def add_rows(src, arm_row, outcome_row, map_row=None):
    with open(src / "arm_curation.csv", "a") as f:
        f.write(arm_row + "\n")
    with open(src / "outcome_curation.csv", "a") as f:
        f.write(outcome_row + "\n")
    if map_row:
        with open(src / "extraction_map.csv", "a") as f:
            f.write(map_row + "\n")


def replace_in(src, name, old, new):
    p = src / name
    text = p.read_text()
    assert old in text
    p.write_text(text.replace(old, new))


def test_a_count_found_in_neither_extraction_is_refused(tmp_path):
    src = copy_sources(tmp_path)
    p = src / "outcome_curation.csv"
    # The Ghoshal conventional-arm conflict: the table says 27, the Discussion says 23.
    p.write_text(p.read_text().replace("O017,Ghoshal_CF,Ghoshal_S_2008,grade_3_4_mucositis,142,27,",
                                       "O017,Ghoshal_CF,Ghoshal_S_2008,grade_3_4_mucositis,142,23,"))
    with pytest.raises(ValueError, match=r"O017: extraction A1#\d+ holds 27/142, not 23/142"):
        load_builder().build(str(tmp_path))


def test_a_count_found_only_in_an_unaccepted_endpoint_row_is_refused(tmp_path):
    # Horiot 1992: only extraction A2 recorded a functional (patient-reported) mucosal G3-G4 row (64/158);
    # it is not an accepted equivalent of the objective grade 3-4 class (review R7). The one-extraction
    # refusal itself is tested with ELAN-RT (test_adjudicated_single_extraction_is_refused_for_a_candidate_arm).
    src = copy_sources(tmp_path)
    add_rows(src, "X_EO,EORTC22791,EORTC 22791,Horiot_JC_1992,CF,candidate,C1,x,none,x,W5 35d[2.0],x",
             "O900,X_EO,Horiot_JC_1992,grade_3_4_mucositis,158,64,,usable,x", "O900,A,,31,x")
    with pytest.raises(ValueError, match="O900: extraction A2#31 endpoint .*not an accepted equivalent"):
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
    o = dict(outcome_id="O068", report_prefix="Ortholan_C_2026", endpoint_class="grade_3_plus_mucositis")
    hits1 = [r for r in e1 if not b.row_problem(r, o, 91, 11)]
    hits2 = [r for r in e2 if not b.row_problem(r, o, 91, 11)]
    assert not hits1 and len(hits2) == 1


# Identity checks (review finding F7, 2026-11-11). Each must fail on identity, not on chance.
def test_a_count_moved_to_another_arms_value_is_refused(tmp_path):
    # Gupta P conventional arm 8/30 changed to 15/30, which is the boost + cisplatin arm's count.
    src = copy_sources(tmp_path)
    replace_in(src, "outcome_curation.csv", "O055,GuptaP_CFCCT,Gupta_P_2020,grade_3_mucositis,30,8,",
               "O055,GuptaP_CFCCT,Gupta_P_2020,grade_3_mucositis,30,15,")
    with pytest.raises(ValueError, match=r"O055: extraction D1#\d+ holds 8/30, not 15/30"):
        load_builder().build(str(tmp_path))


def test_one_extraction_row_cannot_confirm_two_outcomes(tmp_path):
    # A second outcome with the same arm, numbers and extraction rows as O055 (different class).
    src = copy_sources(tmp_path)
    m = {r.split(",")[0]: r for r in (src / "extraction_map.csv").read_text().splitlines()[1:]}
    with open(src / "outcome_curation.csv", "a") as f:
        f.write("O902,GuptaP_CFCCT,Gupta_P_2020,grade_3_4_mucositis,30,8,,usable,x\n")
    with open(src / "extraction_map.csv", "a") as f:
        f.write(",".join(["O902"] + m["O055"].split(",")[1:]) + "\n")
    with pytest.raises(ValueError, match=r"O902: extraction row D1#\d+ already confirms O055"):
        load_builder().build(str(tmp_path))


def test_only_standard_pairings_or_recorded_equivalences_confirm_a_count():
    # Review R7: grade 1, grade 2, any grade, grade 4 alone and pharyngeal rows never confirm a severe
    # count; a confluent class is not confirmed by a grade >= 3 code (or the reverse) without a record.
    b = load_builder()
    b.EQUIVALENCES.clear()
    b.EQUIVALENCES.update(b.load_equivalences(os.path.join(ROOT, "data", "early_tox", "source")))
    o = dict(outcome_id="X", report_prefix="Gupta_P_2020", endpoint_class="grade_3_mucositis")
    row = dict(report_file="Gupta_P_2020_x.pdf", denominator=30, events=8, endpoint="max_grade_ge3_mucositis")
    assert b.row_problem(row, o, 30, 8) == ""
    for bad in ("other:grade 1 mucositis", "other:mucositis grade 2", "other:mucositis any grade",
                "other:mucositis grade 4", "other:grade 3-4 mucositis", "confluent_mucositis"):
        assert "not an accepted equivalent" in b.row_problem(dict(row, endpoint=bad), o, 30, 8), bad
    assert "not mucositis" in b.row_problem(dict(row, endpoint="other:pharyngeal grade 3"), o, 30, 8)
    assert "is from" in b.row_problem(dict(row, report_file="Tallari_RV_2017_x.pdf"), o, 30, 8)
    oc = dict(outcome_id="Y", report_prefix="Gupta_M_2015", endpoint_class="confluent_mucositis")
    assert b.row_problem(dict(row, report_file="Gupta_M_2015_x.pdf", endpoint="confluent_mucositis"), oc, 30, 8) == ""
    assert "not an accepted equivalent" in b.row_problem(dict(row, report_file="Gupta_M_2015_x.pdf"), oc, 30, 8)
    # a recorded equivalence holds for its own report only
    ep = "other:grade>=3 mucositis at treatment completion (CTCAE v4)"
    op = dict(outcome_id="Z", report_prefix="Poddar_J_2017", endpoint_class="grade_3_mucositis_end_of_rt")
    assert b.row_problem(dict(row, report_file="Poddar_J_2017_x.pdf", endpoint=ep), op, 30, 8) == ""
    assert "not an accepted equivalent" in b.row_problem(dict(row, endpoint=ep), dict(op, report_prefix="Gupta_P_2020"), 30, 8)


def test_an_equivalence_without_a_basis_is_refused(tmp_path):
    src = copy_sources(tmp_path)
    with open(src / "endpoint_equivalence.csv", "a") as f:
        f.write("Gupta_P_2020,other:grade 1 mucositis,grade_3_mucositis, ,2026-11-11\n")
    with pytest.raises(ValueError, match="endpoint equivalence without a basis"):
        load_builder().build(str(tmp_path))


def test_duplicate_outcome_ids_and_map_rows_are_refused(tmp_path):
    src = copy_sources(tmp_path)
    lines = (src / "outcome_curation.csv").read_text().splitlines()
    dup = next(l for l in lines if l.startswith("O055,"))
    with open(src / "outcome_curation.csv", "a") as f:
        f.write(dup + "\n")
    with pytest.raises(ValueError, match="O055: duplicate outcome_id"):
        load_builder().build(str(tmp_path))
    src2 = copy_sources(tmp_path / "b")
    mline = next(l for l in (src2 / "extraction_map.csv").read_text().splitlines() if l.startswith("O055,"))
    with open(src2 / "extraction_map.csv", "a") as f:
        f.write(mline + "\n")
    with pytest.raises(ValueError, match="O055: duplicate row in extraction_map.csv"):
        load_builder().build(str(tmp_path / "b"))


def test_every_outcome_records_the_rows_that_confirm_it(built):
    _, outcomes, _, _, _ = built
    for o in outcomes:
        rows = o["extraction_rows"].split(";")
        assert len(rows) == (1 if o["outcome_id"] in ("O067", "O068") else 2), o["outcome_id"]


def test_duplicate_reports_add_no_counts_even_when_both_extractions_agree(tmp_path):
    # Maciejewski 1996 is the first 44 CAIR patients: 11/23 is in both extractions, but the
    # screening log marks the report a duplicate, so the builder must refuse it.
    src = copy_sources(tmp_path)
    add_rows(src, "X_MAC,CAIR1,CAIR,Maciejewski_B_1996,Arm A,candidate,C1,x,none,x,W7 35d[2.0],x",
             "O901,X_MAC,Maciejewski_B_1996,confluent_mucositis,23,11,,usable,x", "O901,A,0,0,x")
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
