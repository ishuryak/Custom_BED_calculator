"""Early-toxicity plan step 2, historical labels (plan tests 1, 2, 3 and 5).

Runs on the committed tables in data/early_tox/, so it needs no copy of the CC BY-NC spreadsheet.
Expected counts were recomputed independently on 2026-11-10 by
docs/development/2026_11_10_early_toxicity/verification/verify_review_claims.py.
"""
import csv
import importlib.util
import os
from collections import Counter, defaultdict

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data", "early_tox")


def rows(name):
    with open(os.path.join(DATA, name), newline="") as f:
        return list(csv.DictReader(f))


@pytest.fixture(scope="module")
def tables():
    sched = {r["schedule_id"]: r for r in rows("schedules.csv")}
    frac = defaultdict(list)
    for r in rows("fractions.csv"):
        frac[r["schedule_id"]].append((int(r["fraction"]), float(r["dose_gy"]), float(r["time_h"])))
    return rows("labels.csv"), sched, frac


def is_event(r):
    return "intol" in r["label_original"].lower()


def test_raw_labels_and_selection_fixtures(tables):
    labels, _, _ = tables
    assert len(labels) == 116
    assert Counter(r["label_original"] for r in labels) == {
        "Tol": 98, "Intol": 13, "Disp_Tol": 2, "Disp_Intol": 2, "Disp": 1}
    complete = [r for r in labels if r["course_scope"] == "complete"]
    sel = {
        "all": labels,
        "drop_unresolved": [r for r in labels if r["label_original"] != "Disp"],
        "drop_all_disputed": [r for r in labels if not r["label_original"].startswith("Disp")],
        "sensitivity": [r for r in labels if r["include_sensitivity"] == "1"],
        "primary": [r for r in labels if r["include_primary"] == "1"],
    }
    got = {k: (len(v), sum(map(is_event, v))) for k, v in sel.items()}
    assert got == {"all": (116, 15), "drop_unresolved": (115, 15), "drop_all_disputed": (111, 13),
                   "sensitivity": (107, 15), "primary": (103, 13)}
    # the selections are exactly the plan's definitions, not just the right size
    assert {r["record_id"] for r in sel["sensitivity"]} == {
        r["record_id"] for r in complete if r["label_original"] != "Disp"}
    assert {r["record_id"] for r in sel["primary"]} == {
        r["record_id"] for r in complete if not r["label_original"].startswith("Disp")}
    for r in sel["primary"]:
        assert r["event_primary"] == str(int(r["label_original"] == "Intol"))


def test_partial_courses_never_enter_a_likelihood(tables):
    labels, _, _ = tables
    partial = [r for r in labels if r["course_scope"] == "partial"]
    assert [int(r["spreadsheet_row"]) for r in partial] == list(range(541, 577, 5))
    assert [r["source_ref"] for r in partial] == ["K64", "K87", "K87", "K87", "K87", "K88", "F15", "K65"]
    assert all(r["label_original"] == "Tol" for r in partial)
    assert all(r["include_primary"] == "0" and r["include_sensitivity"] == "0" for r in partial)


def test_full_vector_identity_separates_the_false_merges(tables):
    labels, sched, frac = tables
    pool = [r for r in labels if r["label_original"] != "Disp"]
    full, coarse = defaultdict(list), defaultdict(list)
    for r in pool:
        s = sched[r["schedule_id"]]
        full[r["schedule_id"]].append(r)
        coarse[(int(s["n_fractions"]), round(float(s["total_dose_gy"]), 2), int(s["inclusive_days"]))].append(r)
    assert (len(coarse), len(full)) == (92, 96)
    assert sum(any(map(is_event, g)) for g in full.values()) == 13
    assert sum(any(map(is_event, g)) for g in coarse.values()) == 13
    merged = sorted(sorted(r["source_ref"] for r in g) for g in coarse.values()
                    if len({r["schedule_id"] for r in g}) > 1)
    assert merged == [["F15", "K88"], ["K17", "K34"], ["K46", "K46"], ["K72", "K73"]]


def test_fraction_vectors_reconcile_with_schedules(tables):
    _, sched, frac = tables
    assert set(sched) == set(frac)
    for sid, s in sched.items():
        f = sorted(frac[sid])
        assert [k for k, _, _ in f] == list(range(1, len(f) + 1))
        assert len(f) == int(s["n_fractions"])
        assert sum(d for _, d, _ in f) == pytest.approx(float(s["total_dose_gy"]), abs=1e-6)
        times = [t for _, _, t in f]
        assert times[0] == 0 and times == sorted(times)
        assert times[-1] == pytest.approx(float(s["last_time_h"]))
        assert int(s["inclusive_days"]) == int(times[-1] // 24) + 1
        assert all(d > 0 for _, d, _ in f)


def test_unknown_label_is_refused():
    spec = importlib.util.spec_from_file_location(
        "build_historical", os.path.join(ROOT, "tools", "early_tox", "build_historical.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.adjudicate("Tol") == ("tolerable", False)
    with pytest.raises(ValueError, match="unknown schedule label"):
        mod.adjudicate("Tolerable?")
