"""Early-toxicity evaluator (plan tests 6, 7 and 8).

Expected values were recomputed independently on 2026-11-10
(docs/development/2026_11_10_early_toxicity/verification/verify_review_claims.py).
"""
import csv
import importlib.util
import math
import os
from collections import defaultdict

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
spec = importlib.util.spec_from_file_location("evaluator", os.path.join(ROOT, "tools", "early_tox", "evaluator.py"))
ev = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ev)


def schedules(prefix):
    fr = defaultdict(list)
    with open(os.path.join(ROOT, "data", "early_tox", f"{prefix}fractions.csv"), newline="") as f:
        for r in csv.DictReader(f):
            fr[r["schedule_id"]].append((int(r["fraction"]), float(r["dose_gy"]), float(r["time_h"])))
    return {k: ([d for _, d, _ in sorted(v)], [t for _, _, t in sorted(v)]) for k, v in fr.items()}


def arm_schedule(arm_id):
    with open(os.path.join(ROOT, "data", "early_tox", "arms.csv"), newline="") as f:
        a = {r["arm_id"]: r for r in csv.DictReader(f)}[arm_id]
    (sid, _), = [v.split(":") for v in a["schedule_variants"].split()]
    return schedules("arm_")[sid]


def test_fowler_no_repair_fixtures_from_the_step2_arm_schedules():
    # 70/35 over 46 elapsed days, and HYPNO's two randomized arms, taken from data/early_tox.
    assert ev.early_bed(*arm_schedule("RTOG9003_SF")) == pytest.approx(53.10544, abs=5e-6)
    hfx = ev.early_bed(*arm_schedule("HYPNO_HFX"))
    nfx = ev.early_bed(*arm_schedule("HYPNO_NFX"))
    assert hfx == pytest.approx(55.86597, abs=5e-6)
    assert nfx == pytest.approx(55.43495, abs=5e-6)
    assert hfx - nfx == pytest.approx(0.43102, abs=1e-5)


def test_repair_fixture_on_the_historical_70_in_35_schedule():
    sched = schedules("")
    hits = [s for s in sched.values() if len(s[0]) == 35 and abs(sum(s[0]) - 70) < 1e-9]
    doses, times = hits[0]
    D = sum(doses)
    assert D + ev.repair_q(doses, times, 3.2) / 10 == pytest.approx(84.12426, abs=5e-6)
    assert D + ev.repair_q(doses, times, None) / 10 == pytest.approx(84.0, abs=1e-9)


def test_repair_stress_fixture_and_limits():
    doses, times = [8.0] * 5, [24.0 * k for k in range(5)]
    extra = (ev.repair_q(doses, times, 3.2) - ev.repair_q(doses, times, None)) / 10
    assert extra == pytest.approx(0.28402, abs=5e-6)
    # infinite separation tends to complete repair
    far = [1e6 * k for k in range(5)]
    assert ev.repair_q(doses, far, 3.2) == pytest.approx(ev.repair_q(doses, far, None))
    # coincident fractions act as one combined dose
    assert ev.repair_q([2.0, 3.0], [0.0, 0.0], 3.2) == pytest.approx(5.0 ** 2)


def test_onset_cancels_from_a_contrast_only_when_both_courses_pass_it():
    a = ev.regular_schedule(35, 2.0)          # 46 elapsed days
    b = ev.regular_schedule(20, 2.75)         # 25 elapsed days
    c = ev.regular_schedule(5, 3.0)           # 4 elapsed days
    diff = lambda x, y, Tk: ev.early_bed(*x, Tk=Tk) - ev.early_bed(*y, Tk=Tk)
    assert diff(a, b, 7.0) == pytest.approx(diff(a, b, 20.0), abs=1e-12)
    assert diff(a, c, 2.0) != pytest.approx(diff(a, c, 7.0), abs=1e-6)


def test_dose_per_fraction_raises_exposure_at_fixed_number_of_fractions():
    vals = [ev.early_bed(*ev.regular_schedule(20, d), halftime_h=3.2) for d in (2.0, 2.5, 3.0, 3.5)]
    assert vals == sorted(vals) and len(set(vals)) == 4


def test_invalid_inputs_are_refused():
    with pytest.raises(ValueError, match="not ordered"):
        ev.repair_q([2.0, 2.0], [24.0, 0.0], 3.2)
    with pytest.raises(ValueError, match="alpha/beta"):
        ev.early_bed([2.0], [0.0], ab=0)
    with pytest.raises(ValueError, match="half-time"):
        ev.repair_q([2.0, 2.0], [0.0, 6.0], 0)


def test_fowler_rate_matches_the_plan():
    assert ev.fowler_rate() == pytest.approx(math.log(2) / (0.35 * 2.5))
    assert ev.fowler_rate() == pytest.approx(0.79217, abs=5e-6)
