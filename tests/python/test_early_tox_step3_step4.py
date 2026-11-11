"""Early-toxicity steps 3 and 4: endpoint mapping (M1-M4) and identification machinery."""
import csv
import importlib.util
import math
import os

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, "tools", "early_tox", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mp = load("build_mapping")
s4 = load("step4_identification")


def committed(name):
    with open(os.path.join(ROOT, "data", "early_tox", name), newline="") as f:
        return list(csv.DictReader(f))


def test_mapping_matches_the_approved_decisions():
    rows = mp.build(ROOT)
    assert [{k: str(v) for k, v in r.items()} for r in rows] == committed("endpoint_mapping.csv")
    prim = [r for r in rows if r["in_primary"]]
    assert len(prim) == 35 and len({r["family_id"] for r in prim}) == 15
    # M1: CTCAE v4 trials pooled (scale group G2)
    assert {r["family_id"] for r in prim if r["scale_group"] == "G2"} == {"Tandon", "Sikdar", "Shivananjappa"}
    # M2: W1-only sensitivity is the nine families with a stated window
    assert {r["family_id"] for r in rows if r["in_w1_sensitivity"]} == {
        "CAIR1", "EORTC22851", "GORTEC9902", "ORO9301", "PMH_HF", "RTOG0129", "RTOG9003", "Tandon", "Vienna"}
    # M3: CAIR primary uses the approximate EORTC III-IV count, severe confluent is secondary
    cair = {r["outcome_id"]: r["analysis_set"] for r in rows if r["family_id"] == "CAIR1"}
    assert cair == {"O018": "primary", "O019": "primary", "O020": "secondary", "O021": "secondary"}
    # M4 and the benchmark
    assert {r["family_id"] for r in rows if r["analysis_set"] == "sensitivity"} == {"IAEA_ACC", "Pinto", "Poddar"}
    assert {r["family_id"] for r in rows if r["analysis_set"] == "benchmark"} == {"HYPNO"}


def test_a_family_without_a_rule_is_refused(monkeypatch):
    monkeypatch.delitem(mp.FAMILY_RULES, "Mani")
    with pytest.raises(ValueError, match="Mani has no mapping rule"):
        mp.build(ROOT)


@pytest.fixture(scope="module")
def design():
    return s4.load_design()


def test_reference_arms_have_zero_contrast_and_every_family_has_one(design):
    rows, fams = design
    dE, _, ref = s4.contrasts(rows, fams, 10.0, s4.ev.fowler_rate(), 7.0)
    assert len(ref) == 15
    assert all(dE[i] == 0 for i in ref.values())
    assert sum(1 for x in dE if x != 0) >= 18


def test_onset_drops_out_at_7_days_but_not_at_20(design):
    # Every primary arm lasts more than 7 elapsed days, so Tk cancels from all contrasts (plan R3).
    rows, fams = design
    J7, _ = s4.jacobian(rows, fams, 10.0, s4.ev.fowler_rate(), 7.0)
    assert np.abs(J7[:, 2]).max() < 1e-8
    J20, _ = s4.jacobian(rows, fams, 10.0, s4.ev.fowler_rate(), 20.0)
    assert np.abs(J20[:, 2]).max() > 0.1          # Vienna V-CHART ends on day 16


def test_higher_alpha_beta_lowers_exposure_differences_of_larger_fractions(design):
    rows, fams = design
    J, keep = s4.jacobian(rows, fams, 10.0, s4.ev.fowler_rate(), 7.0)
    tandon = [k for k, i in enumerate(keep) if rows[i]["arm_id"] == "Tandon_SMART"][0]
    assert J[tandon, 0] < 0                        # 2.4 Gy against 2.0 Gy fractions


def test_map_fit_recovers_the_truth_from_noise_free_counts(design):
    rows, fams = design
    fam_idx = np.array([fams.index(x["family"]) for x in rows])
    N = np.array([x["N"] for x in rows])
    dE, _, _ = s4.contrasts(rows, fams, 10.0, s4.ev.fowler_rate(), 7.0)
    eta = np.full(len(fams), -0.5)
    y = N / (1 + np.exp(-(eta[fam_idx] + 0.1 * dE)))
    start = np.concatenate([np.zeros(len(fams)), [math.log(0.05), math.log(10), math.log(0.8)]])
    res, se, ok = s4.fit_map(rows, fams, 7.0, y, fam_idx, start)
    k = len(fams)
    assert ok and res.success
    # priors pull slightly; the shared parameters must still land near the truth
    assert math.exp(res.x[k]) == pytest.approx(0.1, rel=0.15)
    assert math.exp(res.x[k + 2]) == pytest.approx(s4.ev.fowler_rate(), rel=0.15)
    assert math.exp(res.x[k + 1]) == pytest.approx(10.0, rel=0.35)
