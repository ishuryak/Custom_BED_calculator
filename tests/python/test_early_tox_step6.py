"""Early-toxicity step 6: leave-family-out predictive scoring, bootstrap and the HYPNO benchmark."""
import importlib.util
import math
import os

import numpy as np
import pytest
from scipy.stats import binom

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_spec = importlib.util.spec_from_file_location("s6", os.path.join(ROOT, "tools", "early_tox", "step6_validate.py"))
s6 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s6)


def test_binomial_log_pmf_matches_scipy():
    y, n, p = np.array([0, 3, 10]), np.array([10, 10, 10]), np.array([0.2, 0.3, 0.9])
    assert np.allclose(s6.binom_logpmf(y, n, p), binom.logpmf(y, n, p))


def fake_draws(S=50, eta=-0.5, kappa=0.1, gamma=0.4):
    return dict(m0=np.full(S, eta), delta=np.zeros(S), tau0=np.zeros(S), kappa=np.full(S, kappa),
                ab=np.full(S, 10.0), r=np.full(S, 0.8), gamma=np.full(S, gamma), sigma_arm=np.zeros(S))


def fake_family():
    return dict(rows=[dict(arm_id="R"), dict(arm_id="A"), dict(arm_id="B")], ref=0,
                y=np.array([10, 18, 25]), N=np.array([50, 50, 50]), dD=np.array([0.0, 2.0, -4.0]),
                dQ=np.array([0.0, 30.0, 10.0]), dTe=np.array([0.0, -10.0, -20.0]), chemo=np.array([0.0, 0.0, 1.0]),
                group="G1")


def test_predictive_density_is_exact_without_random_effects_and_conditional_drops_the_reference():
    d, fd = fake_draws(), fake_family()
    lin = -0.5 + 0.1 * (fd["dD"] + fd["dQ"] / 10 - 0.8 * fd["dTe"]) + 0.4 * fd["chemo"]
    lp = binom.logpmf(fd["y"], fd["N"], 1 / (1 + np.exp(-lin)))
    all_, ref = s6.predictive(d, fd, 5, np.random.default_rng(0))
    assert all_ == pytest.approx(lp.sum(), abs=1e-9)
    assert all_ - ref == pytest.approx(lp[1:].sum(), abs=1e-9)


def test_reference_arm_informs_the_baseline_in_the_conditional_task():
    # With a wide baseline, the conditional density of the other arms must exceed their marginal.
    d, fd = fake_draws(S=400), fake_family()
    d["tau0"] = np.full(400, 2.0)
    all_, ref = s6.predictive(d, fd, 200, np.random.default_rng(1))
    d1 = dict(fd, y=fd["y"][1:], N=fd["N"][1:], dD=fd["dD"][1:], dQ=fd["dQ"][1:], dTe=fd["dTe"][1:],
              chemo=fd["chemo"][1:], rows=fd["rows"][1:], ref=0)
    marg_other = s6.predictive(d, d1, 200, np.random.default_rng(2))[0]
    assert all_ - ref > marg_other


def test_bayesian_bootstrap_of_a_constant_difference_is_that_constant():
    b = s6.bayes_boot([0.7] * 10, 1)
    assert b["lo"] == pytest.approx(0.7) and b["hi"] == pytest.approx(0.7) and b["p_positive"] == 1
    b = s6.bayes_boot([1.0, -3.0, 0.5], 1)
    assert b["lo"] < b["mean"] < b["hi"] and b["mean"] == pytest.approx(-0.5)


def test_training_sets_never_contain_the_held_out_family():
    rows, fams = s6.s5.load_rows()
    for f in fams:
        tr, tf = s6.training_set(rows, fams, f)
        assert f not in tf and all(x["family"] != f for x in tr)
        assert len(tr) + sum(x["family"] == f for x in rows) == len(rows)


def test_scale_group_offset_for_held_out_family():
    rng = np.random.default_rng(0)
    fit = dict(kappa=np.zeros(3000), delta_all=np.tile([[0.5]], (3000, 1)))
    assert (s6.group_delta(fit, "G1", ["G2"], rng) == 0).all()
    assert (s6.group_delta(fit, "G2", ["G2"], rng) == 0.5).all()
    prior = s6.group_delta(fit, "other", ["G2"], rng)
    assert abs(prior.std() - 1) < 0.05


def test_hypno_benchmark_uses_its_randomized_control():
    fd = s6.external_family("HYPNO")
    ids = [x["arm_id"] for x in fd["rows"]]
    assert ids[fd["ref"]] == "HYPNO_NFX" and sorted(fd["y"]) == [190, 208]
    hfx = ids.index("HYPNO_HFX")
    assert fd["dD"][hfx] == pytest.approx(-11.0) and fd["chemo"][hfx] == 0
    assert "HYPNO" not in s6.s5.REF_ARM                      # temporary reference removed again


# Conditional predictive check (review F2)
def test_conditional_prediction_responds_to_the_reference_count():
    d, fd = fake_draws(S=400), fake_family()
    d["tau0"] = np.full(400, 1.5)
    lo = s6.contrast_check(dict(fd, y=np.array([2, 18, 25])), d, 1, 50)
    hi = s6.contrast_check(dict(fd, y=np.array([40, 18, 25])), d, 1, 50)
    # the other arm's predicted count moves only through the conditioned baseline
    assert hi[0]["predicted_count_median"] > lo[0]["predicted_count_median"] + 5


def test_uninformative_reference_gives_the_unconditional_prediction():
    d, fd = fake_draws(S=300), fake_family()
    d["tau0"] = np.full(300, 1.0)
    none = dict(fd, y=np.array([0, 18, 25]), N=np.array([0, 50, 50]))
    c = s6.contrast_check(none, d, 3, 20)
    assert c[0]["weight_ess"] == pytest.approx(300 * 20)            # all weights equal


def test_conditional_check_handles_zero_and_all_events_and_shares_the_reference():
    d, fd = fake_draws(S=200), fake_family()
    c = s6.contrast_check(dict(fd, y=np.array([0, 50, 0])), d, 4, 20)
    assert all(np.isfinite([x["observed_lor"], x["predicted_obs_lor_5"], x["predicted_obs_lor_95"]]).all() for x in c)
    assert len({x["weight_ess"] for x in c}) == 1 and len(c) == 2      # two arms, one shared reference


def test_probability_mixture_matches_a_pure_code_and_rejects_bad_weights():
    d, fd = fake_draws(S=100), fake_family()
    a = s6.predictive(d, s6.with_chemo(fd, "B", code=-1.0), 10, np.random.default_rng(5))
    b = s6.predictive(d, s6.with_chemo(fd, "B", mix=[(1.0, -1.0)]), 10, np.random.default_rng(5))
    assert a == pytest.approx(b)
    with pytest.raises(ValueError, match="sum to 1"):
        s6.predictive(d, s6.with_chemo(fd, "B", mix=[(0.5, -1.0)]), 10, np.random.default_rng(5))


def test_hypno_is_scored_under_three_declared_chemotherapy_scenarios():
    sc = s6.hypno_scenarios(s6.external_family("HYPNO"))
    assert set(sc) == {"code_0", "code_minus1", "mixture_75pct"}
    hfx = [x["arm_id"] for x in sc["code_0"]["rows"]].index("HYPNO_HFX")
    assert sc["code_minus1"]["chemo"][hfx] == -1 and sc["mixture_75pct"]["chemo_mix"][hfx] == [(0.75, -1.0), (0.25, 0.0)]


# Cache identity and eligibility (review F4)
def test_cache_without_identity_or_with_another_identity_is_refused(tmp_path):
    rows, fams = s6.s5.load_rows()
    tr, tf = s6.training_set(rows, fams, "Mani")
    np.savez(tmp_path / "x__legacy.npz", families=np.array(tf), gate=np.array([1, 0, 1.0, 900, 900]))
    with pytest.raises(RuntimeError, match="no run identity"):
        s6.fit_cached(None, "x__legacy", tr, tf, s6.MODELS["candidate"], str(tmp_path))
    np.savez(tmp_path / "x__other.npz", families=np.array(tf), identity=np.array("not-this-one"))
    with pytest.raises(RuntimeError, match="identity differs"):
        s6.fit_cached(None, "x__other", tr, tf, s6.MODELS["candidate"], str(tmp_path))
    (tmp_path / "x__broken.npz").write_bytes(b"PK\x03\x04 truncated")
    with pytest.raises(RuntimeError, match="unreadable cache"):
        s6.fit_cached(None, "x__broken", tr, tf, s6.MODELS["candidate"], str(tmp_path))
    with pytest.raises(RuntimeError, match="rescore-only mode does not sample"):
        s6.fit_cached(None, "x__missing", tr, tf, s6.MODELS["candidate"], str(tmp_path), rescore_only=True)


def test_subsets_failed_gates_and_legacy_gates_are_never_eligible():
    good = [dict(gate_pass=True, scope="complete schema gate")] * 3
    assert s6.numerical_rule(True, True, good)
    assert not s6.numerical_rule(True, False, good)                                  # subset run
    assert not s6.numerical_rule(True, True, good + [dict(gate_pass=False, scope="complete schema gate")])
    assert not s6.numerical_rule(True, True, good + [dict(gate_pass=True, scope="legacy narrow gate")])
    assert not s6.numerical_rule(False, True, good)
    assert not s6.numerical_rule(True, True, [])
