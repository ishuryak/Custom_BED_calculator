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
