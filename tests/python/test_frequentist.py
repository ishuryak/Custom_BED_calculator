"""Likelihood-only comparison engine (plan revision 3, review of 2026-11-12 F01-F06).
Deterministic tests only: no optimizer randomness, no sampling. Each test names the defect it
catches; each was observed red against a deliberately broken engine (record in
docs/development/2026_11_12_frequentist_review_response.md)."""
import copy
import gzip
import json
import math
import os
import sys

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import minimize, minimize_scalar
from scipy.special import gammaln, logsumexp

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "fit"))
from frequentist import designs as D  # noqa: E402
from frequentist import likelihood as L  # noqa: E402

TD = D.tumor_design()
ED = D.early_design()


# --- independent oracle: mode-scaled adaptive quad (scipy QUADPACK), written apart from the engine ---
def oracle(logf):
    h = lambda z: float(logf(z)) - 0.5 * z * z - 0.5 * math.log(2 * math.pi)
    m = minimize_scalar(lambda z: -h(z), bounds=(-30, 30), method="bounded", options=dict(xatol=1e-13))
    mx = h(m.x)
    val, _ = quad(lambda z: math.exp(h(z) - mx), -np.inf, np.inf, epsabs=1e-14, epsrel=1e-13, limit=500, points=None)
    return mx + math.log(val)


def grid_family(y, N, base, tau, sa, half=12.0, m=3001):
    """Independent nested reference: dense trapezoid over u and z on [-12, 12] (exponentially
    convergent for these Gaussian-tailed analytic integrands). Nested scipy.quad calls were tried
    first and overflowed at extreme outer values, so they are not used as the reference."""
    z = np.linspace(-half, half, m); dz = z[1] - z[0]
    lphi = -0.5 * z * z - 0.5 * math.log(2 * math.pi)
    tot = np.zeros(m)
    for i in range(len(y)):
        l = lbinom(float(y[i]), float(N[i]), base[i] + tau * z[:, None] + sa * z[None, :], "logit")
        tot += l[:, 0] if sa == 0 else logsumexp(l + lphi[None, :], axis=1) + math.log(dz)
    return float(logsumexp(tot + lphi) + math.log(dz))


def lbinom(y, n, eta, link):
    eta = np.asarray(eta, float)
    if link == "logit":
        lp, lq = -np.logaddexp(0, -eta), -np.logaddexp(0, eta)
    else:                                   # P(control) = exp(-e^eta)
        t = np.exp(np.minimum(eta, 700)); lp = -t
        lq = np.where(eta < -36, eta, np.log(-np.expm1(-t)))
    return gammaln(n + 1) - gammaln(y + 1) - gammaln(n - y + 1) + y * lp + (n - y) * lq


def one_arm_tumor(n, y, eta0, sigma, nodes=40):
    d = dict(TD, J=1, unit=np.array([0]), dsum=np.array([0.0]), d2sum=np.array([0.0]), T=np.array([0.0]),
             npat=np.array([n]), y=np.array([y]))
    th = dict(c=np.array([eta0 - L.smax(0.0) * 0]), alpha=1.0, ab=10.0, lam=0.0, Tk=0.0, g=0.0, sigma=sigma)
    return L.tumor_loglik(th, d, nodes)


# --- T1: identity layers (review F02) ---------------------------------------------------------------
def test_T1_designs_match_the_bayesian_fits():
    ids = D.reproduce_bayesian_identities()
    assert ids["tumor_di"].startswith("2c828af6") and ids["early"].startswith("4747b297")
    assert (TD["N"], TD["J"], int(TD["npat"].sum())) == (45, 21, 8279)
    assert (ED["I"], ED["J"], int(ED["N"].sum()), int(ED["y"].sum())) == (51, 22, 6042, 2945)


@pytest.mark.parametrize("change", ["count", "order", "time", "setting"])
def test_T1_file_design_and_run_identities_are_distinct(change):
    from frequentist import fit_ml as FM
    d = copy.deepcopy(TD)
    base = (d["design_sha256"], FM.run_identity("tumor", d, FM.SUPPORTS["matched"]))
    if change == "count":
        d["y"] = d["y"].copy(); d["y"][3] += 1
    elif change == "order":
        perm = np.r_[1, 0, np.arange(2, d["N"])]
        for k in ("unit", "dsum", "d2sum", "T", "npat", "y"):
            d[k] = d[k][perm]
        d["arms"] = [d["arms"][i] for i in perm]
    elif change == "time":
        d["T"] = d["T"].copy(); d["T"][0] -= 1        # inclusive -> elapsed convention on one arm
    new_design = D.canonical_digest(d, D.TUMOR_KEYS, dict(units=d["units"], arms=d["arms"]))
    sup = dict(FM.SUPPORTS["matched"], Tk=(0.0, 61.0)) if change == "setting" else FM.SUPPORTS["matched"]
    d["design_sha256"] = new_design
    new_run = FM.run_identity("tumor", d, sup)
    assert new_run != base[1]
    assert (new_design != base[0]) == (change != "setting")      # a setting moves the run, not the design


# --- T2: known value (review F02: sign of the binomial constant) -----------------------------------
def test_T2_full_nll_reproduces_the_archived_fit():
    m = json.load(open(os.path.join(ROOT, "docs", "refit_explainer", "main_free.json")))["di"]
    lc = float(L.log_choose(TD["npat"], TD["y"]).sum())
    assert lc == pytest.approx(5211.115489844190, abs=1e-6)
    th = dict(alpha=m["alpha"], ab=10.0, lam=m["lam"], Tk=m["third"], g=m["g"], sigma=0.0, c=np.zeros(TD["J"]))

    def nll(c):
        th["c"] = c
        return -L.tumor_loglik(th, TD)
    r = minimize(nll, np.zeros(TD["J"]), method="BFGS", options=dict(gtol=1e-10))
    full_nll = m["neg_log_lik"] - lc                    # archived value omits log C(N, y)
    assert full_nll == pytest.approx(131.291127555958, abs=1e-9)
    assert r.fun == pytest.approx(full_nll, abs=1e-3)
    assert abs(r.fun - (m["neg_log_lik"] + lc)) > 1e4   # the revision-2 instruction ("add") is refused


# --- T3: integration accuracy (review F01) ---------------------------------------------------------
TUMOR_BIG = int(TD["npat"].max())
_iy = int(np.argmax(TD["npat"]))


@pytest.mark.parametrize("n,y,tcp0,sigma", [(750, 525, 0.7, 0.25), (750, 525, 0.7, 1.0),
                                           (TUMOR_BIG, int(TD["y"][_iy]), None, 0.05),
                                           (TUMOR_BIG, int(TD["y"][_iy]), None, 1.0), (40, 39, 0.9, 2.0)])
def test_T3_tumor_arm_integral_matches_an_independent_integrator(n, y, tcp0, sigma):
    eta0 = math.log(-math.log(tcp0 if tcp0 else y / n))
    exact = oracle(lambda z: lbinom(y, n - 0, eta0 + sigma * z, "cloglog"))
    v20, v40, v80 = (one_arm_tumor(n, y, eta0, sigma, k) for k in (20, 40, 80))
    assert abs(v20 - v40) < 1e-6                # the ladder accepts 40 here ...
    assert v40 == pytest.approx(exact, abs=1e-9) and v80 == pytest.approx(exact, abs=1e-9)   # ... and 40 is exact


def test_T3_fixed_gauss_hermite_is_detected_as_inaccurate():
    """Positive control for the test above: the ordinary 30-node rule the revision-2 plan proposed is
    off by about -6.87 at the review's fixture, so an engine using it would fail T3."""
    eta0, s = math.log(-math.log(0.7)), 1.0
    x, w = np.polynomial.hermite.hermgauss(30)
    naive = float(logsumexp(np.log(w) + lbinom(525, 750, eta0 + s * math.sqrt(2) * x, "cloglog")) - 0.5 * math.log(math.pi))
    exact = oracle(lambda z: lbinom(525, 750, eta0 + s * z, "cloglog"))
    assert naive - exact == pytest.approx(-6.87156, abs=1e-3)


def _family(fid):
    idx = [i for i, f in enumerate(ED["fam"]) if ED["families"][f] == fid]
    d = dict(ED, I=len(idx), J=1, G=0, fam=np.zeros(len(idx), int), N=ED["N"][idx], y=ED["y"][idx],
             dD=np.zeros(len(idx)), dQ=np.zeros(len(idx)), dTe=np.zeros(len(idx)), chemo=np.zeros(len(idx)),
             grp=np.zeros((1, 0)))
    return d


@pytest.mark.parametrize("tau,sa", [(0.5, 0.0), (1.0, 0.0), (2.0, 0.0), (1.0, 0.5), (0.3, 1.5)])
def test_T3_nested_family_integral_matches_an_independent_integrator(tau, sa):
    d = _family("RTOG0129")
    mu = np.log(d["y"] / (d["N"] - d["y"]))
    assert list(d["N"]) == [361, 360] and list(d["y"]) == [141, 119]
    # engine: arm offsets enter through m0 = 0 and a chemo-coded offset per arm
    th = dict(m0=0.0, delta=np.zeros(0), tau=tau, kappa=0.0, ab=10.0, r=1.0, gamma=1.0, sigma_arm=sa)
    d["chemo"] = mu
    got = {n: L.early_loglik(th, d, (n, n)) for n in (20, 40)}

    exact = grid_family(d["y"], d["N"], mu, tau, sa)
    for n in (20, 40):
        assert got[n] == pytest.approx(exact, abs=1e-9)
    if sa == 0 and tau == 1.0:                           # the review's -7.974 outer counterexample
        x, w = np.polynomial.hermite.hermgauss(20)
        naive = float(logsumexp(np.log(w) + sum(lbinom(d["y"][i], d["N"][i], mu[i] + tau * math.sqrt(2) * x, "logit")
                                                  for i in range(2))) - 0.5 * math.log(math.pi))
        assert naive - exact == pytest.approx(-7.97402, abs=1e-3)


def test_T3_zero_variance_branches_are_exact():
    th = dict(c=np.zeros(TD["J"]), alpha=0.044, ab=10.6, lam=0.028, Tk=37.2, g=0.0117, sigma=0.0)
    eta = L.tumor_eta0(th, TD)
    exact = float(lbinom(TD["y"], TD["npat"], eta, "cloglog").sum())
    assert L.tumor_loglik(th, TD) == pytest.approx(exact, abs=1e-10)
    the = dict(m0=-0.5, delta=np.zeros(ED["G"]), tau=0.0, kappa=0.05, ab=12.9, r=1.55, gamma=0.55, sigma_arm=0.0)
    exact = float(lbinom(ED["y"], ED["N"], L.early_base(the, ED), "logit").sum())
    assert L.early_loglik(the, ED) == pytest.approx(exact, abs=1e-10)


@pytest.mark.parametrize("scale", [1.0, 4.0])
def test_T3_total_likelihoods_converge_on_the_ladder(scale):
    th = dict(c=np.full(TD["J"], -0.5), alpha=0.044, ab=10.6, lam=0.028, Tk=37.2, g=0.0117, sigma=0.065 * scale * 4)
    v20, v80 = L.tumor_loglik(th, TD, 20), L.tumor_loglik(th, TD, 80)
    assert abs(v20 - v80) < 1e-6
    the = dict(m0=-0.5, delta=np.zeros(ED["G"]), tau=0.8 * scale, kappa=0.047, ab=12.9, r=1.55, gamma=0.55,
               sigma_arm=0.24 * scale)
    val, nodes, change = L.converged_loglik(L.early_loglik, the, ED, (10, 10))
    assert change < 1e-6 and abs(val - L.early_loglik(the, ED, (80, 80))) < 1e-6


# --- T4: contrasts (review F02: same-row identity, rounded-input tolerance) -------------------------
def test_T4_tumor_contrasts_reproduce_the_posterior_summary_as_a_sanity_check():
    import csv
    rows = list(csv.DictReader(open(os.path.join(ROOT, "draws", "tumor_di_primary.csv"))))
    assert len(rows) == 4000 and len({r["di_draw"] for r in rows}) == 1000        # 1000 draws, each 4 times
    vals = np.array([L.tumor_contrasts({k: float(r[k]) for k in ("alpha", "ab", "Tk", "g")} | {"lam": float(r["lambda"])}, TD)
                     for r in rows])
    summ = json.load(open(os.path.join(ROOT, "fit", "output_refit_2026_10_25c_main", "di_primary", "summary.json")))
    med = {c["regimen"]: c["median"] for c in summ["contrasts_vs_anchor"]}
    diffs = [abs(np.median(vals[:, k]) - med[name]) for k, name in enumerate(TD["regimens"])]
    assert len(diffs) == 10 and max(diffs) < 0.1                 # review measured 0.09225 Gy


def test_T4_early_contrasts_match_the_stored_rounded_generated_quantities():
    import pandas as pd
    dr = pd.read_csv(os.path.join(ROOT, "docs", "development", "2026_11_10_early_toxicity", "step5_review",
                                  "overdispersed", "U02", "draws.csv.gz"))
    th = dict(kappa=dr["kappa"].to_numpy()[:, None], ab=dr["ab"].to_numpy()[:, None], r=dr["r"].to_numpy()[:, None])
    pred = L.early_contrasts(th, ED)
    stored = dr[[f"lor_dec[{k}]" for k in range(1, 5)]].to_numpy()
    err = np.abs(pred - stored)
    assert err.max() < 2e-5          # six-significant-digit inputs and outputs (review: max 1.51e-5)
    assert err.max() > 1e-6          # documents why the revision-2 1e-6 tolerance was unattainable


def test_T4_contrast_formula_on_unrounded_values():
    th = dict(alpha=0.05, ab=8.0, lam=0.03, Tk=30.0, g=0.01)
    n, dpf, T = 20, 2.75, 26.0
    K = n * dpf * (1 + dpf / th["ab"]); K0 = 35 * 2.0 * (1 + 2.0 / th["ab"])
    want = (K - (th["lam"] * L.smax(T - 30) + 0.01 * T) / 0.05) - (K0 - (th["lam"] * L.smax(47 - 30) + 0.01 * 47) / 0.05)
    assert L.tumor_contrasts(th, TD)[0] == pytest.approx(want, abs=1e-12)
    te = dict(kappa=0.05, ab=12.0, r=1.5)
    assert L.early_contrasts(te, ED)[1] == pytest.approx(0.05 * (ED["kD"][1] + ED["kQ"][1] / 12 - 1.5 * ED["kTe"][1]), abs=1e-14)


def test_T3_inner_mode_search_does_not_cycle():
    """Regression: the first engine draft used a clipped Newton that cycled between z = 0 and z = 5 on
    a real arm (398 patients, 269 events, eta0 = -2, spread 0.96) and returned the integral 114 log
    units too low. The bracketed search finds the mode, and the integral matches the reference."""
    y, n, e0, s = 269.0, 398.0, -2.0, 0.96
    z, _ = L._mode(np.array([e0]), s, L.logit_terms, np.array([y]), np.array([n]))
    _, l1, _ = L.logit_terms(e0 + s * z, y, n)
    assert abs(s * l1[0] - z[0]) < 1e-8 and z[0] > 0.5
    v, _, _ = L.agq_1d(np.array([e0]), s, L.logit_terms, (np.array([y]), np.array([n])), 20)
    exact = oracle(lambda t: lbinom(y, n, e0 + s * t, "logit")) - float(L.log_choose(n, y))
    assert v[0] == pytest.approx(exact, abs=1e-9)


def test_T3_full_early_likelihood_matches_a_grid_reference_at_large_spreads():
    th = dict(m0=-0.5, delta=np.zeros(ED["G"]), tau=3.2, kappa=0.047, ab=12.9, r=1.55, gamma=0.55, sigma_arm=0.96)
    base = L.early_base(th, ED)
    fam_eng = L.early_loglik(th, ED, (20, 20), per_family=True)
    for j in range(ED["J"]):
        idx = np.flatnonzero(ED["fam"] == j)
        ref = grid_family(ED["y"][idx], ED["N"][idx], base[idx], th["tau"], th["sigma_arm"])   # includes log C
        assert fam_eng[j] == pytest.approx(ref, abs=1e-8), ED["families"][j]


# --- T5: analytic gradients against central differences (deterministic code check, review F03) -----
def _fd(f, th, k, h=1e-6, j=None):
    a, b = copy.deepcopy(th), copy.deepcopy(th)
    if j is None:
        hh = h * max(1, abs(th[k])); a[k] += hh; b[k] -= hh
    else:
        hh = h; a[k][j] += hh; b[k][j] -= hh
    return (f(a) - f(b)) / (2 * hh)


def test_T5_tumor_gradient_matches_finite_differences():
    th = dict(c=np.linspace(-1, 0.3, TD["J"]), alpha=0.044, ab=10.6, lam=0.028, Tk=37.2, g=0.0117, sigma=0.15)
    _, g = L.tumor_loglik_grad(th, TD, 40)
    f = lambda t: L.tumor_loglik(t, TD, 40)
    for k in L.TUMOR_SCALARS:
        assert g[k] == pytest.approx(_fd(f, th, k), rel=1e-5, abs=1e-5), k
    assert g["c"][5] == pytest.approx(_fd(f, th, "c", j=5), rel=1e-5, abs=1e-6)


def test_T5_early_gradient_matches_finite_differences():
    th = dict(m0=-0.5, delta=np.array([0.3]), tau=0.8, kappa=0.047, ab=12.9, r=1.55, gamma=0.55, sigma_arm=0.24)
    _, g = L.early_loglik_grad(th, ED, (30, 30))
    f = lambda t: L.early_loglik(t, ED, (30, 30))
    for k in ("m0", "tau", "kappa", "ab", "r", "gamma", "sigma_arm"):
        assert g[k] == pytest.approx(_fd(f, th, k), rel=1e-5, abs=1e-5), k
    assert g["delta"][0] == pytest.approx(_fd(f, th, "delta", j=0), rel=1e-5)


# --- T8: inference controls on likelihoods with known answers (review F05) --------------------------
from frequentist import fit_ml as FM  # noqa: E402


class Toy:
    """Two parameters named like real targets (alpha, gamma) with a chosen log-likelihood."""
    kind = "toy"
    nodes = 20

    def __init__(self, ll, lo=(-10, -10), hi=(10, 10), edge=(("comp", "comp"), ("comp", "comp"))):
        self.ll, self.names = ll, ["alpha", "gamma"]
        self.lo, self.hi, self.edge, self.sc = np.array(lo, float), np.array(hi, float), list(edge), np.ones(2)
        self.targets = dict(alpha=lambda th, d: th["alpha"], gamma=lambda th, d: th["gamma"],
                            lam_plus_g=lambda th, d: th["alpha"] + th["gamma"])
        self.d = None

    def unpack(self, x): return dict(alpha=float(x[0]), gamma=float(x[1]))
    def pack(self, th): return np.array([th["alpha"], th["gamma"]])
    def loglik(self, th, nodes=None): return float(self.ll(np.array([th["alpha"], th["gamma"]])))

    def loglik_grad(self, x, nodes=None):
        h = 1e-6
        g = np.array([(self.ll(x + h * e) - self.ll(x - h * e)) / (2 * h) for e in np.eye(2)])
        return float(self.ll(x)), g

    def target(self, name, th): return float(self.targets[name](th, None))
    def all_targets(self, th): return {k: self.target(k, th) for k in self.targets}


def _toy_profile(toy, name, x_hat):
    mle = dict(x=list(x_hat), ll=toy.ll(np.array(x_hat)), nodes=20)
    return FM.profile(toy, name, mle, np.random.default_rng(1))


Z = 1.959963984540054


def test_T8_quadratic_profiles_recover_known_intervals():
    toy = Toy(lambda x: -0.5 * ((x[0] - 1) / 0.2) ** 2 - 0.5 * ((x[1] - 3) / 0.5) ** 2)
    p = _toy_profile(toy, "alpha", [1, 3])
    (piece,) = p["confidence_set"]
    assert piece["lower"]["closed"] and piece["upper"]["closed"]
    assert piece["lower"]["value"] == pytest.approx(1 - Z * 0.2, abs=1e-5)
    assert piece["upper"]["value"] == pytest.approx(1 + Z * 0.2, abs=1e-5)
    p = _toy_profile(toy, "lam_plus_g", [1, 3])                      # derived target, equality-constrained
    (piece,) = p["confidence_set"]
    sd = math.sqrt(0.2 ** 2 + 0.5 ** 2)
    assert piece["lower"]["value"] == pytest.approx(4 - Z * sd, abs=1e-4)
    assert piece["upper"]["value"] == pytest.approx(4 + Z * sd, abs=1e-4)


def test_T8_an_estimate_at_a_support_edge_gives_an_open_end():
    toy = Toy(lambda x: -0.5 * ((x[0] - 1) / 0.2) ** 2 - 0.5 * (x[1] - 3) ** 2, lo=(0.8, -10),
              edge=(("support", "comp"), ("comp", "comp")))
    (piece,) = _toy_profile(toy, "alpha", [1, 3])["confidence_set"]
    assert not piece["lower"]["closed"] and "support limit" in piece["lower"]["why"]
    assert piece["lower"]["value"] == pytest.approx(0.8)
    assert piece["upper"]["closed"]


def test_T8_a_flat_profile_is_reported_flat_and_open_at_computational_limits():
    toy = Toy(lambda x: -0.5 * ((x[0] - 1) / 0.2) ** 2, lo=(-10, -5), hi=(10, 5))
    p = _toy_profile(toy, "gamma", [1, 0.0])
    assert p["flat"]
    (piece,) = p["confidence_set"]
    assert "computational limit" in piece["lower"]["why"] and "computational limit" in piece["upper"]["why"]


def test_T8_a_bimodal_profile_gives_a_disconnected_set():
    toy = Toy(lambda x: float(logsumexp([-0.5 * ((x[0] - 2) / 0.3) ** 2, -0.5 * ((x[0] + 2) / 0.3) ** 2]))
              - 0.5 * x[1] ** 2, lo=(-5, -5), hi=(5, 5))
    p = _toy_profile(toy, "alpha", [2, 0])
    assert p["disconnected"] and len(p["confidence_set"]) == 2
    lows = sorted(c["lower"]["value"] for c in p["confidence_set"])
    assert lows[0] < -2 < lows[1]


# --- T8: restart, corruption and publication controls (review F06) ----------------------------------
@pytest.fixture(scope="module")
def tumor_mle():
    m = FM.Model("tumor", "matched")
    r = FM.multistart(m, np.random.default_rng([FM.SEED, 1]), 3)
    r["identity"] = FM.run_identity("tumor", m.d, "matched")
    return r


def _mtimes(d):
    return {f: os.stat(os.path.join(d, f)).st_mtime_ns for f in os.listdir(d) if f.startswith("rep_")}


def test_T8_resume_reruns_nothing_and_reruns_damaged_records_identically(tmp_path, tumor_mle):
    bdir = str(tmp_path / "boot")
    first = FM.bootstrap("tumor", "matched", tumor_mle, 3, bdir, workers=3)
    assert [r["status"] for r in first] == ["ok"] * 3
    before = _mtimes(bdir)
    FM.bootstrap("tumor", "matched", tumor_mle, 3, bdir, workers=3)          # resume: nothing to do
    assert _mtimes(bdir) == before
    p2, p3 = os.path.join(bdir, "rep_00002.json"), os.path.join(bdir, "rep_00003.json")
    txt = open(p2).read(); open(p2, "w").write(txt.replace('"attempts": 1', '"attempts": 1 '))   # same JSON, still valid
    rec2 = json.loads(txt); rec2["fit"]["ll"] += 1.0                          # an edited value: digest fails
    open(p2, "w").write(json.dumps(rec2))
    open(p3, "w").write(txt[: len(txt) // 2])                                 # truncated
    again = FM.bootstrap("tumor", "matched", tumor_mle, 3, bdir, workers=3)
    for b in (1, 2):
        assert again[b]["y_sha256"] == first[b]["y_sha256"] and again[b]["seed"] == first[b]["seed"]
        assert again[b]["fit"]["ll"] == pytest.approx(first[b]["fit"]["ll"], abs=1e-6)
    assert os.stat(os.path.join(bdir, "rep_00001.json")).st_mtime_ns == before["rep_00001.json"]


def test_T8_a_stale_identity_is_not_reused(tmp_path, tumor_mle):
    bdir = str(tmp_path / "boot")
    FM.bootstrap("tumor", "matched", tumor_mle, 1, bdir, workers=1)
    rec = json.load(open(os.path.join(bdir, "rep_00001.json")))
    assert FM.read_record(os.path.join(bdir, "rep_00001.json"), tumor_mle["identity"]) is not None
    assert FM.read_record(os.path.join(bdir, "rep_00001.json"), "0" * 64) is None
    other = dict(tumor_mle, identity="f" * 64)
    out = FM.bootstrap("tumor", "matched", other, 1, bdir, workers=1)
    assert out[0]["identity"] == "f" * 64 and out[0]["y_sha256"] == rec["y_sha256"]


def test_T8_a_failed_record_is_retried_with_its_seed(tmp_path, tumor_mle):
    bdir = str(tmp_path / "boot"); os.makedirs(bdir)
    fake = dict(identity=tumor_mle["identity"], replicate=1, seed=FM.replicate_seed(1), attempts=1, status="failed",
                why="forced", generating="x", y_sha256="", seconds=0.0, fit=None, boundary_solution=False)
    fake["record_sha256"] = FM.record_digest(fake)
    json.dump(fake, open(os.path.join(bdir, "rep_00001.json"), "w"))
    (rec,) = FM.bootstrap("tumor", "matched", tumor_mle, 1, bdir, workers=1)
    assert rec["attempts"] == 2 and rec["seed"] == FM.replicate_seed(1) and rec["status"] == "ok"


def test_T8_incomplete_or_mixed_generations_are_refused(tmp_path, tumor_mle):
    bdir = str(tmp_path / "boot")
    recs = FM.bootstrap("tumor", "matched", tumor_mle, 2, bdir, workers=2)
    with pytest.raises(RuntimeError, match="no MANIFEST"):
        FM.load_generation(bdir)                                    # interrupted before publication
    FM.publish_bootstrap("tumor", "matched", tumor_mle, 2, bdir, 2, ["test"])
    man, got, s = FM.load_generation(bdir)
    assert len(got) == 2 and s["ok"] == 2
    r = json.load(open(os.path.join(bdir, "rep_00002.json")))
    r["identity"] = "e" * 64; r["record_sha256"] = FM.record_digest(r)       # a record from another generation
    json.dump(r, open(os.path.join(bdir, "rep_00002.json"), "w"))
    with pytest.raises(RuntimeError, match="generation refused"):
        FM.load_generation(bdir)


# --- T6: independent R cross-checks (required: an unavailable R fails, review F01/F06) --------------
def test_T6_r_cross_checks_pass():
    from frequentist import r_check as RC
    if not os.path.exists(RC.RSCRIPT) and os.environ.get("FREQ_ALLOW_NO_R") == "1":
        pytest.skip("R unavailable and FREQ_ALLOW_NO_R=1: the required cross-check has NOT passed")
    out = os.path.join(ROOT, "fit", "output_frequentist", "_t6_test")   # Windows R cannot use a WSL /tmp folder
    RC.main(["--out", out])                              # sys.exit(1) on any failed comparison
    res = json.load(open(os.path.join(out, "r_check.json")))
    assert res["passed"] and res["tumor"]["objective"]["ok"] and res["early_glm"]["kappa_profile"]["ok"]


# --- T7: golden values (objective and identifiable predictions, not ridge coordinates) --------------
GOLDEN = dict(
    versions=dict(numpy="2.4.3", scipy="1.17.1"),
    tumor_ll=-131.09788481544865,
    tumor_contrasts=[-3.040579625, -10.883852583, -5.169377070, 0.234248024, 2.326401972, -2.519630101,
                     -4.224767120, 4.392892379, 7.264355208, 8.466886621],
    early_ll=-179.77568720592626,
    early_lor=[1.393028203, 0.960939760, 0.640470443, 0.364575204])


def test_T7_golden_estimates():
    import scipy
    assert (np.__version__, scipy.__version__) == tuple(GOLDEN["versions"].values()), \
        "dependency changed: re-derive the golden values in the same commit as the upgrade"
    m = FM.Model("tumor", "matched")
    r = FM.multistart(m, np.random.default_rng([FM.SEED, 1]), 3)
    assert r["ll"] == pytest.approx(GOLDEN["tumor_ll"], abs=1e-6)
    assert [r["targets"][f"contrast[{k}]"] for k in range(10)] == pytest.approx(GOLDEN["tumor_contrasts"], abs=1e-4)
    e = FM.Model("early", "early_matched")
    r = FM.multistart(e, np.random.default_rng([FM.SEED, 1]), 2)
    assert r["ll"] == pytest.approx(GOLDEN["early_ll"], abs=1e-6)
    assert [r["targets"][f"lor[{k}]"] for k in range(4)] == pytest.approx(GOLDEN["early_lor"], abs=1e-4)
