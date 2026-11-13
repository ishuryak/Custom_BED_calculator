"""Likelihood-only fits, profile likelihoods and the parametric bootstrap for the comparison of plan
revision 3 (docs/development/2026_11_12_frequentist_comparison_plan.md, sections 4-6).

Usage (from the repository root; nothing here samples with Stan):
    python3 fit/frequentist/fit_ml.py validate                      identities, accuracy checks; no fitting
    python3 fit/frequentist/fit_ml.py mle      --model tumor|early --support matched|expanded
    python3 fit/frequentist/fit_ml.py profiles --model ... --support ... [--targets a,b]
    python3 fit/frequentist/fit_ml.py bootstrap --model ... --B 20 [--out DIR]    (pilot: B = 20)
Env: FREQ_WORKERS (default cores - 2), FREQ_OUT (default fit/output_frequentist).

Every result is written atomically (temporary file, then rename) into one generation folder named by
the run identity. A bootstrap replicate is one JSON record holding its seed, status, attempts and a
sha256 of its own content; a resumed run accepts only records of the same identity whose digest
verifies, retries failed ones under the declared policy, and never reruns a completed one. Assembled
summaries and MANIFEST.json are written only after every required record verifies.
"""
import argparse
import concurrent.futures as cf
import hashlib
import json
import math
import os
import platform
import sys
import time
import zlib

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import scipy
from scipy.optimize import brentq, minimize
from scipy.special import expit

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from frequentist import designs as D  # noqa: E402
from frequentist import likelihood as L  # noqa: E402

ROOT = D.ROOT
OUT = os.environ.get("FREQ_OUT", os.path.join(ROOT, "fit", "output_frequentist"))
SEED = 20261112
CUTOFF = 0.5 * 3.841458820694124          # regular one-parameter 95% profile cutoff (chi-square_1 / 2)
QUAD_TOL = 1e-6                            # |change in total log-likelihood| when nodes double
CONTRAST_TOL = 1e-4                        # natural units
RESIDUAL_TOL = 1e-7                        # derived-target equality constraint
FAIL_LIMIT = 0.05                          # unresolved bootstrap failures that block a final interval
MULTI_TOL = 1e-3                           # starts within this log-likelihood of the best count as ties

# --- supports: (lower, upper, kind of each edge). "support" edges are the Bayesian model's own
# parameter support (or, for expanded g >= 0, a stated physical assumption); "comp" edges are
# computational limits, never evidence: an estimate or interval reaching one is reported as open.
SUPPORTS = {
    "matched": dict(c=(-30, 30, "comp", "comp"), alpha=(1e-6, 2.0, "support", "comp"), ab=(0.05, 1000.0, "comp", "comp"),
                    lam=(0.0, 2.0, "support", "comp"), Tk=(0.0, 60.0, "support", "support"),
                    g=(0.0046, 0.021, "support", "support"), sigma=(0.0, 5.0, "support", "comp")),
    "expanded": dict(c=(-30, 30, "comp", "comp"), alpha=(1e-6, 2.0, "support", "comp"), ab=(0.05, 1000.0, "comp", "comp"),
                     lam=(0.0, 2.0, "support", "comp"), Tk=(0.0, 120.0, "support", "comp"),
                     g=(0.0, 0.2, "support", "comp"), sigma=(0.0, 5.0, "support", "comp")),
    "early_matched": dict(m0=(-20, 20, "comp", "comp"), delta=(-20, 20, "comp", "comp"), tau=(0.0, 10.0, "support", "comp"),
                          kappa=(0.0, 1.0, "support", "comp"), ab=(0.05, 1000.0, "comp", "comp"),
                          r=(0.0, 20.0, "support", "comp"), gamma=(-20, 20, "comp", "comp"),
                          sigma_arm=(0.0, 10.0, "support", "comp")),
}
EDGE_WORD = dict(comp="computational", support="support")
MODEL_SUPPORTS = {"tumor": ("matched", "expanded"), "early": ("early_matched",)}


# --- models ------------------------------------------------------------------------------------------
class Model:
    def __init__(self, kind, support_name, d=None):
        self.kind, self.support_name = kind, support_name
        self.d = d if d is not None else (D.tumor_design() if kind == "tumor" else D.early_design())
        self.sup = SUPPORTS[support_name]
        if kind == "tumor":
            self.blocks = [("c", self.d["J"])] + [(k, 1) for k in L.TUMOR_SCALARS]
            self.scale = dict(c=1.0, alpha=0.05, ab=10.0, lam=0.03, Tk=30.0, g=0.01, sigma=0.1)
            self.nodes = 20
        else:
            self.blocks = [("m0", 1), ("delta", self.d["G"])] + [(k, 1) for k in ("tau", "kappa", "ab", "r", "gamma", "sigma_arm")]
            self.scale = dict(m0=1.0, delta=1.0, tau=1.0, kappa=0.05, ab=10.0, r=1.0, gamma=1.0, sigma_arm=0.3)
            self.nodes = (20, 20)
        self.names = [k if n == 1 and k not in ("c", "delta") else f"{k}[{i}]" for k, n in self.blocks for i in range(n)]
        self.sc = np.concatenate([np.full(n, self.scale[k]) for k, n in self.blocks])
        lo = np.concatenate([np.full(n, self.sup[k][0]) for k, n in self.blocks])
        hi = np.concatenate([np.full(n, self.sup[k][1]) for k, n in self.blocks])
        self.lo, self.hi = lo, hi
        self.edge = [(self.sup[k][2], self.sup[k][3]) for k, n in self.blocks for _ in range(n)]
        self.targets = TUMOR_TARGETS if kind == "tumor" else EARLY_TARGETS

    def unpack(self, x):
        th, i = {}, 0
        for k, n in self.blocks:
            th[k] = np.array(x[i:i + n], float) if k in ("c", "delta") else float(x[i])
            i += n
        return th

    def pack(self, th):
        return np.concatenate([np.atleast_1d(np.asarray(th[k], float)) for k, _ in self.blocks])

    def loglik(self, th, nodes=None):
        nodes = nodes or self.nodes
        return (L.tumor_loglik if self.kind == "tumor" else L.early_loglik)(th, self.d, nodes)

    def loglik_grad(self, x, nodes=None):
        th = self.unpack(x)
        ll, g = (L.tumor_loglik_grad if self.kind == "tumor" else L.early_loglik_grad)(th, self.d, nodes or self.nodes)
        return ll, self.pack(g)

    def target(self, name, th):
        return float(self.targets[name](th, self.d))

    def all_targets(self, th):
        return {k: self.target(k, th) for k in self.targets}


def _contrast(k):
    return lambda th, d: L.tumor_contrasts(th, d)[k]


def _lor(k):
    return lambda th, d: L.early_contrasts(th, d)[k]


TUMOR_TARGETS = dict(alpha=lambda th, d: th["alpha"], ab=lambda th, d: th["ab"], lam=lambda th, d: th["lam"],
                     Tk=lambda th, d: th["Tk"], g=lambda th, d: th["g"], sigma=lambda th, d: th["sigma"],
                     lam_plus_g=lambda th, d: th["lam"] + th["g"],
                     rate_bed=lambda th, d: (th["lam"] + th["g"]) / th["alpha"],
                     rate_eqd2=lambda th, d: (th["lam"] + th["g"]) / th["alpha"] / (1 + 2 / th["ab"]),
                     **{f"contrast[{k}]": _contrast(k) for k in range(10)})
EARLY_TARGETS = dict(kappa=lambda th, d: th["kappa"], ab=lambda th, d: th["ab"], r=lambda th, d: th["r"],
                     gamma=lambda th, d: th["gamma"], tau=lambda th, d: th["tau"], sigma_arm=lambda th, d: th["sigma_arm"],
                     **{f"lor[{k}]": _lor(k) for k in range(4)})
PARAMETER_TARGETS = {"alpha", "ab", "lam", "Tk", "g", "sigma", "kappa", "r", "gamma", "tau", "sigma_arm"}
# Units of every reported target.
UNITS = dict(alpha="per Gy", ab="Gy", lam="per day", Tk="days", g="per day", sigma="cloglog scale",
             lam_plus_g="per day", rate_bed="Gy of BED per day", rate_eqd2="Gy EQD2 per day",
             kappa="log odds per Gy of early BED", r="Gy of BED per day", gamma="log odds",
             tau="log odds", sigma_arm="log odds", **{f"contrast[{k}]": "Gy of tumor BED" for k in range(10)},
             **{f"lor[{k}]": "log odds ratio" for k in range(4)})


# --- identity ----------------------------------------------------------------------------------------
def code_digest():
    return {f: D.sha_file(os.path.join(HERE, f)) for f in ("likelihood.py", "fit_ml.py", "designs.py")}


def run_identity(kind, d, support, extra=None):
    payload = dict(model=kind, design_sha256=d["design_sha256"], file_sha256=d["file_sha256"],
                   support={k: list(v) for k, v in (support.items() if isinstance(support, dict) else SUPPORTS[support].items())},
                   quad=dict(start=list(np.atleast_1d(20)), tol=QUAD_TOL), cutoff=CUTOFF, seed=SEED,
                   optimizer=dict(method="L-BFGS-B, SLSQP for derived profiles", starts=START_POLICY, multi_tol=MULTI_TOL),
                   fail_limit=FAIL_LIMIT, code=code_digest(),
                   software=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                   **(extra or {}))
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


# --- optimization ------------------------------------------------------------------------------------
START_POLICY = dict(mle_random=15, boot_random=2, retry_random=8, retry_node_factor=2)


def deterministic_start(m):
    """Posterior medians of the deployed fits (the natural reference point), baselines from the data."""
    if m.kind == "tumor":
        d = m.d
        p = np.clip(np.bincount(d["unit"], d["y"]) / np.bincount(d["unit"], d["npat"]), 0.02, 0.98)
        th = dict(c=np.log(-np.log(p)), alpha=0.044, ab=10.6, lam=0.028, Tk=37.2, g=0.0117, sigma=0.065)
    else:
        th = dict(m0=-0.54, delta=np.zeros(m.d["G"]), tau=0.80, kappa=0.048, ab=12.9, r=1.55, gamma=0.55, sigma_arm=0.24)
    return np.clip(m.pack(th), m.lo, m.hi)


def random_start(m, rng):
    x = deterministic_start(m)
    th = m.unpack(x)
    lu = lambda a, b: math.exp(rng.uniform(math.log(a), math.log(b)))
    if m.kind == "tumor":
        th.update(c=th["c"] + rng.normal(0, 0.3, len(th["c"])), alpha=lu(0.005, 0.5), ab=lu(1, 100), lam=lu(0.001, 0.3),
                  Tk=rng.uniform(5, min(55, m.hi[m.names.index("Tk")])),
                  g=rng.uniform(m.sup["g"][0], min(m.sup["g"][1], 0.05)), sigma=rng.uniform(0, 0.5))
    else:
        th.update(m0=rng.normal(-0.5, 1), delta=rng.normal(0, 0.5, m.d["G"]), tau=rng.uniform(0.1, 2), kappa=lu(0.005, 0.3),
                  ab=lu(1, 100), r=lu(0.1, 5), gamma=rng.normal(0, 1), sigma_arm=rng.uniform(0.01, 1))
    return np.clip(m.pack(th), m.lo, m.hi)


def _objective(m, nodes, fixed=None):
    def f(xs):
        x = xs * m.sc
        try:
            ll, g = m.loglik_grad(x, nodes)
        except (ArithmeticError, FloatingPointError, ValueError):
            return 1e12, np.zeros_like(xs)
        if not np.isfinite(ll):
            return 1e12, np.zeros_like(xs)
        g = -g * m.sc
        if fixed is not None:
            g[fixed] = 0.0
        return -ll, g
    return f


def projected_gradient(m, x, g):
    """Gradient components that could still improve the objective given the bounds (KKT check)."""
    pg = g.copy()
    at_lo, at_hi = x <= m.lo + 1e-10 * (1 + abs(m.lo)), x >= m.hi - 1e-10 * (1 + abs(m.hi))
    pg[at_lo & (g < 0)] = 0.0         # g = d loglik / dx: at a lower bound, a negative slope is optimal
    pg[at_hi & (g > 0)] = 0.0
    return pg, at_lo, at_hi


def local_fit(m, x0, nodes=None, lo=None, hi=None):
    lo = m.lo if lo is None else lo
    hi = m.hi if hi is None else hi
    fixed = np.flatnonzero(lo == hi)
    r = minimize(_objective(m, nodes, fixed if len(fixed) else None), np.clip(x0, lo, hi) / m.sc, jac=True,
                 method="L-BFGS-B", bounds=list(zip(lo / m.sc, hi / m.sc)),
                 options=dict(maxiter=20000, maxfun=40000, ftol=1e-15, gtol=1e-9))
    x = np.clip(r.x * m.sc, lo, hi)
    ll, g = m.loglik_grad(x, nodes)
    pg, _, _ = projected_gradient(m, x, g)
    if len(fixed):
        pg[fixed] = 0.0
    return dict(x=x, ll=ll, pgrad=float(np.max(np.abs(pg * m.sc))), success=bool(r.success) or float(np.max(np.abs(pg * m.sc))) < 1e-3,
                message=str(r.message))


def boundary_report(m, x):
    out = []
    for i, name in enumerate(m.names):
        for side, val, kind in (("lower", m.lo[i], m.edge[i][0]), ("upper", m.hi[i], m.edge[i][1])):
            if abs(x[i] - val) <= 1e-8 * (1 + abs(val)):
                out.append(dict(parameter=name, side=side, value=float(val), kind=kind))
    return out


def verify_nodes(m, x, nodes):
    """Re-evaluate at doubled nodes; returns (ok, change, finer nodes)."""
    n2 = tuple(2 * k for k in nodes) if isinstance(nodes, tuple) else 2 * nodes
    a, b = m.loglik(m.unpack(x), nodes), m.loglik(m.unpack(x), n2)
    return abs(a - b) < QUAD_TOL, abs(a - b), n2


def multistart(m, rng, n_random, extra_starts=(), nodes=None):
    nodes = nodes or m.nodes
    starts = [deterministic_start(m), *extra_starts] + [random_start(m, rng) for _ in range(n_random)]
    fits = [local_fit(m, s, nodes) for s in starts]
    for _ in range(4):                      # escalate nodes until the optimum is accurate
        best = max((f for f in fits if np.isfinite(f["ll"])), key=lambda f: f["ll"])
        ok, change, n2 = verify_nodes(m, best["x"], nodes)
        if ok:
            break
        nodes = n2
        fits = [local_fit(m, f["x"], nodes) for f in fits]
    else:
        raise ArithmeticError(f"quadrature did not reach {QUAD_TOL} at the optimum")
    best = max(fits, key=lambda f: f["ll"])
    ties = [f for f in fits if best["ll"] - f["ll"] < MULTI_TOL]
    tie_targets = [m.all_targets(m.unpack(f["x"])) for f in ties]
    spread = {k: float(max(t[k] for t in tie_targets) - min(t[k] for t in tie_targets)) for k in m.targets}
    th = m.unpack(best["x"])
    return dict(x=best["x"].tolist(), ll=best["ll"], nodes=nodes, quad_change=change, pgrad=best["pgrad"],
                converged=best["pgrad"] < 1e-3, n_starts=len(starts), n_ties=len(ties),
                tie_spread=spread, all_ll=sorted([f["ll"] for f in fits], reverse=True)[:8],
                boundary=boundary_report(m, best["x"]), theta={k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in th.items()},
                targets=m.all_targets(th))


# --- profiles ----------------------------------------------------------------------------------------
def constrained_fit(m, name, q0, x0, nodes, rng, n_random=2):
    """max loglik subject to target(name) = q0, from x0 and n_random perturbed starts."""
    starts = [x0] + [np.clip(x0 + rng.normal(0, 0.05, len(x0)) * m.sc * (1 + np.abs(x0 / m.sc)), m.lo, m.hi)
                     for _ in range(n_random)]
    best = None
    if name in PARAMETER_TARGETS:
        i = m.names.index(name)
        if not (m.lo[i] <= q0 <= m.hi[i]):
            return None
        lo, hi = m.lo.copy(), m.hi.copy()
        lo[i] = hi[i] = q0
        for s in starts:
            s = s.copy(); s[i] = q0
            f = local_fit(m, s, nodes, lo, hi)
            if best is None or f["ll"] > best["ll"]:
                best = f
        best["residual"] = 0.0
        return best
    obj = _objective(m, nodes)
    con = dict(type="eq", fun=lambda xs: (m.target(name, m.unpack(xs * m.sc)) - q0) / max(1.0, abs(q0)))
    for s in starts:
        try:
            r = minimize(obj, s / m.sc, jac=True, method="SLSQP", bounds=list(zip(m.lo / m.sc, m.hi / m.sc)),
                         constraints=[con], options=dict(maxiter=2000, ftol=1e-13))
        except (ValueError, ArithmeticError):
            continue
        x = np.clip(r.x * m.sc, m.lo, m.hi)
        res = abs(m.target(name, m.unpack(x)) - q0)
        if res > RESIDUAL_TOL * max(1.0, abs(q0)):
            continue
        ll = m.loglik(m.unpack(x), nodes)
        if best is None or ll > best["ll"]:
            best = dict(x=x, ll=ll, residual=res, success=bool(r.success), message=str(r.message))
    return best


def target_range(m, name):
    """The values the target can take inside the support (parameters: their bounds; derived: None)."""
    if name in PARAMETER_TARGETS:
        i = m.names.index(name)
        return float(m.lo[i]), float(m.hi[i]), m.edge[i]
    return None, None, ("derived", "derived")


def profile(m, name, mle, rng, max_steps=40, scan=25):
    """Profile log-likelihood of one target, walked outward from the estimate in both directions until
    it falls 0.5 below the cutoff and stays below for four further, doubling steps (so a second,
    disconnected region above the cutoff would be seen), or the support / computational limit is
    reached. Crossings are refined by root finding. Returns the curve and the confidence set."""
    x_hat, ll_hat, nodes = np.array(mle["x"]), mle["ll"], mle["nodes"]
    q_hat = m.target(name, m.unpack(x_hat))
    lo_lim, hi_lim, edges = target_range(m, name)
    step0 = max(abs(q_hat) * 0.04, {"Tk": 1.0, "g": 0.0005, "sigma": 0.01, "tau": 0.05, "sigma_arm": 0.02,
                                       "gamma": 0.05, "lam": 0.002}.get(name, 1e-3 * (1 + abs(q_hat))))
    curve = [(q_hat, ll_hat, 0.0)]
    better = None
    for direction in (-1, 1):
        q, x, step, below, n = q_hat, x_hat.copy(), step0, 0, 0
        while n < max_steps:
            n += 1
            qn = q + direction * step
            lim = lo_lim if direction < 0 else hi_lim
            if lim is not None and (qn - lim) * direction > 0:
                qn = lim
            f = constrained_fit(m, name, qn, x, nodes, rng)
            if f is None:
                curve.append((qn, None, None)); break
            if f["ll"] > ll_hat + 1e-4:
                better = f                      # the "MLE" was not global: report, caller refits
            curve.append((qn, f["ll"], f.get("residual", 0.0)))
            drop = ll_hat - f["ll"]
            dprev = ll_hat - curve[-2][1] if curve[-2][1] is not None else 0.0
            step = step * (2.0 if abs(drop - dprev) < 0.15 else 0.5 if abs(drop - dprev) > 1.0 else 1.0)
            q, x = qn, f["x"]
            below = below + 1 if drop > CUTOFF + 0.5 else 0
            if below >= 1:
                step *= 2
            if below >= 5 or (lim is not None and q == lim):
                break
    # Coarse scan over the whole range, so a second region above the cutoff is seen wherever it lies
    # (a doubling walk can step over a narrow one). Parameters: their support (log-spaced when the
    # lower end is positive); derived targets: three times the walked span on each side.
    walked = [c[0] for c in curve if c[1] is not None]
    if lo_lim is not None:
        grid = (np.geomspace(lo_lim, hi_lim, scan) if lo_lim > 0 else np.linspace(lo_lim, hi_lim, scan))
    else:
        span = max(max(walked) - q_hat, q_hat - min(walked), step0)
        grid = np.linspace(q_hat - 3 * span, q_hat + 3 * span, scan)
    near = lambda q: min(abs(q - w) for w in walked) < 1e-9 * (1 + abs(q))
    for q in grid:
        if near(q):
            continue
        f = constrained_fit(m, name, float(q), x_hat, nodes, rng)
        curve.append((float(q), None if f is None else f["ll"], None if f is None else f.get("residual", 0.0)))
        if f is not None and f["ll"] > ll_hat + 1e-4:
            better = f
    pts = sorted((c for c in curve if c[1] is not None), key=lambda c: c[0])
    qs, lls = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
    inside = lls >= ll_hat - CUTOFF
    # Confidence set: maximal runs of grid points above the cutoff, ends refined by root finding.
    pieces, i = [], 0
    while i < len(qs):
        if not inside[i]:
            i += 1; continue
        j = i
        while j + 1 < len(qs) and inside[j + 1]:
            j += 1
        ends = []
        for a, b, side in ((i - 1, i, "lower"), (j, j + 1, "upper")):
            if 0 <= a and b < len(qs):
                xa = pts[a if side == "lower" else b]
                g = lambda q: (constrained_fit(m, name, q, x_hat, nodes, rng, 1) or dict(ll=-np.inf))["ll"] - (ll_hat - CUTOFF)
                try:
                    ends.append(dict(value=float(brentq(g, qs[a], qs[b], xtol=1e-6 * (1 + abs(q_hat)), maxiter=60)), closed=True))
                except ValueError:
                    ends.append(dict(value=float(qs[i] if side == "lower" else qs[j]), closed=False, why="root not bracketed"))
            else:
                k = i if side == "lower" else j
                lim = lo_lim if side == "lower" else hi_lim
                kind = edges[0] if side == "lower" else edges[1]
                at_lim = lim is not None and abs(qs[k] - lim) < 1e-12 * (1 + abs(lim))
                ends.append(dict(value=float(qs[k]), closed=False,
                                 why=(f"reaches the {EDGE_WORD.get(kind, kind)} limit" if at_lim else "not crossed within the search")))
        pieces.append(dict(lower=ends[0], upper=ends[1]))
        i = j + 1
    flat = bool(np.all(ll_hat - lls < CUTOFF))
    return dict(target=name, units=UNITS[name], estimate=q_hat, ll_hat=ll_hat, curve=[list(p) for p in pts],
                confidence_set=pieces, disconnected=len(pieces) > 1, flat=flat,
                better_point=None if better is None else dict(x=better["x"].tolist(), ll=better["ll"]))


# --- parametric bootstrap ------------------------------------------------------------------------
def simulate(m, x_hat, rng):
    """New counts from the fitted model: fresh random effects, fixed design and denominators.
    F-T: fresh arm effects with the fitted unit intercepts. F-E: fresh family and arm effects."""
    th, d = m.unpack(x_hat), m.d
    if m.kind == "tumor":
        eta = L.tumor_eta0(th, d) + th["sigma"] * rng.standard_normal(d["N"])
        y = rng.binomial(d["npat"], np.exp(-np.exp(eta)))
        return dict(d, y=y)
    u = rng.standard_normal(d["J"])
    eta = L.early_base(th, d) + th["tau"] * u[d["fam"]] + th["sigma_arm"] * rng.standard_normal(d["I"])
    return dict(d, y=rng.binomial(d["N"], expit(eta)))


def record_digest(rec):
    body = {k: v for k, v in rec.items() if k != "record_sha256"}
    return hashlib.sha256(json.dumps(body, sort_keys=True, default=float).encode()).hexdigest()


def write_atomic(path, obj):
    tmp = f"{path}.tmp{os.getpid()}"
    with open(tmp, "w") as f:
        json.dump(obj, f, sort_keys=True, default=float)
        f.flush(); os.fsync(f.fileno())
    os.replace(tmp, path)


def read_record(path, identity):
    """A completed record of this identity, or None (missing, malformed, stale or edited)."""
    try:
        rec = json.load(open(path))
    except (OSError, ValueError):
        return None
    if not isinstance(rec, dict) or rec.get("identity") != identity or rec.get("record_sha256") != record_digest(rec):
        return None
    if rec.get("status") not in ("ok", "failed"):
        return None
    return rec


def replicate_seed(b):
    return [SEED, 7, b]


def run_replicate(kind, support_name, x_hat, b, identity, out_dir, previous=None, max_attempts=3):
    """One bootstrap replicate: simulate (seeded by its index), refit, write its record. A failed
    previous record is retried with more starts and doubled nodes; its seed never changes."""
    m = Model(kind, support_name)
    attempt = (previous or {}).get("attempts", 0) + 1
    rng = np.random.default_rng(replicate_seed(b))
    dsim = simulate(m, np.array(x_hat), rng)
    ms = Model(kind, support_name, d=dsim)
    t0 = time.time()
    fit_rng = np.random.default_rng(replicate_seed(b) + [attempt])
    n_rand = START_POLICY["boot_random"] if attempt == 1 else START_POLICY["retry_random"]
    nodes = ms.nodes if attempt == 1 else (tuple(2 * k for k in ms.nodes) if isinstance(ms.nodes, tuple) else 2 * ms.nodes)
    try:
        f = multistart(ms, fit_rng, n_rand, extra_starts=[np.array(x_hat)], nodes=nodes)
        status = "ok" if f["converged"] else "failed"
        why = "" if f["converged"] else f"projected gradient {f['pgrad']:.2e}"
    except (ArithmeticError, FloatingPointError, ValueError) as err:
        f, status, why = None, "failed", f"{type(err).__name__}: {err}"
    rec = dict(identity=identity, replicate=b, seed=replicate_seed(b), attempts=attempt, status=status, why=why,
               generating="fitted MLE, fresh random effects, fixed design and denominators",
               y_sha256=hashlib.sha256(np.asarray(dsim["y"]).astype("<i8").tobytes()).hexdigest(),
               seconds=round(time.time() - t0, 2),
               fit=None if f is None else {k: f[k] for k in ("x", "ll", "nodes", "quad_change", "pgrad", "converged",
                                                              "n_starts", "n_ties", "boundary", "targets")},
               boundary_solution=bool(f and f["boundary"]))
    rec["record_sha256"] = record_digest(rec)
    write_atomic(os.path.join(out_dir, f"rep_{b:05d}.json"), rec)
    return rec


def bootstrap(kind, support_name, mle, B, out_dir, workers):
    os.makedirs(out_dir, exist_ok=True)
    identity = mle["identity"]
    todo = []
    done = failed_final = 0
    for b in range(1, B + 1):
        rec = read_record(os.path.join(out_dir, f"rep_{b:05d}.json"), identity)
        if rec is None:
            todo.append((b, None))
        elif rec["status"] == "ok":
            done += 1
        elif rec["attempts"] < 3:
            todo.append((b, rec))
        else:
            failed_final += 1
    print(f"bootstrap {kind}/{support_name}: {done} complete, {failed_final} failed after retries, {len(todo)} to run",
          flush=True)
    with cf.ProcessPoolExecutor(max_workers=workers) as ex:
        while todo:
            futs = {ex.submit(run_replicate, kind, support_name, mle["x"], b, identity, out_dir, prev): b for b, prev in todo}
            todo = []
            for fu in cf.as_completed(futs):
                rec = fu.result()
                print(f"  replicate {rec['replicate']} {rec['status']} attempt {rec['attempts']} {rec['seconds']} s", flush=True)
                if rec["status"] == "failed" and rec["attempts"] < 3:
                    todo.append((rec["replicate"], rec))
    recs = [read_record(os.path.join(out_dir, f"rep_{b:05d}.json"), identity) for b in range(1, B + 1)]
    if any(r is None for r in recs):
        raise RuntimeError("a replicate record is missing or does not verify; nothing assembled")
    return recs


def summarize_bootstrap(recs, mle):
    ok = [r for r in recs if r["status"] == "ok"]
    n_fail = len(recs) - len(ok)
    out = dict(B=len(recs), ok=len(ok), unresolved_failures=n_fail, failure_rate=n_fail / len(recs),
               provisional=n_fail / len(recs) > FAIL_LIMIT or n_fail > 0,
               boundary_solutions=sum(r["boundary_solution"] for r in ok), targets={})
    from scipy.stats import binom
    for k in mle["targets"]:
        v = np.sort(np.array([r["fit"]["targets"][k] for r in ok]))
        n = len(v)
        if n < 10:
            continue
        def q_ci(p):        # distribution-free 95% interval for the p-quantile from order statistics
            lo_i, hi_i = int(binom.ppf(0.025, n, p)), int(binom.ppf(0.975, n, p))
            return [float(v[max(lo_i - 1, 0)]), float(v[min(hi_i, n - 1)])]
        out["targets"][k] = dict(q025=float(np.quantile(v, 0.025)), q975=float(np.quantile(v, 0.975)),
                                 q025_mc95=q_ci(0.025), q975_mc95=q_ci(0.975),
                                 mean=float(v.mean()), bias=float(v.mean() - mle["targets"][k]),
                                 bias_mcse=float(v.std(ddof=1) / math.sqrt(n)), sd=float(v.std(ddof=1)))
    return out


def load_generation(bdir):
    """A complete, verified bootstrap generation: MANIFEST.json present, every listed record present
    with its recorded digest and the manifest's identity, and the summary matching. Anything else
    (an interrupted publication, a mixed or edited generation) is refused."""
    mp = os.path.join(bdir, "MANIFEST.json")
    if not os.path.exists(mp):
        raise RuntimeError(f"{bdir}: no MANIFEST.json; the generation is incomplete")
    man = json.load(open(mp))
    recs = []
    for fn, digest in man["records"].items():
        r = read_record(os.path.join(bdir, fn), man["identity"])
        if r is None or r["record_sha256"] != digest:
            raise RuntimeError(f"{bdir}/{fn}: missing, edited or of another identity; generation refused")
        recs.append(r)
    s = json.load(open(os.path.join(bdir, "summary.json")))
    if s.get("record_sha256") != man["summary_sha256"] or record_digest(s) != s["record_sha256"] or s.get("identity") != man["identity"]:
        raise RuntimeError(f"{bdir}: summary does not match its manifest; generation refused")
    return man, recs, s


def publish_bootstrap(kind, support_name, mle, B, bdir, workers, invocation):
    recs = bootstrap(kind, support_name, mle, B, bdir, workers)
    s = summarize_bootstrap(recs, mle)
    s.update(identity=mle["identity"], B=B, seconds_total=sum(r["seconds"] for r in recs))
    s["record_sha256"] = record_digest(s)
    write_atomic(os.path.join(bdir, "summary.json"), s)
    man = dict(identity=mle["identity"], records={f"rep_{r['replicate']:05d}.json": r["record_sha256"] for r in recs},
               summary_sha256=s["record_sha256"], invocation=invocation)
    write_atomic(os.path.join(bdir, "MANIFEST.json"), man)       # last: marks a complete generation
    return s


# --- CLI ---------------------------------------------------------------------------------------------
def gen_dir(kind, support_name, identity):
    return os.path.join(OUT, f"{kind}_{support_name}_{identity[:16]}")


def do_mle(kind, support_name):
    m = Model(kind, support_name)
    ident = run_identity(kind, m.d, support_name)
    out = gen_dir(kind, support_name, ident)
    os.makedirs(out, exist_ok=True)
    p = os.path.join(out, "mle.json")
    if os.path.exists(p):
        r = json.load(open(p))
        if r.get("identity") == ident and r.get("record_sha256") == record_digest(r):
            return r, out
    t0 = time.time()
    r = multistart(m, np.random.default_rng([SEED, 1]), START_POLICY["mle_random"])
    r.update(identity=ident, model=kind, support=support_name, seconds=round(time.time() - t0, 1),
             design_sha256=m.d["design_sha256"], file_sha256=m.d["file_sha256"])
    r["record_sha256"] = record_digest(r)
    write_atomic(p, r)
    return r, out


def do_profiles(kind, support_name, targets, workers):
    mle, out = do_mle(kind, support_name)
    m = Model(kind, support_name)
    targets = targets or list(m.targets)
    pdir = os.path.join(out, "profiles")
    os.makedirs(pdir, exist_ok=True)
    todo = []
    for t in targets:
        p = os.path.join(pdir, f"{t}.json")
        if os.path.exists(p):
            r = json.load(open(p))
            if r.get("identity") == mle["identity"] and r.get("record_sha256") == record_digest(r):
                continue
        todo.append(t)
    with cf.ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_profile_job, kind, support_name, mle, t, pdir): t for t in todo}
        for fu in cf.as_completed(futs):
            r = fu.result()
            print(f"  profile {r['target']}: {len(r['curve'])} points, {r['seconds']} s, "
                  f"set {[(p['lower']['value'], p['upper']['value']) for p in r['confidence_set']]}", flush=True)


def _profile_job(kind, support_name, mle, t, pdir, sigma_fixed=None):
    m = Model(kind, support_name)
    t0 = time.time()
    mm = mle
    if sigma_fixed is not None:      # hybrid sensitivity: sigma held at an externally chosen value
        i = m.names.index("sigma")
        m.lo[i] = m.hi[i] = sigma_fixed
        x = np.array(mle["x"]); x[i] = sigma_fixed
        f = local_fit(m, x, mle["nodes"])
        mm = dict(mle, x=f["x"].tolist(), ll=f["ll"])
    r = profile(m, t, mm, np.random.default_rng([SEED, 3, zlib.crc32(t.encode())]))
    r.update(identity=mle["identity"], seconds=round(time.time() - t0, 1), sigma_fixed=sigma_fixed)
    r["record_sha256"] = record_digest(r)
    name = t if sigma_fixed is None else f"{t}__sigma_{sigma_fixed:g}"
    write_atomic(os.path.join(pdir, f"{name}.json"), r)
    return r


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["validate", "mle", "profiles", "sigma_sensitivity", "bootstrap"])
    ap.add_argument("--model", choices=["tumor", "early"], default="tumor")
    ap.add_argument("--support", default=None)
    ap.add_argument("--targets", default="")
    ap.add_argument("--B", type=int, default=20)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    support = a.support or MODEL_SUPPORTS[a.model][0]
    workers = int(os.environ.get("FREQ_WORKERS", max(1, (os.cpu_count() or 2) - 2)))
    print(f"fit_ml {a.action} model={a.model} support={support} workers={workers} argv={sys.argv[1:]}", flush=True)
    if a.action == "validate":
        print(json.dumps(D.reproduce_bayesian_identities()))
        return
    if a.action == "mle":
        r, out = do_mle(a.model, support)
        print(json.dumps({k: r[k] for k in ("ll", "nodes", "quad_change", "pgrad", "n_ties", "boundary", "targets", "seconds")},
                         indent=1, default=float))
        print("written to", out)
        return
    if a.action == "profiles":
        do_profiles(a.model, support, [t for t in a.targets.split(",") if t], workers)
        return
    if a.action == "sigma_sensitivity":
        mle, out = do_mle(a.model, support)
        pdir = os.path.join(out, "profiles")
        os.makedirs(pdir, exist_ok=True)
        jobs = [(t, s) for s in (0.065, 0.17) for t in [f"contrast[{k}]" for k in range(10)]]
        with cf.ProcessPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(_profile_job, a.model, support, mle, t, pdir, s) for t, s in jobs]
            for fu in cf.as_completed(futs):
                r = fu.result()
                print(f"  sigma={r['sigma_fixed']} {r['target']} {r['seconds']} s", flush=True)
        return
    if a.action == "bootstrap":
        mle, out = do_mle(a.model, support)
        bdir = a.out or os.path.join(out, f"bootstrap_B{a.B}")
        s = publish_bootstrap(a.model, support, mle, a.B, bdir, workers, sys.argv)
        print(json.dumps({k: s[k] for k in ("B", "ok", "unresolved_failures", "boundary_solutions", "seconds_total")}))


if __name__ == "__main__":
    main()
