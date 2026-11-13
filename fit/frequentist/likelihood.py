"""Marginal likelihoods for the likelihood-only (frequentist) comparison, plan revision 3
(docs/development/2026_11_12_frequentist_comparison_plan.md). Pure numerical functions, no I/O.

F-T (tumor DI, fit/stan/di.stan without priors). Arm i of fitting unit j, failures f = N - y:
    eta_i = c_j - alpha (B_i - Bbar_j) + lam smax(T_i - Tk) + g (T_i - Tbar_j) + sigma z_i,
    B = dsum + d2sum / ab,  TCP = exp(-exp(eta)),  z_i ~ N(0, 1) integrated out.
F-E (early, tools/early_tox/step5_model.stan unit U02 without priors). Arm i of family j:
    logit p_i = m0 + grp_j . delta + tau u_j + kappa (dD + dQ / ab - r dTe) + gamma chemo + s_arm z_i,
    u_j, z_i ~ N(0, 1) integrated out (z nested inside u).

Every log-likelihood here is the FULL log-likelihood, including log C(N, y), as Stan's binomial_lpmf.

Integration: adaptive Gauss-Hermite quadrature centered at the mode of each integrand and scaled by
its curvature (Liu and Pierce 1994; review F01: ordinary fixed-node Gauss-Hermite is off by up to
-6.9 log units on these trial sizes). The F-E outer integrand is the SUM of the inner log integrals
of a family's arms. Zero variances use the exact closed form. Accuracy is not assumed: callers
evaluate at n and 2n nodes (`converged_loglik`) and accept only a change below the registered
tolerance.
"""
import math

import numpy as np
from scipy.special import gammaln, log_expit, logsumexp

LOG_SQRT_2PI = 0.5 * math.log(2 * math.pi)
SMOOTH_K = 3.0          # smax sharpness (1/day), as fit/stan/common_functions.stanfunctions and R SMOOTH_K
NEWTON_ITERS = 50


def log_choose(n, k):
    n, k = np.asarray(n, float), np.asarray(k, float)
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


def smax(x):
    """log(1 + exp(3x)) / 3, evaluated stably (the tumor fit's smooth hinge)."""
    return np.logaddexp(0.0, SMOOTH_K * np.asarray(x, float)) / SMOOTH_K


def log1mexp(a):
    """log(1 - exp(-a)) for a > 0, stable at both ends (Maechler 2012)."""
    a = np.asarray(a, float)
    return np.where(a < math.log(2), np.log(-np.expm1(-a)), np.log1p(-np.exp(-a)))


# --- observation models: log density and first two derivatives in the linear predictor -------------
def cloglog_terms(eta, y, f):
    """Tumor arm: y controlled, f failed, P(control) = exp(-e^eta). Returns (l, l1, l2) without log C."""
    mu = np.exp(np.minimum(eta, 700.0))
    em1 = np.expm1(mu)
    l = -y * mu + f * log1mexp(mu)
    r = np.where(mu < 1e-300, 1.0, mu / np.where(em1 == 0, 1.0, em1))     # mu / (e^mu - 1)
    l1 = -y * mu + f * r
    q = np.where(mu < 1e-300, 0.0, mu * mu / (em1 * (-np.expm1(-mu))))    # mu^2 e^mu / (e^mu - 1)^2
    l2 = -y * mu + f * (r - q)
    return l, l1, l2


def logit_terms(eta, y, n):
    """Early arm: y events of n, logit link. Returns (l, l1, l2) without log C."""
    l = y * log_expit(eta) + (n - y) * log_expit(-eta)
    p = np.exp(log_expit(eta))
    return l, y - n * p, -n * p * (1 - p)


# --- one-dimensional adaptive quadrature over z ~ N(0, 1) entering as eta0 + s z ----------------------
def bracketed_newton(grad_curv, shape, max_iter=200, tol=1e-12):
    """Root of a strictly decreasing derivative g(z) (the integrand is log-concave), elementwise.
    grad_curv(z) -> (g, g') with g' < 0. Newton steps inside a bracket [lo, hi] with g(lo) > 0 > g(hi),
    bisection whenever Newton would leave it. A clipped plain Newton cycled between two points on
    real arms (2026-11-12, the engine's first draft); this cannot."""
    lo, hi = np.full(shape, -1.0), np.full(shape, 1.0)
    for _ in range(80):                                   # expand until the root is bracketed
        glo, _ = grad_curv(lo)
        ghi, _ = grad_curv(hi)
        bad_lo, bad_hi = glo <= 0, ghi >= 0
        if not (bad_lo.any() or bad_hi.any()):
            break
        lo = np.where(bad_lo, 2 * lo, lo)
        hi = np.where(bad_hi, 2 * hi, hi)
    else:
        raise ArithmeticError("mode not bracketed within 2^80")
    z = np.zeros(shape)
    for _ in range(max_iter):
        g, gp = grad_curv(z)
        lo = np.where(g > 0, z, lo)
        hi = np.where(g <= 0, z, hi)
        zn = z - g / gp
        zn = np.where((zn > lo) & (zn < hi), zn, 0.5 * (lo + hi))
        done = np.abs(zn - z) < tol * (1 + np.abs(z))
        z = zn
        if np.all(done) or np.all(hi - lo < tol):
            break
    return z


def _mode(eta0, s, obs, *args):
    """Mode and curvature scale of h(z) = l(eta0 + s z) - z^2/2 (log-concave: l'' <= 0)."""
    shape = np.broadcast(eta0, s, *args).shape

    def gc(z):
        _, l1, l2 = obs(eta0 + s * z, *args)
        return s * l1 - z, s * s * l2 - 1.0
    z = bracketed_newton(gc, shape)
    _, _, l2 = obs(eta0 + s * z, *args)
    return z, 1.0 / np.sqrt(1.0 - s * s * l2)


def agq_1d(eta0, s, obs, args, nodes, with_z=False):
    """log of the integral over z of exp(l(eta0 + s z)) phi(z), elementwise, with adaptive Gauss-Hermite.
    Also returns d/d eta0 of that log integral (posterior mean of l'), used for the outer Newton steps."""
    x, w = np.polynomial.hermite.hermgauss(nodes)
    eta0, s = np.broadcast_arrays(np.asarray(eta0, float), np.asarray(s, float))
    zero = s == 0
    zh, sh = _mode(eta0, np.where(zero, 0.0, s), obs, *args)
    zk = zh[..., None] + math.sqrt(2) * sh[..., None] * x                  # nodes per element
    l, l1, l2 = obs(eta0[..., None] + s[..., None] * zk, *[np.asarray(a, float)[..., None] for a in args])
    logterm = np.log(w) + x * x + l - 0.5 * zk * zk - LOG_SQRT_2PI
    out = logsumexp(logterm, axis=-1) + np.log(math.sqrt(2) * sh)
    pw = np.exp(logterm - logsumexp(logterm, axis=-1, keepdims=True))
    m1 = (pw * l1).sum(-1)
    m2 = (pw * (l2 + l1 * l1)).sum(-1) - m1 * m1                         # d2 log L / d eta0^2
    l0, l10, l20 = obs(eta0, *args)
    out = np.where(zero, l0, out)
    if with_z:          # E[l' z]: d log L / d s (Fisher identity); 0 at s = 0, where the posterior is the prior
        return out, np.where(zero, l10, m1), np.where(zero, l20, m2), np.where(zero, 0.0, (pw * l1 * zk).sum(-1))
    return out, np.where(zero, l10, m1), np.where(zero, l20, m2)


# --- F-T ------------------------------------------------------------------------------------------
def tumor_eta0(th, d):
    """Linear predictor without the arm effect. th: dict(c, alpha, ab, lam, Tk, g)."""
    B = d["dsum"] + d["d2sum"] / th["ab"]
    J = d["J"]
    cnt = np.bincount(d["unit"], minlength=J)
    Bbar = np.bincount(d["unit"], B, minlength=J) / cnt
    Tbar = np.bincount(d["unit"], d["T"], minlength=J) / cnt
    return (th["c"][d["unit"]] - th["alpha"] * (B - Bbar[d["unit"]]) + th["lam"] * smax(d["T"] - th["Tk"])
            + th["g"] * (d["T"] - Tbar[d["unit"]]))


def tumor_loglik(th, d, nodes=20, per_arm=False):
    """Full marginal log-likelihood of F-T (sum over arms, or per arm)."""
    y, n = d["y"].astype(float), d["npat"].astype(float)
    ll, _, _ = agq_1d(tumor_eta0(th, d), th["sigma"], cloglog_terms, (y, n - y), nodes)
    ll = ll + log_choose(n, y)
    return ll if per_arm else float(ll.sum())


TUMOR_SCALARS = ("alpha", "ab", "lam", "Tk", "g", "sigma")


def tumor_loglik_grad(th, d, nodes=20):
    """(log-likelihood, gradient dict) of F-T. The gradient is the posterior mean of the complete-data
    score at the quadrature weights (Fisher identity); tests compare it with finite differences."""
    y, n = d["y"].astype(float), d["npat"].astype(float)
    ll, m1, _, m1z = agq_1d(tumor_eta0(th, d), th["sigma"], cloglog_terms, (y, n - y), nodes, with_z=True)
    u, J = d["unit"], d["J"]
    cnt = np.bincount(u, minlength=J)
    cen = lambda v: v - (np.bincount(u, v, minlength=J) / cnt)[u]
    B = d["dsum"] + d["d2sum"] / th["ab"]
    deta = dict(alpha=-cen(B), ab=th["alpha"] * cen(d["d2sum"]) / th["ab"] ** 2, lam=smax(d["T"] - th["Tk"]),
                Tk=-th["lam"] * np.exp(log_expit(SMOOTH_K * (d["T"] - th["Tk"]))), g=cen(d["T"]))
    grad = {k: float((m1 * v).sum()) for k, v in deta.items()}
    grad["sigma"] = float(m1z.sum())
    grad["c"] = np.bincount(u, m1, minlength=J)
    return float((ll + log_choose(n, y)).sum()), grad


def tumor_contrasts(th, d):
    """Effective tumor BED (Gy at the parameters' alpha/beta) of each evaluation regimen minus the
    anchor 35 x 2 Gy in 47 days: delta_eval of fit/stan/di.stan."""
    def eff(dsum, d2sum, T):
        return dsum + d2sum / th["ab"] - (th["lam"] * smax(T - th["Tk"]) + th["g"] * T) / th["alpha"]
    return eff(d["eval_dsum"], d["eval_d2sum"], d["eval_T"]) - eff(d["anchor_dsum"], d["anchor_d2sum"], d["anchor_T"])


# --- F-E ------------------------------------------------------------------------------------------
def early_base(th, d):
    """Arm linear predictor without the family and arm effects."""
    lin = th["m0"] + d["grp"][d["fam"]] @ th["delta"] if d["G"] else np.full(d["I"], th["m0"])
    E = d["dD"] + d["dQ"] / th["ab"] - th["r"] * d["dTe"]
    return lin + th["kappa"] * E + th["gamma"] * d["chemo"]


def early_loglik(th, d, nodes=(20, 20), per_family=False):
    """Full marginal log-likelihood of F-E: for each family, the outer adaptive integral over u of
    exp(sum over its arms of the inner log integral over z)."""
    n_out, n_in = nodes
    y, N = d["y"].astype(float), d["N"].astype(float)
    base, fam, J, tau, sa = early_base(th, d), d["fam"], d["J"], th["tau"], th["sigma_arm"]

    def family_terms(u):                       # u: (J,) or (J, K) -> sum over arms of inner log L, + derivs
        uu = u[fam]
        li, g1, g2 = agq_1d(base.reshape((-1,) + (1,) * (uu.ndim - 1)) + tau * uu, sa, logit_terms,
                            (y.reshape(base.shape + (1,) * (uu.ndim - 1)), N.reshape(base.shape + (1,) * (uu.ndim - 1))),
                            n_in)
        S = np.zeros((J,) + uu.shape[1:])
        np.add.at(S, fam, li)
        D1 = np.zeros_like(S); np.add.at(D1, fam, tau * g1)
        D2 = np.zeros_like(S); np.add.at(D2, fam, tau * tau * g2)
        return S, D1, D2

    if tau == 0:
        S, _, _ = family_terms(np.zeros(J))
        fam_ll = S
    else:
        def gc(u):                             # H(u) = S(u) - u^2/2 is concave (Prekopa), so H' decreases
            _, D1, D2 = family_terms(u)
            return D1 - u, D2 - 1.0
        u = bracketed_newton(gc, (J,), tol=1e-10)
        _, _, D2 = family_terms(u)
        sh = 1.0 / np.sqrt(np.maximum(1.0 - D2, 1e-12))
        x, w = np.polynomial.hermite.hermgauss(n_out)
        uk = u[:, None] + math.sqrt(2) * sh[:, None] * x
        S, _, _ = family_terms(uk)
        fam_ll = logsumexp(np.log(w) + x * x + S - 0.5 * uk * uk - LOG_SQRT_2PI, axis=1) + np.log(math.sqrt(2) * sh)
    fam_ll = fam_ll + np.bincount(fam, log_choose(N, y), minlength=J)
    return fam_ll if per_family else float(fam_ll.sum())


def early_loglik_grad(th, d, nodes=(20, 20)):
    """(log-likelihood, gradient dict) of F-E by the Fisher identity over the nested quadrature."""
    n_out, n_in = nodes
    y, N = d["y"].astype(float), d["N"].astype(float)
    base, fam, J, tau, sa = early_base(th, d), d["fam"], d["J"], th["tau"], th["sigma_arm"]
    E = d["dD"] + d["dQ"] / th["ab"] - th["r"] * d["dTe"]
    dbase = dict(m0=np.ones(d["I"]), kappa=E, ab=-th["kappa"] * d["dQ"] / th["ab"] ** 2,
                 r=-th["kappa"] * d["dTe"], gamma=d["chemo"])
    if tau == 0:
        uk, wfam = np.zeros((J, 1)), np.ones((J, 1))
    else:
        def gc(u):
            uu = u[fam]
            _, g1, g2 = agq_1d(base + tau * uu, sa, logit_terms, (y, N), n_in)
            D1 = np.zeros(J); np.add.at(D1, fam, tau * g1)
            D2 = np.zeros(J); np.add.at(D2, fam, tau * tau * g2)
            return D1 - u, D2 - 1.0
        uh = bracketed_newton(gc, (J,), tol=1e-10)
        sh = 1.0 / np.sqrt(np.maximum(-gc(uh)[1], 1e-12))
        x, w = np.polynomial.hermite.hermgauss(n_out)
        uk = uh[:, None] + math.sqrt(2) * sh[:, None] * x
    K = uk.shape[1]
    li, m1, _, m1z = agq_1d(base[:, None] + tau * uk[fam], sa, logit_terms, (y[:, None], N[:, None]), n_in, with_z=True)
    S = np.zeros((J, K)); np.add.at(S, fam, li)
    if tau == 0:
        fam_ll, pk = S[:, 0], np.ones((J, 1))
    else:
        lt = np.log(w) + x * x + S - 0.5 * uk * uk - LOG_SQRT_2PI
        fam_ll = logsumexp(lt, axis=1) + np.log(math.sqrt(2) * sh)
        pk = np.exp(lt - logsumexp(lt, axis=1, keepdims=True))      # outer posterior weights per family
    pa = pk[fam]                                                      # (I, K)
    em1 = (pa * m1).sum(1)                                            # E[l'] per arm
    grad = {k: float((em1 * v).sum()) for k, v in dbase.items()}
    grad["delta"] = (d["grp"][fam] * em1[:, None]).sum(0) if d["G"] else np.zeros(0)
    grad["tau"] = float((pa * m1 * uk[fam]).sum())
    grad["sigma_arm"] = float((pa * m1z).sum())
    ll = float(fam_ll.sum() + log_choose(N, y).sum())
    return ll, grad


def early_contrasts(th, d):
    """Log odds ratio of each decision schedule against 70 Gy / 35 (lor_dec of step5_model.stan)."""
    return th["kappa"] * (d["kD"] + d["kQ"] / th["ab"] - th["r"] * d["kTe"])


# --- accuracy ladder --------------------------------------------------------------------------------
def converged_loglik(fn, th, d, start, tol=1e-6, max_nodes=320):
    """Evaluate fn at n and 2n nodes from `start`, doubling until the total changes by less than tol.
    Returns (value at the finer level, nodes used, change). Raises if max_nodes is reached first."""
    n = start
    prev = fn(th, d, n)
    while True:
        n2 = tuple(2 * k for k in n) if isinstance(n, tuple) else 2 * n
        if (max(n2) if isinstance(n2, tuple) else n2) > max_nodes:
            raise ArithmeticError(f"quadrature did not converge to {tol} by {max_nodes} nodes")
        cur = fn(th, d, n2)
        if abs(cur - prev) < tol:
            return cur, n2, abs(cur - prev)
        n, prev = n2, cur
