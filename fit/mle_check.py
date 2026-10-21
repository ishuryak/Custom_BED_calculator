"""Maximum-likelihood check of overdispersion (binomial, per-trial intercepts,
alpha/beta 10, g 0.0154 as in the 2018 paper). Writes fit/output/mle_check.json."""
import csv, json, numpy as np
from scipy.optimize import minimize
R = list(csv.DictReader(open("data/shuryak2018_table1.csv")))
d = np.array([float(r["dose_per_fraction_gy"]) for r in R]); n = np.array([float(r["n_fractions"]) for r in R])
T = np.array([float(r["total_time_days"]) for r in R]); N = np.array([int(r["n_patients"]) for r in R])
y = np.round(N * np.array([float(r["tcp_actual_pct"]) for r in R]) / 100)
tr = [r["ref_no"] for r in R]; ut = sorted(set(tr), key=tr.index); j = np.array([ut.index(t) for t in tr]); J = len(ut)
sp = lambda x: np.where(3 * x > 30, x, np.log1p(np.exp(np.minimum(3 * x, 30))) / 3)
K = n * d * (1 + d / 10)
def eta(th, model):
    a, lg, third = th[:3]; c = th[3:]
    if model == "di": rep = np.exp(lg) * sp(T - third)
    else: rep = np.exp(lg) * (1 - np.exp(-a * K / T)) * sp(T - np.exp(third) / a * T / K)
    return c[j] - a * K + rep + 0.0154 * T
def nll(th, model):
    P = np.clip(np.exp(-np.exp(eta(th, model))), 1e-12, 1 - 1e-12)
    return -(y * np.log(P) + (N - y) * np.log1p(-P)).sum()
out = {}
for model, x0 in (("di", [0.07, np.log(0.035), 28.6]), ("dd", [0.224, np.log(1.17), np.log(14.5)])):
    best = None
    for s in range(30):
        rng = np.random.default_rng(s)
        t0 = np.array(x0) * (1 + (0.3 * rng.normal(size=3) if s else 0))
        c0 = [np.log(-np.log(np.clip(y[j == k].sum() / N[j == k].sum(), .05, .95))) + t0[0] * K[j == k].mean() - 0.0154 * T[j == k].mean() for k in range(J)]
        r = minimize(nll, np.r_[t0, c0], args=(model,), method="Nelder-Mead", options=dict(maxiter=200000, maxfev=200000, xatol=1e-8, fatol=1e-10))
        r = minimize(nll, r.x, args=(model,), method="BFGS")
        if best is None or r.fun < best.fun: best = r
    P = np.exp(-np.exp(eta(best.x, model)))
    X2 = float(((y - N * P) ** 2 / (N * P * (1 - P))).sum()); df = len(y) - (J + 3)
    out[model] = dict(neg_log_lik=float(best.fun), alpha=float(best.x[0]), lam=float(np.exp(best.x[1])),
                      third=float(best.x[2] if model == "di" else np.exp(best.x[2])), pearson_x2=X2, df=df, ratio=X2 / df)
out["loglik_dd_minus_di"] = out["di"]["neg_log_lik"] - out["dd"]["neg_log_lik"]
json.dump(out, open("fit/output/mle_check.json", "w"), indent=1)
print(json.dumps(out, indent=1))
