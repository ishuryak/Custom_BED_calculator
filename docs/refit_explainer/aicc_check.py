"""AICc of DI vs DD by maximum likelihood (binomial, one intercept per fitting unit,
alpha/beta 10 Gy, g fixed 0.0154/day as in the 2018 paper; optional g free).
Same model form as fit/mle_check.py, applied to a given dataset with fraction components."""
import csv, json, sys, numpy as np
from scipy.optimize import minimize
path, gmode = sys.argv[1], sys.argv[2]   # gmode: fixed | free
R = list(csv.DictReader(open(path, encoding="utf-8")))
T = np.array([float(r["total_time_days"]) for r in R]); N = np.array([int(r["n_patients"]) for r in R])
y = np.round(N * np.array([float(r["tcp_actual_pct"]) for r in R]) / 100)
if "total_dose_gy" in R[0]:
    ds = np.array([float(r["total_dose_gy"]) for r in R]); d2 = np.array([float(r["sum_d2_gy2"]) for r in R])
else:
    d = np.array([float(r["dose_per_fraction_gy"]) for r in R]); n = np.array([float(r["n_fractions"]) for r in R]); ds, d2 = n*d, n*d*d
K = ds + d2 / 10
tr = [r["ref_no"] for r in R]; ut = sorted(set(tr), key=tr.index); j = np.array([ut.index(t) for t in tr]); J = len(ut)
sp = lambda x: np.where(3 * x > 30, x, np.log1p(np.exp(np.minimum(3 * x, 30))) / 3)
def g_of(th): return 0.0154 if gmode == "fixed" else 0.0046 + (0.021-0.0046)/(1+np.exp(-th[3]))
NS = 3 if gmode == "fixed" else 4
def eta(th, model):
    a, lg, third = th[:3]; c = th[NS:]
    if model == "di": rep = np.exp(lg) * sp(T - third)
    else: rep = np.exp(lg) * (1 - np.exp(-a * K / T)) * sp(T - np.exp(third) / a * T / K)
    return c[j] - a * K + rep + g_of(th) * T
def nll(th, model):
    if th[0] <= 0: return 1e10
    P = np.clip(np.exp(-np.exp(eta(th, model))), 1e-12, 1 - 1e-12)
    return -(y * np.log(P) + (N - y) * np.log1p(-P)).sum()
out = dict(data=path, g=gmode, n_arms=len(y), n_units=J)
for model, x0 in (("di", [0.07, np.log(0.035), 28.6]), ("dd", [0.224, np.log(1.17), np.log(14.5)])):
    best = None
    for s in range(40):
        rng = np.random.default_rng(s)
        t0 = np.array(x0) * (1 + (0.3 * rng.normal(size=3) if s else 0))
        if gmode == "free": t0 = np.r_[t0, rng.normal() if s else 0.0]
        g0 = g_of(t0)
        c0 = [np.log(-np.log(np.clip(y[j == k].sum() / N[j == k].sum(), .05, .95))) + t0[0] * K[j == k].mean() - g0 * T[j == k].mean() for k in range(J)]
        r = minimize(nll, np.r_[t0, c0], args=(model,), method="Nelder-Mead", options=dict(maxiter=400000, maxfev=400000, xatol=1e-9, fatol=1e-11))
        r = minimize(nll, r.x, args=(model,), method="BFGS")
        if best is None or r.fun < best.fun: best = r
    k = NS + J; n = len(y)
    P = np.exp(-np.exp(eta(best.x, model)))
    X2 = float(((y - N * P) ** 2 / (N * P * (1 - P))).sum())
    aic = 2 * best.fun + 2 * k; aicc = aic + 2 * k * (k + 1) / (n - k - 1)
    out[model] = dict(neg_log_lik=float(best.fun), k=k, AIC=float(aic), AICc=float(aicc), pearson_x2=X2, df=n-k, chat=X2/(n-k),
                      alpha=float(best.x[0]), lam=float(np.exp(best.x[1])), third=float(best.x[2] if model=="di" else np.exp(best.x[2])), g=float(g_of(best.x)))
out["dAICc_DI_minus_DD"] = out["di"]["AICc"] - out["dd"]["AICc"]
chat = max(out["di"]["chat"], out["dd"]["chat"])
out["chat_used"] = chat
for m in ("di", "dd"):
    k = out[m]["k"]; n = out["n_arms"]
    out[m]["QAICc"] = 2*out[m]["neg_log_lik"]/chat + 2*(k+1) + 2*(k+1)*(k+2)/(n-k-2)
out["dQAICc_DI_minus_DD"] = out["di"]["QAICc"] - out["dd"]["QAICc"]
print(json.dumps(out, indent=1))
