"""Exploratory Bayesian refit of the Fowler acute-mucosa boundary to the
Pardo-Montero 2021 supporting data. P(intolerable) = logistic((BED10(ab) -
r*max(0, Te - Tk) - B0)/s), Te = elapsed days. Random-walk Metropolis."""
import numpy as np, pickle, sys, json
S = pickle.load(open(sys.argv[1], 'rb')); prior_kind = sys.argv[2]; seed = int(sys.argv[3])
S = [s for s in S if s['cls'] != 'Disp']
if prior_kind.endswith('_nodisp'): S = [s for s in S if not s['cls'].startswith('Disp')]
mu = np.log(2)/3.2
D = np.array([s['d'].sum() for s in S])
Q = np.array([s['d'] @ np.exp(-mu*np.abs(s['t'][:,None]-s['t'][None,:])) @ s['d'] for s in S])
Te = np.array([s['T']-1 for s in S]); y = np.array([s['y'] for s in S])
print(f"n={len(y)} intolerable={y.sum()}", file=sys.stderr)
def smax(x, k=3.0): return np.where(k*x > 30, x, np.log1p(np.exp(np.minimum(k*x, 30)))/k)
# theta = B0, Tk, log r, log s, log ab
def logpost(th):
    B0, Tk, lr, ls, lab = th
    if not (30 < B0 < 120 and 0 <= Tk <= 40): return -np.inf
    r, s, ab = np.exp(lr), np.exp(ls), np.exp(lab)
    lp = 0.0
    lp += -0.5*((lab-np.log(10))/0.354)**2                    # ab: median 10, 95% 5-20 (Denham 1995: 9.3, 5.8-17.9)
    if prior_kind.startswith('informed'):
        # Bentzen 2001 CHART: Dprolif 0.80 (0.7-1.1) Gy/day in 2-Gy-equivalents -> BED10 x1.2: 0.96 (0.84-1.32)
        lp += -0.5*((lr-np.log(0.96))/0.12)**2
    else:
        lp += -0.5*((lr-np.log(0.8))/1.0)**2                   # weak: 95% 0.11-5.7 Gy10/day
    lp += -0.5*(s/5.0)**2 + ls                                 # half-normal(5) on s, Jacobian
    z = (D + Q/ab - r*smax(Te-Tk) - B0)/s
    ll = np.sum(y*(-np.logaddexp(0, -z)) + (1-y)*(-np.logaddexp(0, z)))
    return lp + ll
rng = np.random.default_rng(seed)
def run(th0, n=60000, burn=20000):
    th = np.array(th0, float); lp = logpost(th); cov = np.diag([2.0, 3.0, 0.2, 0.3, 0.2])**2
    out = []; acc = 0; hist = []
    for i in range(n):
        if i > 2000 and i % 1000 == 0 and i < burn:
            cov = np.cov(np.array(hist[-5000:]).T)*2.38**2/5 + 1e-8*np.eye(5)
        prop = rng.multivariate_normal(th, cov); lpp = logpost(prop)
        if np.log(rng.uniform()) < lpp - lp: th, lp = prop, lpp; acc += 1
        hist.append(th.copy())
        if i >= burn and i % 10 == 0: out.append(th.copy())
    return np.array(out), acc/n
chains = []
for c, init in enumerate([[67, 7, np.log(0.8), 0.5, np.log(10)], [70, 18, np.log(1.0), 1.0, np.log(10)],
                          [62, 3, np.log(0.6), 1.5, np.log(8)], [75, 25, np.log(1.3), 0.0, np.log(12)]]):
    ch, a = run(init); chains.append(ch); print(f"chain {c} accept {a:.2f}", file=sys.stderr)
C = np.array(chains)  # chains x draws x 5
# split R-hat
def rhat(x):
    m, n = x.shape; h = n//2; x = np.concatenate([x[:, :h], x[:, h:2*h]]); m, n = x.shape
    W = x.var(axis=1, ddof=1).mean(); B = n*x.mean(axis=1).var(ddof=1); return np.sqrt(((n-1)/n*W + B/n)/W)
names = ['B0', 'Tk_elapsed', 'r=lambda/alpha', 's', 'ab']
draws = C.reshape(-1, 5).copy(); draws[:, 2:] = np.exp(draws[:, 2:])
for j, nm in enumerate(names):
    q = np.percentile(draws[:, j], [2.5, 50, 97.5])
    print(f"{nm:16s} median {q[1]:7.2f}  95% {q[0]:7.2f}-{q[2]:7.2f}  Rhat {rhat(C[:, :, j]):.3f}")
np.save(sys.argv[4], draws)
