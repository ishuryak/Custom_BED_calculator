"""Test T6: the independent R cross-checks of plan revision 3 (fit/frequentist/glmer_check.R).

    python3 fit/frequentist/r_check.py [--out DIR]

Writes the fixture, runs Windows Rscript (FREQ_RSCRIPT overrides the path) in that folder, compares
and writes r_check.json there. Exit 1 if any required comparison fails, and also if R or lme4 is
unavailable: a check that could not run has not passed.
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np
from scipy.special import expit

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from frequentist import fit_ml as FM  # noqa: E402
from frequentist import likelihood as L  # noqa: E402

RSCRIPT = os.environ.get("FREQ_RSCRIPT", "/mnt/c/Program Files/R/R-4.5.1/bin/Rscript.exe")
COEF_ABS, COEF_REL, OBJ_TOL = 1e-6, 1e-4, 1e-4
SYNTH = dict(m0=-0.5, delta=np.array([0.3]), tau=0.0, kappa=0.05, ab=10.0, r=1.0, gamma=0.5, sigma_arm=0.0)


def close(a, b):
    return abs(a - b) <= COEF_ABS + COEF_REL * abs(b)


def fixed_tumor_fit(m, ab, Tk):
    """Python fit at fixed (ab, Tk) on the glmer support: alpha, lambda and g unconstrained."""
    m.lo, m.hi = m.lo.copy(), m.hi.copy()
    for k, (lo, hi) in dict(alpha=(-2, 2), lam=(-2, 2), g=(-1, 1)).items():
        i = m.names.index(k); m.lo[i], m.hi[i] = lo, hi
    for k, v in dict(ab=ab, Tk=Tk).items():
        i = m.names.index(k); m.lo[i] = m.hi[i] = v
    x0 = FM.deterministic_start(m); x0[m.names.index("ab")], x0[m.names.index("Tk")] = ab, Tk
    return FM.local_fit(m, x0, 40)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(FM.OUT, "r_check"))
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    mle, _ = FM.do_mle("tumor", "matched")
    m = FM.Model("tumor", "matched")
    d = m.d
    ab, Tk = mle["theta"]["ab"], mle["theta"]["Tk"]
    py = fixed_tumor_fit(m, ab, Tk)
    th = m.unpack(py["x"])
    tumor = dict(unit=d["unit"].tolist(), B=(d["dsum"] + d["d2sum"] / ab).tolist(), hinge=L.smax(d["T"] - Tk).tolist(),
                 T=d["T"].tolist(), fail=(d["npat"] - d["y"]).tolist(), ctrl=d["y"].tolist())
    # Early GLM fixture: synthetic counts from known interior parameters on the real design.
    me = FM.Model("early", "early_matched")
    rng = np.random.default_rng([FM.SEED, 6])
    y_syn = rng.binomial(me.d["N"], expit(L.early_base(SYNTH, me.d)))
    eg = dict(y=y_syn.tolist(), N=me.d["N"].tolist(), grp=me.d["grp"][me.d["fam"], 0].tolist(), dD=me.d["dD"].tolist(),
              dQ=me.d["dQ"].tolist(), dTe=me.d["dTe"].tolist(), chemo=me.d["chemo"].tolist())
    er = dict(eg, y=me.d["y"].tolist(), fam=me.d["fam"].tolist())
    json.dump(dict(tumor=tumor, early_glm=eg, early_real=er), open(os.path.join(a.out, "r_fixture.json"), "w"))
    if not os.path.exists(RSCRIPT):
        sys.exit(f"REQUIRED CHECK NOT RUN: Rscript not found at {RSCRIPT}")
    script = os.path.join(HERE, "glmer_check.R")
    with open(script) as f, open(os.path.join(a.out, "glmer_check.R"), "w") as g:
        g.write(f.read())
    p = subprocess.run([RSCRIPT, "glmer_check.R"], cwd=a.out, capture_output=True, text=True)
    if p.returncode != 0 or not os.path.exists(os.path.join(a.out, "r_result.json")):
        sys.exit(f"REQUIRED CHECK NOT RUN: R failed ({p.returncode}): {p.stderr[-2000:]}")
    r = json.load(open(os.path.join(a.out, "r_result.json")))
    res, ok = dict(R=r["R"], lme4=r["lme4"], python=dict(ab=ab, Tk=Tk)), True

    # Tumor coefficients and the objective at R's coefficients, on Python's own likelihood.
    rt = r["tumor"]
    comp = {k: dict(python=float(th[k]), R=float(rt[k]), ok=close(float(th[k]), float(rt[k]))) for k in ("alpha", "lam", "g")}
    units = [int(u) for u in rt["unit_levels"]]
    J = d["J"]
    cnt = np.bincount(d["unit"], minlength=J)
    Bbar = np.bincount(d["unit"], d["dsum"] + d["d2sum"] / ab, minlength=J) / cnt
    Tbar = np.bincount(d["unit"], d["T"], minlength=J) / cnt
    c = np.zeros(J)
    for u, v in zip(units, rt["intercepts"]):
        c[u] = v - rt["alpha"] * Bbar[u] + rt["g"] * Tbar[u]
    th_r = dict(th, alpha=rt["alpha"], lam=rt["lam"], g=rt["g"], sigma=rt["sigma"], c=c)
    ll_py, ll_at_r = L.tumor_loglik(th, d, 40), L.tumor_loglik(th_r, d, 40)
    comp["objective"] = dict(python_at_python=ll_py, python_at_R=ll_at_r, R_reported=rt["logLik"],
                             ok=abs(ll_py - ll_at_r) < OBJ_TOL and ll_py >= ll_at_r - OBJ_TOL)
    comp["sigma"] = dict(python=float(th["sigma"]), R=float(rt["sigma"]), ok=abs(th["sigma"] - rt["sigma"]) < 1e-3)
    res["tumor"] = comp
    ok &= all(v["ok"] for v in comp.values())

    # Early GLM: coefficients and the profile interval of kappa against R's profile of dD.
    me.d = dict(me.d, y=y_syn)
    for k in ("tau", "sigma_arm"):
        i = me.names.index(k); me.lo[i] = me.hi[i] = 0.0
    fe = FM.multistart(me, np.random.default_rng([FM.SEED, 8]), 4)
    the = me.unpack(np.array(fe["x"]))
    rc = r["early_glm"]["coef"]
    coefs = dict(dD=(the["kappa"], rc["dD"]), dQ=(the["kappa"] / the["ab"], rc["dQ"]),
                 dTe=(-the["kappa"] * the["r"], rc["dTe"]), chemo=(the["gamma"], rc["chemo"]))
    ce = {k: dict(python=float(p_), R=float(r_), ok=close(float(p_), float(r_))) for k, (p_, r_) in coefs.items()}
    prof = FM.profile(me, "kappa", fe, np.random.default_rng(9))
    (piece,) = prof["confidence_set"]
    lo_ci, hi_ci = r["early_glm"]["dD_ci"]
    ce["kappa_profile"] = dict(python=[piece["lower"]["value"], piece["upper"]["value"]], R=[lo_ci, hi_ci],
                               ok=abs(piece["lower"]["value"] - lo_ci) < 1e-3 * abs(lo_ci) + 1e-6
                               and abs(piece["upper"]["value"] - hi_ci) < 1e-3 * abs(hi_ci) + 1e-6)
    ce["truth"] = dict(kappa=SYNTH["kappa"], ab=SYNTH["ab"], r=SYNTH["r"], gamma=SYNTH["gamma"])
    res["early_glm"] = ce
    ok &= all(v["ok"] for k, v in ce.items() if k != "truth")
    res["early_laplace_descriptive"] = r["early_laplace"]
    res["passed"] = bool(ok)
    json.dump(res, open(os.path.join(a.out, "r_check.json"), "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
