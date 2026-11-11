"""Plan step 5: first refit of early (acute mucosal) toxicity to randomized-arm counts.

Usage:
    python3 tools/early_tox/step5_fit.py [--units U02,U03 | --lofo] [--od] [--smoke] [--out DIR]

Decisions applied (Igor Shuryak, 2026-11-11, "C1-C3 as recommended, go"):
    C1  within-trial chemotherapy term gamma ~ N(0, 1) on the log-odds scale, coded per arm against
        its family reference by cumulative planned concurrent dose (+1 more, -1 less, 0 same);
        sensitivity fits drop the five families or use an on/off coding
    C2  Tk fixed at 7 (Fowler), 14 and 21 days as named scenarios
    C3  arm-level overdispersion sigma_arm ~ HalfNormal(0.5). Every unit is fitted twice, binomial
        (OUT/binomial) and overdispersed (--od, OUT/overdispersed); the overdispersed set is primary
        only if the binomial U02 posterior predictive check fails (chi-square p < 0.05)
Model and priors (D3, M1-M4): see tools/early_tox/step5_model.stan and the step 5 record.

Each unit (one fit) is written to OUT/<unit>/ as soon as it finishes and skipped on rerun if its
diagnostics.json exists (plan rule 12). One chain runs at a time, one thread, and the script waits
while Windows free memory is below 3.5 GB.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import csv
import gzip
import importlib.util
import json
import math
import subprocess
import sys
import time

import numpy as np
from cmdstanpy import CmdStanModel, cmdstan_path

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HERE = os.path.join(ROOT, "tools", "early_tox")


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


s4 = _load("step4_identification")
ev = s4.ev

SEED = 20261111
OUT_DEFAULT = os.path.join(ROOT, "docs", "development", "2026_11_10_early_toxicity", "step5")
STAN_FILE = os.path.join(HERE, "step5_model.stan")

# Reference arms of the M4 sensitivity trials (the arm closest to conventional fractionation).
REF_ARM = dict(s4.REF_ARM, IAEA_ACC="IAEA_ACC5", Pinto="Pinto_CF", Poddar="Poddar_5F")

# C1: arm minus its family reference, by cumulative planned concurrent chemotherapy dose.
CHEMO_CUMULATIVE = {
    "GORTEC9902_ACCCT": -1,   # 2 carboplatin/5-FU cycles against 3
    "GORTEC9902_VACC": -1,    # none against 3 cycles
    "RTOG0129_AFXC": -1,      # 2 cisplatin cycles against 3
    "Vienna_VCHART_MMC": +1,  # adds mitomycin C
    "Olmi_CFCT": +1,          # adds carboplatin/5-FU
    "Bartelink_MFD": -1,      # 10 mg/m2 on 15 days (150) against 6 mg/m2 on 35 days (210); daily dose is higher
}
# Sensitivity coding: only a whole regimen added or removed counts.
CHEMO_ONOFF = {"GORTEC9902_VACC": -1, "Vienna_VCHART_MMC": +1, "Olmi_CFCT": +1}
CHEMO_FAMILIES = {"GORTEC9902", "RTOG0129", "Vienna", "ORO9301", "Bartelink2002"}

GROUP_COLS = ("G2", "other")          # G1 is the reference scale group (M1)

# Decision schedules (5 days a week, one fraction a day), each compared with 70 Gy / 35.
DECISION = [(15, 3.5), (20, 2.75), (25, 2.4), (30, 2.2)]

FOWLER = dict(ab=ev.FOWLER_2003["ab"], r=ev.fowler_rate())
PRIOR = dict(Lstar_sd=1.0, ab_mu=math.log(10), ab_sd=0.354, r_mu=math.log(0.8), r_sd=0.5,
             gamma_sd=1.0, od_sd=0.5)

# name: (description, settings)
UNITS = {
    "U01": ("benchmark: Fowler exposure fixed, Tk 7", dict(fixed_E=1)),
    "U02": ("candidate, Tk 7 (primary)", dict()),
    "U03": ("candidate, Tk 14 (C2)", dict(Tk=14.0)),
    "U04": ("candidate, Tk 21 (C2)", dict(Tk=21.0)),
    "U05": ("candidate, L* ~ HalfNormal(0.5)", dict(Lstar_sd=0.5)),
    "U06": ("candidate, L* ~ HalfNormal(2)", dict(Lstar_sd=2.0)),
    "U07": ("candidate, alpha/beta prior sd doubled (0.708)", dict(ab_sd=0.708)),
    "U08": ("candidate, five chemotherapy-differing families dropped, no gamma (C1)", dict(drop_chemo=True, use_gamma=0)),
    "U09": ("candidate, chemotherapy coded on/off only (C1)", dict(chemo="onoff")),
    "U10": ("candidate, W1 windows only (M2)", dict(w1_only=True)),
    "U11": ("candidate, plus IAEA-ACC, Pinto, Poddar (M4)", dict(m4=True)),
    "U14": ("candidate without the chemotherapy term (plan as written)", dict(use_gamma=0)),
}
# Leave-one-family-out influence of the primary candidate (U02 settings), one unit per family.
LOFO_UNITS = {f"L_{f}": (f"candidate without family {f}", dict(drop_family=f)) for f in sorted(s4.REF_ARM)}
REPORT_PARAMS = ["Lstar", "kappa", "ab", "r", "gamma[1]", "sigma_arm[1]", "m0", "tau0"]


def read(name):
    with open(os.path.join(ROOT, "data", "early_tox", name), newline="") as f:
        return list(csv.DictReader(f))


def load_rows(w1_only=False, m4=False, drop_chemo=False, drop_family=None):
    """Arm rows for one analysis set, with D, Q (3.2 h repair) and elapsed-time terms per arm."""
    fr = {}
    for r in read("arm_fractions.csv"):
        fr.setdefault(r["schedule_id"], []).append((int(r["fraction"]), float(r["dose_gy"]), float(r["time_h"])))
    fr = {k: ([d for _, d, _ in sorted(v)], [t for _, _, t in sorted(v)]) for k, v in fr.items()}
    arms = {a["arm_id"]: a for a in read("arms.csv")}
    out = {o["outcome_id"]: o for o in read("outcomes.csv")}
    rows = []
    for m in read("endpoint_mapping.csv"):
        keep = m["in_primary"] == "1" or (m4 and m["analysis_set"] == "sensitivity")
        if not keep or (w1_only and m["in_primary"] == "1" and m["window"] != "W1"):
            continue
        if (drop_chemo and m["family_id"] in CHEMO_FAMILIES) or m["family_id"] == drop_family:
            continue
        o, a = out[m["outcome_id"]], arms[m["arm_id"]]
        variants = []
        for v in a["schedule_variants"].split():
            sid, w = v.split(":")
            variants.append((fr[sid], float(w)))
        rows.append(dict(outcome_id=m["outcome_id"], arm_id=m["arm_id"], family=m["family_id"],
                         scale_group=m["scale_group"], N=int(o["denominator"]), y=int(o["events"]),
                         variants=variants, chemo_text=a["planned_systemic_therapy"]))
    fams = sorted({r["family"] for r in rows})
    for f in fams:
        if sum(1 for r in rows if r["arm_id"] == REF_ARM.get(f)) != 1:
            raise ValueError(f"family {f}: reference arm {REF_ARM.get(f)} not found exactly once")
    return rows, fams


def schedule_terms(doses, times, Tk):
    T = ev.elapsed_days(times)
    return sum(doses), ev.repair_q(doses, times, s4.HALFTIME_H), max(0.0, T - Tk)


def stan_data(rows, fams, Tk, settings):
    terms = np.array([s4.arm_terms(x, Tk)[:3] for x in rows])            # D, Q, Te
    ref = {f: next(i for i, x in enumerate(rows) if x["arm_id"] == REF_ARM[f]) for f in fams}
    d = np.array([terms[i] - terms[ref[x["family"]]] for i, x in enumerate(rows)])
    coding = CHEMO_ONOFF if settings.get("chemo") == "onoff" else CHEMO_CUMULATIVE
    chemo = [float(coding.get(x["arm_id"], 0)) for x in rows]
    grp = [[1.0 if any(x["family"] == f and x["scale_group"] == g for x in rows) else 0.0 for g in GROUP_COLS]
           for f in fams]
    used = [k for k in range(len(GROUP_COLS)) if any(row[k] for row in grp)]
    grp = [[row[k] for k in used] for row in grp]
    ref70 = schedule_terms(*ev.regular_schedule(35, 2.0), Tk)
    dec = np.array([np.subtract(schedule_terms(*ev.regular_schedule(n, dd), Tk), ref70) for n, dd in DECISION])
    pr = dict(PRIOR, **{k: v for k, v in settings.items() if k in PRIOR})
    return dict(
        I=len(rows), J=len(fams), fam=[fams.index(x["family"]) + 1 for x in rows],
        N=[x["N"] for x in rows], y=[x["y"] for x in rows],
        dD=d[:, 0].tolist(), dQ=d[:, 1].tolist(), dTe=d[:, 2].tolist(), chemo=chemo,
        G=len(used), grp=grp if used else np.zeros((len(fams), 0)).tolist(),
        fixed_E=int(settings.get("fixed_E", 0)), ab_fix=FOWLER["ab"], r_fix=FOWLER["r"],
        use_gamma=int(settings.get("use_gamma", 1)), use_od=int(settings.get("use_od", 0)),
        Lstar_sd=pr["Lstar_sd"], ab_mu=pr["ab_mu"], ab_sd=pr["ab_sd"], r_mu=pr["r_mu"], r_sd=pr["r_sd"],
        gamma_sd=pr["gamma_sd"], od_sd=pr["od_sd"],
        K=len(DECISION), kD=dec[:, 0].tolist(), kQ=dec[:, 1].tolist(), kTe=dec[:, 2].tolist(),
        prior_only=int(settings.get("prior_only", 0)))


def windows_free_gb():
    """Free physical memory on the Windows host in GB, or None if the probe gives no number."""
    try:
        out = subprocess.run(["powershell.exe", "-NoProfile", "-Command",
                              "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"],
                             capture_output=True, text=True, timeout=60).stdout.strip()
        return int(out) / 1048576 if out.isdigit() else None
    except (OSError, subprocess.SubprocessError, ValueError):
        return None


def wait_for_memory(limit_gb=3.5, max_wait_s=3600):
    t0 = time.time()
    while True:
        free = windows_free_gb()
        if free is not None and free >= limit_gb:
            return free
        if time.time() - t0 > max_wait_s:
            raise RuntimeError(f"Windows free memory {free} GB stayed below {limit_gb} GB (or probe failed)")
        print(f"  waiting: Windows free memory {free} GB", flush=True)
        time.sleep(60)


def compile_model():
    return CmdStanModel(stan_file=STAN_FILE)


def diagnostics(fit, names):
    s = fit.summary(percentiles=(5, 50, 95))
    rows = [n for n in s.index if n in names or n.split("[")[0] in ("dE_dec", "lor_dec", "eta", "delta")]
    s = s.loc[rows]
    div = int(np.sum(fit.divergences))
    worst = dict(rhat=float(s["R_hat"].max()), ess_bulk=float(s["ESS_bulk"].min()), ess_tail=float(s["ESS_tail"].min()))
    gate = div == 0 and worst["rhat"] < 1.01 and worst["ess_bulk"] > 400 and worst["ess_tail"] > 400
    return s, dict(divergences=div, **worst, gate_pass=bool(gate))


def run_unit(model, name, out_dir, smoke=False, od=False):
    desc, settings = {**UNITS, **LOFO_UNITS}[name]
    settings = dict(settings, use_od=1) if od else settings
    udir = os.path.join(out_dir, name)
    if os.path.exists(os.path.join(udir, "diagnostics.json")):
        print(f"{name}: done, skipped", flush=True)
        return
    os.makedirs(udir, exist_ok=True)
    rows, fams = load_rows(settings.get("w1_only", False), settings.get("m4", False), settings.get("drop_chemo", False),
                           settings.get("drop_family"))
    Tk = settings.get("Tk", 7.0)
    data = stan_data(rows, fams, Tk, settings)
    free = wait_for_memory()
    t0 = time.time()
    fit = model.sample(data=data, chains=4, parallel_chains=1, threads_per_chain=1,
                       iter_warmup=200 if smoke else 1000, iter_sampling=200 if smoke else 1000,
                       adapt_delta=0.95, seed=SEED, show_progress=False, show_console=False)
    secs = time.time() - t0
    s, diag = diagnostics(fit, set(REPORT_PARAMS))
    s.to_csv(os.path.join(udir, "summary.csv"), float_format="%.5g")
    dr = fit.draws_pd()
    keep = [c for c in dr.columns if c in REPORT_PARAMS or c.split("[")[0] in ("dE_dec", "lor_dec", "T_obs", "T_rep")]
    with gzip.open(os.path.join(udir, "draws.csv.gz"), "wt") as f:
        dr[keep].to_csv(f, index=False, float_format="%.6g")
    yrep = fit.stan_variable("y_rep")
    with open(os.path.join(udir, "ppc_arms.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["arm_id", "family", "N", "y", "chemo", "yrep_median", "yrep_5", "yrep_95", "p_upper", "p_lower"])
        for i, x in enumerate(rows):
            col = yrep[:, i]
            w.writerow([x["arm_id"], x["family"], x["N"], x["y"], data["chemo"][i],
                        np.median(col), np.percentile(col, 5), np.percentile(col, 95),
                        round(float(np.mean(col >= x["y"])), 4), round(float(np.mean(col <= x["y"])), 4)])
    T_obs, T_rep = fit.stan_variable("T_obs"), fit.stan_variable("T_rep")
    ll = fit.stan_variable("log_lik")
    lppd = float(np.sum(np.log(np.mean(np.exp(ll), axis=0))))
    p_waic = float(np.sum(np.var(ll, axis=0, ddof=1)))
    diag.update(unit=name, description=desc, settings=settings, Tk=Tk, n_arms=len(rows), n_families=len(fams),
                families=fams, smoke=smoke, seconds=round(secs, 1), windows_free_gb_at_start=round(free, 1),
                ppc_chisq_p=round(float(np.mean(T_rep >= T_obs)), 4),
                waic=round(-2 * (lppd - p_waic), 3), p_waic=round(p_waic, 3),
                decision_schedules=[f"{n} x {d} Gy" for n, d in DECISION],
                cmdstan=os.path.basename(cmdstan_path()),
                revision=subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"],
                                        capture_output=True, text=True).stdout.strip())
    with open(os.path.join(udir, "diagnostics.json"), "w") as f:
        json.dump(diag, f, indent=1, default=str)
    print(f"{name}: {desc}: {secs:.0f} s, divergences {diag['divergences']}, Rhat {diag['rhat']:.4f}, "
          f"ESS {diag['ess_bulk']:.0f}/{diag['ess_tail']:.0f}, PPC p {diag['ppc_chisq_p']}, "
          f"gate {'PASS' if diag['gate_pass'] else 'FAIL'}", flush=True)


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--units", default=",".join(UNITS))
    ap.add_argument("--lofo", action="store_true", help="run the leave-one-family-out units instead")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--od", action="store_true", help="add arm-level overdispersion to every unit (C3)")
    ap.add_argument("--out", default=OUT_DEFAULT)
    a = ap.parse_args(argv[1:])
    out = os.path.join(a.out + ("_smoke" if a.smoke else ""), "overdispersed" if a.od else "binomial")
    os.makedirs(out, exist_ok=True)
    model = compile_model()
    for u in (list(LOFO_UNITS) if a.lofo else a.units.split(",")):
        run_unit(model, u, out, smoke=a.smoke, od=a.od)


if __name__ == "__main__":
    main(sys.argv)
