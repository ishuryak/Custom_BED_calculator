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

Each unit (one fit) is written to OUT/<unit>.partial-<pid>/ and published to OUT/<unit>/ by an
atomic rename after its manifest.json (run identity, output hashes, gate, raw-chain archive) is
written. A rerun reuses a unit only when its identity and output hashes match, and refuses a
mismatch, a missing or unreadable manifest, or changed outputs (review F4). Raw CmdStan chains
are archived under $EARLY_TOX_ARCHIVE (default ~/early_tox_chain_archive), outside the repository. One chain runs at a time, one thread, and the script waits
while Windows free memory is below 3.5 GB.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import shutil
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
REF_ARM = dict(s4.REF_ARM, IAEA_ACC="IAEA_ACC5", Pinto="Pinto_CF", Poddar="Poddar_5F",
               Sanchiz="Sanchiz_A", Dwivedi="Dwivedi_CF", ELAN="ELAN_SF")
SENSITIVITY_FAMILIES = {"IAEA_ACC", "Pinto", "Poddar", "Sanchiz", "Dwivedi", "ELAN"}

# C1: arm minus its family reference, by cumulative planned concurrent chemotherapy dose.
CHEMO_CUMULATIVE = {
    "GORTEC9902_ACCCT": -1,   # 2 carboplatin/5-FU cycles against 3
    "GORTEC9902_VACC": -1,    # none against 3 cycles
    "RTOG0129_AFXC": -1,      # 2 cisplatin cycles against 3
    "Vienna_VCHART_MMC": +1,  # adds mitomycin C
    "Olmi_CFCT": +1,          # adds carboplatin/5-FU
    "Bartelink_MFD": -1,      # 10 mg/m2 on 15 days (150) against 6 mg/m2 on 35 days (210); daily dose is higher
    "GhoshLaskar_CRT": +1,    # adds weekly cisplatin
    "GuptaP_CBRT": -1,        # none against 3 x cisplatin 75 mg/m2
    "Sanchiz_C": +1,          # adds 5-FU
}
# Arms whose planned therapy text differs from the reference arm's but whose cumulative planned
# concurrent dose is the same (code 0).
CHEMO_SAME_DOSE = {"GuptaP_CBRTCCT": "cisplatin 3 x 75 mg/m2 in both arms; only the days differ"}
# Sensitivity coding: only a whole regimen added or removed counts.
# External (never fitted) arms whose planned chemotherapy differs from their reference arm: scored
# under declared scenarios, never a single silent code (review F1, 2026-11-11). HYPNO planned 4 vs 5
# weekly cisplatin courses, prescribed in about 75% of patients.
EXTERNAL_CHEMO_SCENARIOS = {
    "HYPNO_HFX": {"code_0": 0.0, "code_minus1": -1.0, "mixture_75pct": ("mixture", 0.75, -1.0)},
}
# Decision D5 sensitivity: Tandon's 25-fraction arm coded -1 (planned courses not stated).
CHEMO_SENSITIVITY = {"Tandon_SMART": -1}


def planned_therapy():
    """Structured planned systemic therapy per arm (data/early_tox/source/planned_therapy.csv)."""
    out = {}
    for r in read_source("planned_therapy.csv"):
        if r["arm_id"] in out:
            raise ValueError(f"{r['arm_id']}: duplicate planned-therapy row")
        out[r["arm_id"]] = r
    return out


def read_source(name):
    with open(os.path.join(ROOT, "data", "early_tox", "source", name), newline="") as f:
        return list(csv.DictReader(f))


CHEMO_ONOFF = {"GORTEC9902_VACC": -1, "Vienna_VCHART_MMC": +1, "Olmi_CFCT": +1,
               "GhoshLaskar_CRT": +1, "GuptaP_CBRT": -1, "Sanchiz_C": +1}
CHEMO_FAMILIES = {"GORTEC9902", "RTOG0129", "Vienna", "ORO9301", "Bartelink2002", "GhoshLaskar", "GuptaP",
                  "Sanchiz"}

GROUP_COLS = ("G2", "other")          # G1 is the reference scale group (M1)

# Decision schedules (5 days a week, one fraction a day), each compared with 70 Gy / 35.
DECISION = [(15, 3.5), (20, 2.75), (25, 2.4), (30, 2.2)]

FOWLER = dict(ab=ev.FOWLER_2003["ab"], r=ev.fowler_rate())
PRIOR = dict(Lstar_sd=1.0, ab_mu=math.log(10), ab_sd=0.354, r_mu=math.log(0.8), r_sd=0.5,
             gamma_sd=1.0, od_sd=0.5)

# name: (description, settings)
COMPLETE_REPAIR = "complete"
UNITS = {
    "U01": ("benchmark: Fowler exposure fixed, Tk 7", dict(fixed_E=1)),
    "U02": ("candidate, Tk 7 (primary)", dict()),
    "U03": ("candidate, Tk 14 (C2)", dict(Tk=14.0)),
    "U04": ("candidate, Tk 21 (C2)", dict(Tk=21.0)),
    "U05": ("candidate, L* ~ HalfNormal(0.5)", dict(Lstar_sd=0.5)),
    "U06": ("candidate, L* ~ HalfNormal(2)", dict(Lstar_sd=2.0)),
    "U07": ("candidate, alpha/beta prior sd doubled (0.708)", dict(ab_sd=0.708)),
    "U08": ("candidate, chemotherapy-differing families dropped, no gamma (C1)", dict(drop_chemo=True, use_gamma=0)),
    "U09": ("candidate, chemotherapy coded on/off only (C1)", dict(chemo="onoff")),
    "U10": ("candidate, nominal W1 windows only (M2; stated window, not visit schedule)", dict(w1_only=True)),
    "U11": ("candidate, plus the sensitivity trials (M4 and 2026-11-11 additions)", dict(m4=True)),
    "U14": ("candidate without the chemotherapy term (plan as written)", dict(use_gamma=0)),
    "U15": ("candidate, adequate peak ascertainment only (decision D3)", dict(adequate_only=True)),
    "U16": ("benchmark, adequate peak ascertainment only (decision D3)", dict(adequate_only=True, fixed_E=1)),
    "U17": ("candidate, repair half-time 2 h (F8)", dict(halftime_h=2.0)),
    "U18": ("candidate, repair half-time 4 h (F8)", dict(halftime_h=4.0)),
    "U19": ("candidate, complete repair (F8)", dict(halftime_h=COMPLETE_REPAIR)),
    "U20": ("benchmark, repair half-time 2 h (F8)", dict(halftime_h=2.0, fixed_E=1)),
    "U21": ("benchmark, repair half-time 4 h (F8)", dict(halftime_h=4.0, fixed_E=1)),
    "U22": ("benchmark, complete repair (F8)", dict(halftime_h=COMPLETE_REPAIR, fixed_E=1)),
}
# Leave-one-family-out influence of the primary candidate (U02 settings), one unit per family.
LOFO_UNITS = {f"L_{f}": (f"candidate without family {f}", dict(drop_family=f))
              for f in sorted(REF_ARM) if f not in SENSITIVITY_FAMILIES}
REPORT_PARAMS = ["Lstar", "kappa", "ab", "r", "gamma[1]", "sigma_arm[1]", "m0", "tau0"]


def read(name):
    with open(os.path.join(ROOT, "data", "early_tox", name), newline="") as f:
        return list(csv.DictReader(f))


def load_rows(w1_only=False, m4=False, drop_chemo=False, drop_family=None, adequate_only=False):
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
        if adequate_only and m["in_adequate_sensitivity"] != "1":
            continue
        o, a = out[m["outcome_id"]], arms[m["arm_id"]]
        variants = []
        for v in a["schedule_variants"].split():
            sid, w = v.split(":")
            variants.append((fr[sid], float(w)))
        rows.append(dict(outcome_id=m["outcome_id"], arm_id=m["arm_id"], family=m["family_id"],
                         scale_group=m["scale_group"], N=int(o["denominator"]), y=int(o["events"]),
                         variants=variants, chemo_text=a["planned_systemic_therapy"]))
    ids = [r["arm_id"] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError(f"arm(s) entered twice: {sorted({a for a in ids if ids.count(a) > 1})}")
    fams = sorted({r["family"] for r in rows})
    for f in fams:
        if sum(1 for r in rows if r["arm_id"] == REF_ARM.get(f)) != 1:
            raise ValueError(f"family {f}: reference arm {REF_ARM.get(f)} not found exactly once")
    return rows, fams


def analysis_totals(rows):
    """Totals of an analysis set, derived from the loaded rows (review F6: never typed by hand)."""
    ids = [r["arm_id"] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("an arm is counted twice")
    return dict(arms=len(rows), families=len({r["family"] for r in rows}),
                evaluable=sum(r["N"] for r in rows), events=sum(r["y"] for r in rows))


def halftime_of(settings):
    """Repair half-time in hours from the settings: 3.2 by default, None for complete repair."""
    h = settings.get("halftime_h", s4.HALFTIME_H)
    return None if h == COMPLETE_REPAIR else float(h)


def schedule_terms(doses, times, Tk, halftime_h=s4.HALFTIME_H):
    T = ev.elapsed_days(times)
    return sum(doses), ev.repair_q(doses, times, halftime_h), max(0.0, T - Tk)


def arm_terms(row, Tk, halftime_h=s4.HALFTIME_H):
    """Patient-weighted D, Q and T_eff = max(0, T - Tk) over an arm's schedule variants."""
    tot = sum(w for _, w in row["variants"])
    acc = np.zeros(3)
    for (doses, times), w in row["variants"]:
        acc += w * np.array(schedule_terms(doses, times, Tk, halftime_h))
    return tuple(acc / tot)


def stan_data(rows, fams, Tk, settings):
    h = halftime_of(settings)
    terms = np.array([arm_terms(x, Tk, h) for x in rows])            # D, Q, Te
    ref = {f: next(i for i, x in enumerate(rows) if x["arm_id"] == REF_ARM[f]) for f in fams}
    d = np.array([terms[i] - terms[ref[x["family"]]] for i, x in enumerate(rows)])
    coding = CHEMO_ONOFF if settings.get("chemo") == "onoff" else CHEMO_CUMULATIVE
    chemo = [float(coding.get(x["arm_id"], 0)) for x in rows]
    grp = [[1.0 if any(x["family"] == f and x["scale_group"] == g for x in rows) else 0.0 for g in GROUP_COLS]
           for f in fams]
    used = [k for k in range(len(GROUP_COLS)) if any(row[k] for row in grp)]
    grp = [[row[k] for k in used] for row in grp]
    ref70 = schedule_terms(*ev.regular_schedule(35, 2.0), Tk, h)
    dec = np.array([np.subtract(schedule_terms(*ev.regular_schedule(n, dd), Tk, h), ref70) for n, dd in DECISION])
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


SAMPLER = dict(chains=4, parallel_chains=1, threads_per_chain=1, iter_warmup=1000, iter_sampling=1000,
               adapt_delta=0.95, seed=SEED)
SMOKE_SAMPLER = dict(SAMPLER, iter_warmup=200, iter_sampling=200)
EBFMI_MIN = 0.3
ARCHIVE_ROOT = os.environ.get("EARLY_TOX_ARCHIVE", os.path.expanduser("~/early_tox_chain_archive"))
OUTPUT_FILES = ("summary.csv", "ppc_arms.csv", "diagnostics.json")


def _normalize(x):
    if isinstance(x, float):
        return float(f"{x:.12g}")
    if isinstance(x, dict):
        return {k: _normalize(v) for k, v in sorted(x.items())}
    if isinstance(x, (list, tuple)):
        return [_normalize(v) for v in x]
    if isinstance(x, np.ndarray):
        return _normalize(x.tolist())
    if isinstance(x, (np.floating, np.integer)):
        return _normalize(x.item())
    return x


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_identity(data, settings, sampler, extra=None):
    """One hash over everything that defines a fit (review F4): data in order, settings (repair, Tk,
    priors), Stan source, sampler configuration and tool versions."""
    import cmdstanpy
    payload = dict(data=_normalize(data), settings=_normalize(settings), sampler=_normalize(sampler),
                   stan_sha256=sha256_file(STAN_FILE), cmdstan=os.path.basename(cmdstan_path()),
                   cmdstanpy=cmdstanpy.__version__, numpy=np.__version__, extra=_normalize(extra))
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


def required_quantities(data):
    """Every sampled or decision quantity whose diagnostics the gate must see (review F5)."""
    J, I, G, K = data["J"], data["I"], data["G"], data["K"]
    req = ["m0", "tau0", "Lstar", "kappa"] + [f"z[{j}]" for j in range(1, J + 1)]
    req += [f"eta[{j}]" for j in range(1, J + 1)] + [f"delta[{g}]" for g in range(1, G + 1)]
    req += [f"lor_dec[{k}]" for k in range(1, K + 1)]
    if not data["fixed_E"]:
        req += ["log_ab[1]", "log_r[1]", "ab", "r"] + [f"dE_dec[{k}]" for k in range(1, K + 1)]
    if data["use_gamma"]:
        req += ["gamma[1]"]
    if data["use_od"]:
        req += ["sigma_arm[1]"] + [f"zarm[{i}]" for i in range(1, I + 1)]
    return req


def gate_from_summary(s, data, divergences, ebfmi, treedepth_hits):
    """Exact-schema gate: every required quantity present with finite R-hat and ESS; fixed
    benchmark constants checked as constants; no divergences; E-BFMI >= 0.3 in every chain."""
    req = required_quantities(data)
    missing = [q for q in req if q not in s.index]
    present = [q for q in req if q in s.index]
    vals = s.loc[present, ["R_hat", "ESS_bulk", "ESS_tail"]].astype(float)
    nonfinite = [q for q in present if not np.isfinite(vals.loc[q]).all()]
    fixed_ok = True
    if data["fixed_E"]:
        for q, v in (("ab", data["ab_fix"]), ("r", data["r_fix"])):
            fixed_ok &= q in s.index and abs(float(s.loc[q, "Mean"]) - v) < 1e-9 and float(s.loc[q, "StdDev"]) < 1e-9
    finite = vals.loc[[q for q in present if q not in nonfinite]]
    out = dict(divergences=int(divergences), missing=missing, nonfinite=nonfinite, fixed_constants_ok=bool(fixed_ok),
               rhat=float(finite["R_hat"].max()) if len(finite) else float("nan"),
               ess_bulk=float(finite["ESS_bulk"].min()) if len(finite) else float("nan"),
               ess_tail=float(finite["ESS_tail"].min()) if len(finite) else float("nan"),
               ebfmi_min=float(min(ebfmi)) if len(ebfmi) else float("nan"),
               treedepth_saturated=int(treedepth_hits), n_required=len(req))
    out["gate_pass"] = bool(not missing and not nonfinite and fixed_ok and out["divergences"] == 0
                            and out["rhat"] < 1.01 and out["ess_bulk"] > 400 and out["ess_tail"] > 400
                            and np.isfinite(out["ebfmi_min"]) and out["ebfmi_min"] >= EBFMI_MIN)
    return out


def chain_energy_checks(fit, max_treedepth=10):
    mv = fit.method_variables()
    energy, depth = mv["energy__"], mv["treedepth__"]                  # (draws, chains)
    ebfmi = [float(np.sum(np.diff(energy[:, c]) ** 2) / len(energy[:, c]) / np.var(energy[:, c]))
             for c in range(energy.shape[1])]
    return ebfmi, int(np.sum(depth >= max_treedepth))


def archive_chains(fit, identity):
    d = os.path.join(ARCHIVE_ROOT, identity[:16])
    os.makedirs(d, exist_ok=True)
    fit.save_csvfiles(d)
    return {x: sha256_file(os.path.join(d, x)) for x in sorted(os.listdir(d)) if x.endswith(".csv")}


def unit_state(udir, identity):
    """'complete' (reusable), 'failed' (same identity, gate failed), 'absent', or raise on mismatch."""
    man = os.path.join(udir, "manifest.json")
    if not os.path.exists(udir):
        return "absent"
    if not os.path.exists(man):
        raise RuntimeError(f"{udir}: outputs without a manifest (legacy or interrupted); refusing to reuse")
    try:
        m = json.load(open(man))
    except (OSError, ValueError) as err:
        raise RuntimeError(f"{udir}: unreadable manifest ({err}); refusing to reuse")
    if m.get("identity") != identity:
        raise RuntimeError(f"{udir}: run identity differs from the current inputs; refusing to reuse")
    for fn, digest in m.get("outputs", {}).items():
        path = os.path.join(udir, fn)
        if not os.path.exists(path) or sha256_file(path) != digest:
            raise RuntimeError(f"{udir}: output {fn} missing or changed since completion; refusing to reuse")
    return "complete" if m.get("gate_pass") else "failed"


def run_unit(model, name, out_dir, smoke=False, od=False):
    desc, settings = {**UNITS, **LOFO_UNITS}[name]
    settings = dict(settings, use_od=1) if od else settings
    udir = os.path.join(out_dir, name)
    rows, fams = load_rows(settings.get("w1_only", False), settings.get("m4", False), settings.get("drop_chemo", False),
                           settings.get("drop_family"), settings.get("adequate_only", False))
    Tk = settings.get("Tk", 7.0)
    data = stan_data(rows, fams, Tk, settings)
    sampler = SMOKE_SAMPLER if smoke else SAMPLER
    identity = run_identity(data, settings, sampler, extra=dict(unit=name, families=fams))
    state = unit_state(udir, identity)
    if state != "absent":
        print(f"{name}: {state} with identical identity, not rerun", flush=True)
        return
    tmp = f"{udir}.partial-{os.getpid()}"
    if os.path.exists(tmp):
        shutil.rmtree(tmp)
    os.makedirs(tmp)
    free = wait_for_memory()
    t0 = time.time()
    fit = model.sample(data=data, show_progress=False, show_console=False, **sampler)
    secs = time.time() - t0
    s = fit.summary(percentiles=(5, 50, 95))
    ebfmi, depth_hits = chain_energy_checks(fit)
    diag = gate_from_summary(s, data, np.sum(fit.divergences), ebfmi, depth_hits)
    keep_rows = [n for n in s.index if n in REPORT_PARAMS or n.split("[")[0] in ("dE_dec", "lor_dec", "eta", "delta")]
    s.loc[keep_rows].to_csv(os.path.join(tmp, "summary.csv"), float_format="%.5g")
    dr = fit.draws_pd()
    keep = [c for c in dr.columns if c in REPORT_PARAMS or c.split("[")[0] in ("dE_dec", "lor_dec", "T_obs", "T_rep")]
    with gzip.open(os.path.join(tmp, "draws.csv.gz"), "wt") as f:
        dr[keep].to_csv(f, index=False, float_format="%.6g")
    yrep = fit.stan_variable("y_rep")
    with open(os.path.join(tmp, "ppc_arms.csv"), "w", newline="") as f:
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
    diag.update(unit=name, description=desc, settings=settings, Tk=Tk, halftime_h=halftime_of(settings),
                n_arms=len(rows), n_families=len(fams), families=fams, totals=analysis_totals(rows), smoke=smoke,
                sampler=sampler, seconds=round(secs, 1), windows_free_gb_at_start=round(free, 1),
                ppc_chisq_p=round(float(np.mean(T_rep >= T_obs)), 4),
                waic=round(-2 * (lppd - p_waic), 3), p_waic=round(p_waic, 3),
                decision_schedules=[f"{n} x {d} Gy" for n, d in DECISION],
                cmdstan=os.path.basename(cmdstan_path()), identity=identity,
                revision=subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"],
                                        capture_output=True, text=True).stdout.strip())
    with open(os.path.join(tmp, "diagnostics.json"), "w") as f:
        json.dump(diag, f, indent=1, default=str)
    chains = archive_chains(fit, identity)
    manifest = dict(identity=identity, unit=name, gate_pass=diag["gate_pass"], smoke=smoke,
                    outputs={fn: sha256_file(os.path.join(tmp, fn)) for fn in OUTPUT_FILES + ("draws.csv.gz",)},
                    chain_archive=os.path.join(ARCHIVE_ROOT, identity[:16]), chain_csv_sha256=chains)
    with open(os.path.join(tmp, "manifest.json.tmp"), "w") as f:
        json.dump(manifest, f, indent=1)
    os.replace(os.path.join(tmp, "manifest.json.tmp"), os.path.join(tmp, "manifest.json"))   # manifest last
    os.replace(tmp, udir)                                                                    # then publish
    print(f"{name}: {desc}: {secs:.0f} s, divergences {diag['divergences']}, Rhat {diag['rhat']:.4f}, "
          f"ESS {diag['ess_bulk']:.0f}/{diag['ess_tail']:.0f}, E-BFMI {diag['ebfmi_min']:.2f}, "
          f"missing {len(diag['missing'])}, PPC p {diag['ppc_chisq_p']}, "
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
