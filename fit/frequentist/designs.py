"""Designs of the two Bayesian fits, imported from the builders that produced them (never rewritten),
and the three identity layers of review F02:

1. input files: sha256 of the bytes of every data file the design is read from;
2. canonical design: sha256 of the ordered arrays and the ordered unit/family and arm IDs that
   enter the likelihood (independent of file formatting, sensitive to any count, order or time);
3. run identity (fit_ml.py): design digest plus model, support, integration and optimizer policy,
   seeds, code hashes and software versions.

The Bayesian fits' own identities are reproduced through tools/build_model_methods.py, which refuses
unless the rebuilt Stan data hash (tumor) and the U02 run identity (early) match the retained fits.
"""
import hashlib
import importlib.util
import json
import os
import sys
import types

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TUMOR_DATA = os.path.join(ROOT, "data", "refit_trials_2026_10_25c_main.csv")
EARLY_DATA_FILES = ["arm_fractions.csv", "arms.csv", "outcomes.csv", "endpoint_mapping.csv"]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha_file(p):
    with open(p, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def canonical_digest(d, keys, ids):
    payload = {k: np.asarray(d[k]).tolist() for k in keys}
    payload["ids"] = ids
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=float).encode()).hexdigest()


TUMOR_KEYS = ["unit", "dsum", "d2sum", "T", "npat", "y", "eval_dsum", "eval_d2sum", "eval_T",
              "anchor_dsum", "anchor_d2sum", "anchor_T"]
EARLY_KEYS = ["fam", "N", "y", "dD", "dQ", "dTe", "chemo", "grp", "kD", "kQ", "kTe"]


def tumor_design(path=TUMOR_DATA):
    F = _load("fit_refit_fq", os.path.join(ROOT, "fit", "fit_refit.py"))
    rows, trials, base = F.load(path)
    d = dict(N=base["N"], J=base["J"], unit=np.array(base["trial"]) - 1,
             dsum=np.array(base["dsum"], float), d2sum=np.array(base["d2sum"], float), T=np.array(base["T"], float),
             npat=np.array(base["npat"]), y=np.array(base["y"]),
             eval_dsum=np.array(base["eval_dsum"], float), eval_d2sum=np.array(base["eval_d2sum"], float),
             eval_T=np.array(base["eval_T"], float), anchor_dsum=float(base["anchor_dsum"]),
             anchor_d2sum=float(base["anchor_d2sum"]), anchor_T=float(base["anchor_T"]),
             regimens=[r[0] for r in F.REGIMENS], units=list(trials), arms=[r["row"] for r in rows])
    d["file_sha256"] = {os.path.basename(path): sha_file(path)}
    d["design_sha256"] = canonical_digest(d, TUMOR_KEYS, dict(units=d["units"], arms=d["arms"]))
    d["stan_base"] = base
    return d


def _step5():
    if "cmdstanpy" not in sys.modules:
        try:
            import cmdstanpy  # noqa: F401
        except ImportError:          # only the data builder is used here; nothing is sampled
            stub = types.ModuleType("cmdstanpy")
            stub.CmdStanModel, stub.cmdstan_path = object, lambda: ""
            sys.modules["cmdstanpy"] = stub
    sys.path.insert(0, os.path.join(ROOT, "tools", "early_tox"))
    return _load("step5_fit_fq", os.path.join(ROOT, "tools", "early_tox", "step5_fit.py"))


def early_design():
    S = _step5()
    unit = os.path.join(ROOT, "docs", "development", "2026_11_10_early_toxicity", "step5_review", "overdispersed", "U02")
    diag = json.load(open(os.path.join(unit, "diagnostics.json")))
    settings = dict(S.UNITS["U02"][1], **diag["settings"])
    rows, fams = S.load_rows(False, False, False, None, False)
    sd = S.stan_data(rows, fams, diag["Tk"], settings)
    if not (sd["use_od"] == 1 and sd["use_gamma"] == 1 and sd["fixed_E"] == 0):
        raise RuntimeError("U02 settings are not the primary candidate (overdispersion, gamma, free exposure)")
    d = dict(I=sd["I"], J=sd["J"], G=sd["G"], fam=np.array(sd["fam"]) - 1, N=np.array(sd["N"]), y=np.array(sd["y"]),
             dD=np.array(sd["dD"]), dQ=np.array(sd["dQ"]), dTe=np.array(sd["dTe"]), chemo=np.array(sd["chemo"], float),
             grp=np.array(sd["grp"], float).reshape(sd["J"], sd["G"]), kD=np.array(sd["kD"]), kQ=np.array(sd["kQ"]),
             kTe=np.array(sd["kTe"]), families=list(fams), arms=[x["arm_id"] for x in rows],
             decision=[f"{n} x {dd:g} Gy vs 70 Gy / 35" for n, dd in S.DECISION], Tk=diag["Tk"],
             halftime_h=diag["halftime_h"])
    d["file_sha256"] = {f: sha_file(os.path.join(ROOT, "data", "early_tox", f)) for f in EARLY_DATA_FILES}
    d["design_sha256"] = canonical_digest(d, EARLY_KEYS, dict(families=d["families"], arms=d["arms"]))
    d["stan_data"] = sd
    return d


def reproduce_bayesian_identities():
    """Raise unless the retained Bayesian fits' identities reproduce from the current data and code
    (the checks of tools/build_model_methods.py). Returns the two reproduced identities."""
    bm = _load("build_model_methods_fq", os.path.join(ROOT, "tools", "build_model_methods.py"))
    try:
        t, e = bm.tumor(ROOT), bm.early(ROOT)
    except SystemExit as err:
        raise RuntimeError(str(err)) from None
    return dict(tumor_di=t["di"]["identity"], tumor_stan_data=t["di"]["stan_data_sha256"], early=e["identity"])
