"""Recompute the complete-schema gate of every retained fit from its archived raw CmdStan chains.

Usage:
    python3 tools/early_tox/recheck_chains.py [--out CSV]

Review R6 (v4.6): per-chain E-BFMI, divergences, treedepth and every required quantity's R-hat and
ESS are recomputed from the raw chain CSVs (not from stored aggregates), for the step 5 review
units (manifest chain_archive) and the step 6 review caches (chain_archive = ARCHIVE/identity[:16]).
Each archive's CSV hashes are checked against the hashes recorded when it was written.
Writes one row per fit; exit status 1 if any fit fails or cannot be checked.
"""
import csv
import glob
import importlib.util
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_spec = importlib.util.spec_from_file_location("step6_validate", os.path.join(ROOT, "tools", "early_tox", "step6_validate.py"))
s6 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s6)
s5 = s6.s5
BASE = os.path.join(ROOT, "docs", "development", "2026_11_10_early_toxicity")
OUT_DEFAULT = os.path.join(BASE, "chain_recheck_2026_11_11.csv")


def gate_from_chains(archive, expected_hashes, data):
    from cmdstanpy import from_csv
    files = sorted(f for f in os.listdir(archive) if f.endswith(".csv"))
    got = {f: s5.sha256_file(os.path.join(archive, f)) for f in files}
    if expected_hashes is not None and got != expected_hashes:
        return dict(status="archive hashes differ from those recorded")
    fit = from_csv([os.path.join(archive, f) for f in files])
    summ = fit.summary(percentiles=(5, 50, 95))
    ebfmi, depth = s5.chain_energy_checks(fit)
    g = s5.gate_from_summary(summ, data, int(np.sum(fit.divergences)), ebfmi, depth, n_chains=fit.chains)
    g["status"] = "checked"
    return g


def step5_units():
    for man_path in sorted(glob.glob(os.path.join(BASE, "step5_review", "*", "*", "manifest.json"))):
        udir = os.path.dirname(man_path)
        man = json.load(open(man_path))
        diag = json.load(open(os.path.join(udir, "diagnostics.json")))
        st = diag["settings"]
        rows, fams = s5.load_rows(st.get("w1_only", False), st.get("m4", False), st.get("drop_chemo", False),
                                  st.get("drop_family"), st.get("adequate_only", False))
        data = s5.stan_data(rows, fams, diag["Tk"], st)
        yield os.path.relpath(udir, BASE), man["chain_archive"], man["chain_csv_sha256"], data, man["identity"]


def step6_caches():
    prim, fams = s5.load_rows()
    cache = os.path.join(BASE, "step6_review", "cache")
    for path in sorted(glob.glob(os.path.join(cache, "*.npz"))):
        z = np.load(path, allow_pickle=False)
        name = os.path.basename(path)[:-4]
        model, rest = name.split("__", 1)
        st = dict(s6.MODELS[model])
        if "__h" in rest:
            h = rest.split("__h")[1]
            st["halftime_h"] = h if h == s5.COMPLETE_REPAIR else float(h)
        if rest.startswith("without_"):
            rows, tf = s6.training_set(prim, fams, rest[len("without_"):])
        else:
            rows, tf = prim, fams
        data = s5.stan_data(rows, tf, s6.TK, st)
        ident = str(z["identity"])
        yield os.path.relpath(path, BASE), os.path.join(s5.ARCHIVE_ROOT, ident[:16]), \
            json.loads(str(z["chain_archive"])), data, ident


def main(argv):
    out = argv[1] if len(argv) > 1 else OUT_DEFAULT
    rows, bad = [], 0
    for source in (step5_units(), step6_caches()):
        for name, archive, hashes, data, ident in source:
            try:
                g = gate_from_chains(archive, hashes, data)
            except Exception as err:                     # unreadable archive: recorded, counted as failure
                g = dict(status=f"not checked: {err}")
            ok = g.get("status") == "checked" and g.get("gate_pass") is True
            bad += not ok
            rows.append(dict(fit=name, identity=ident[:16], status=g.get("status"), gate_pass=g.get("gate_pass"),
                             divergences=g.get("divergences"), rhat=g.get("rhat"), ess_bulk=g.get("ess_bulk"),
                             ess_tail=g.get("ess_tail"), ebfmi_per_chain=g.get("ebfmi_per_chain"),
                             treedepth_saturated=g.get("treedepth_saturated"), missing=len(g.get("missing", []) or []),
                             nonfinite=len(g.get("nonfinite", []) or [])))
            print(f"{name}: {g.get('status')} gate {g.get('gate_pass')}", flush=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} fits rechecked from raw chains; {bad} failed or not checked")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
