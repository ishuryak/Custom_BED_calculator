"""Write the draws files the app reads (draws/) and their manifest, from the
fits in fit/output/. Each file holds 4000 draws (every second post-warmup draw
of the 8000), in the units the app uses: alpha (1/Gy), ab (Gy), lambda (1/day),
Tk (inclusive days, DI) or c_over_alpha (Gy, DD), g (1/day).

Run from the repository root after fit/fit_refit.py:  python3 fit/export_draws.py
"""
import hashlib, json, os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "fit", "output")
DRAWS = os.path.join(ROOT, "draws")
COLS = {"di": ["alpha", "ab", "lambda", "Tk", "g"], "dd": ["alpha", "ab", "lambda", "c_over_alpha", "g"]}
SHIP = ["primary", "paper_restrictions"]
N_KEEP = 4000

os.makedirs(DRAWS, exist_ok=True)
prov = json.load(open(os.path.join(OUT, "provenance.json")))
man = dict(created_by="fit/export_draws.py", data_sha256=prov["data_sha256"], fit=prov, files={})
for model, cols in COLS.items():
    for variant in SHIP:
        s = json.load(open(os.path.join(OUT, f"{model}_{variant}", "summary.json")))
        dr = pd.read_csv(os.path.join(OUT, f"{model}_{variant}", "draws_params.csv"))
        step = len(dr) // N_KEEP
        x = dr.iloc[::step][cols].iloc[:N_KEEP].reset_index(drop=True)
        name = f"tumor_{model}_{variant}.csv"
        path = os.path.join(DRAWS, name)
        x.to_csv(path, index=False, float_format="%.10g", lineterminator="\n")
        man["files"][name] = dict(
            sha256=hashlib.sha256(open(path, "rb").read()).hexdigest(), n_draws=len(x), thinning=step,
            min_ess_bulk=s["min_ess_bulk"], max_rhat=s["max_rhat"], divergences=s["divergences"],
            ppc_inside_90=s["ppc_inside_90"], variant=variant, model=model)
        print(name, len(x), man["files"][name]["sha256"][:12])
json.dump(man, open(os.path.join(DRAWS, "manifest.json"), "w"), indent=1)

# Point values for the deterministic "refit" presets: medians of the
# quantities the data identify (alpha, lambda/alpha, g/alpha and Tk or C/alpha),
# with lambda and g rebuilt from them, so the preset's daily BED losses equal
# the posterior medians. alpha/beta is not identified by these data (its
# posterior equals its prior), so the preset keeps 10 Gy.
pts = {}
for model in ("di", "dd"):
    dr = pd.read_csv(os.path.join(DRAWS, f"tumor_{model}_primary.csv"))
    a = float(dr["alpha"].median())
    pt = dict(alpha=a, ab=10.0, lambda_=float((dr["lambda"] / dr["alpha"]).median()) * a,
              g=float((dr["g"] / dr["alpha"]).median()) * a,
              loss_per_day_after_onset=float((dr["lambda"] / dr["alpha"]).median()),
              growth_loss_per_day=float((dr["g"] / dr["alpha"]).median()))
    pt["lambda"] = pt.pop("lambda_")
    if model == "di":
        pt["Tk"] = float(dr["Tk"].median())
    else:
        pt["c_over_alpha"] = float(dr["c_over_alpha"].median())
    pts[model] = pt
json.dump(dict(source="draws/tumor_*_primary.csv", rule="medians of identified quantities", values=pts),
          open(os.path.join(DRAWS, "refit_point_values.json"), "w"), indent=1)
print(json.dumps(pts, indent=1))
