"""Markdown tables of the comparison, generated from fit/output_frequentist/comparison/comparison.json
(itself assembled only from verified results by compare.py), so no number is typed by hand.

    python3 fit/frequentist/report_tables.py > fit/output_frequentist/comparison/tables.md
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from frequentist import fit_ml as FM  # noqa: E402

LABEL = dict(alpha="α", ab="α/β", lam="λ", Tk="Tk", g="g", sigma="σ (arm spread)", lam_plus_g="λ + g",
             rate_bed="(λ + g)/α", rate_eqd2="(λ + g)/α in EQD2", kappa="κ (slope)", r="r", gamma="γ (chemotherapy)",
             tau="τ (trial spread)", sigma_arm="σ_arm (arm spread)")
TUMOR_REG = ["HYPNO 55 Gy/20 (26 d)", "HYHOPE-type 46.5 Gy/15 (19 d)", "17 × 3.05 Gy (23 d)", "60 Gy/25 (33 d)",
             "65 Gy/30 (40 d)", "66 Gy/33 (45 d)", "30 × 2 Gy (40 d)", "33 × 2 Gy, 6 per week (38 d)",
             "70 Gy/35, 6 per week (40 d)", "68 × 1.2 Gy twice daily (47 d)"]
EARLY_REG = ["15 × 3.5 Gy", "20 × 2.75 Gy", "25 × 2.4 Gy", "30 × 2.2 Gy"]


def n(x, d=3):
    if x is None:
        return "–"
    a = abs(x)
    if a >= 100:
        return f"{x:.0f}".replace("-", "−")
    s = f"{x:.{d}g}" if (a >= 1e-3 or a == 0) else f"{x:.2e}"
    return s.replace("-", "−")


def label(t):
    if t.startswith("contrast["):
        return f"ΔBED, {TUMOR_REG[int(t[9:-1])]}"
    if t.startswith("lor["):
        return f"log OR, {EARLY_REG[int(t[4:-1])]}"
    return LABEL[t]


def prof(r):
    ps = r.get("profile_set") or []
    if not ps:
        return "not run"
    out = []
    for p in ps:
        lo = n(p["lower"]["value"]) + ("" if p["lower"]["closed"] else " (open)")
        hi = n(p["upper"]["value"]) + ("" if p["upper"]["closed"] else " (open)")
        out.append(f"{lo} to {hi}")
    return "; ".join(out) + (" (flat)" if r.get("profile_flat") else "")


def table(mm, with_boot):
    head = "| Target | Units | ML estimate | Profile 95% set | " + ("Bootstrap 2.5–97.5% | " if with_boot else "") \
        + "Posterior median (95%) | Prior 95% | Width ratio |"
    lines = [head, "|" + "---|" * (head.count("|") - 1)]
    for r in mm["rows"]:
        b = r.get("bootstrap")
        boot = f"{n(b['q025'])} to {n(b['q975'])}" if b else "–"
        post = r.get("posterior")
        post = f"{n(post[1])} ({n(post[0])} to {n(post[2])})" if post else "–"
        pri = r.get("prior")
        pri = f"{n(pri[0])} to {n(pri[2])}" if pri else "–"
        cells = [label(r["target"]), r["units"], n(r["mle"], 4), prof(r)] + ([boot] if with_boot else []) + \
            [post, pri, n(r.get("width_ratio"), 2)]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main():
    r = json.load(open(os.path.join(FM.OUT, "comparison", "comparison.json")))
    for key, title in (("tumor_matched", "Tumor DI, matched support"), ("tumor_expanded", "Tumor DI, expanded support"),
                       ("early_early_matched", "Acute mucositis (U02)")):
        mm = r["models"][key]
        print(f"### {title}\n")
        print("Log-likelihood at the estimate " + f"{mm['mle']['ll']:.6f}".replace("-", "−") + "; estimate on a limit: "
              + (", ".join(f"{b['parameter']} at its {b['side']} {'computational' if b['kind'] == 'comp' else 'support'} limit"
                           for b in mm["mle"]["boundary"]) or "none") + ".\n")
        print(table(mm, mm["bootstrap"] is not None) + "\n")
        if mm["bootstrap"]:
            b = mm["bootstrap"]
            print(f"Bootstrap: {b['B']} replicates, {b['ok']} converged, {b['unresolved_failures']} unresolved failures, "
                  f"{b['boundary_solutions']} with at least one parameter on a limit; {n(b['seconds_total'] / b['B'], 3)} "
                  "CPU-seconds per replicate.\n")
        if mm["sigma_sensitivity"]:
            print("Hybrid sensitivity (σ held at an externally chosen value; not a likelihood-only result):\n")
            print("| Regimen | σ = 0.065 (posterior median): estimate, profile 95% set | σ = 0.17 (posterior 97.5%): estimate, profile 95% set |")
            print("|---|---|---|")
            sens = mm["sigma_sensitivity"]
            for k in range(10):
                cells = []
                for s in ("0.065", "0.17"):
                    v = sens[f"contrast[{k}]__sigma_{s}"]
                    lo, hi, wl, wh = v["set"][0]
                    cells.append(f"{n(v['estimate'])}: {n(lo) if lo is not None else 'open'} to {n(hi) if hi is not None else 'open'}")
                print(f"| {TUMOR_REG[k]} | " + " | ".join(cells) + " |")
            print()


if __name__ == "__main__":
    main()
