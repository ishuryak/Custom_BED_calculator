"""Posterior predictive check figure for the refit (offline). Writes
docs/development/2026_10_21_v3/fit_ppc.png (600 dpi). Run from the repo root."""
import json, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
plt.style.use(os.path.expanduser("~/.claude/figure-style/publication.mplstyle"))
OUT = "fit/output"
fig, axes = plt.subplots(1, 2, figsize=(13, 6.2), sharey=True)
for ax, model, title in zip(axes, ["di", "dd"], ["DI model (refit)", "DD model (refit)"]):
    s = json.load(open(f"{OUT}/{model}_primary/summary.json"))
    p = s["ppc"]
    obs = np.array([x["obs_pct"] for x in p]); med = np.array([x["pred_median_pct"] for x in p])
    lo = np.array([x["pred_lo_pct"] for x in p]); hi = np.array([x["pred_hi_pct"] for x in p])
    n = np.array([x["n"] for x in p])
    ax.errorbar(med, obs, xerr=[med - lo, hi - med], fmt="none", ecolor="#999999", elinewidth=1, alpha=0.8)
    ax.scatter(med, obs, s=12 + 0.12 * n, c="#0072B2" if model == "di" else "#D55E00", edgecolor="black", linewidth=0.5, zorder=3)
    ax.plot([0, 100], [0, 100], color="black", linewidth=1)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.set_xlabel("Predicted tumor control (%)\nmedian and 90% predictive interval")
    ax.set_title(f"{title}: {s['ppc_inside_90']} of 35 arms inside", fontsize=15, fontweight="bold")
axes[0].set_ylabel("Observed tumor control (%)")
fig.tight_layout(w_pad=3)
os.makedirs("docs/development/2026_10_21_v3", exist_ok=True)
fig.savefig("docs/development/2026_10_21_v3/fit_ppc.png", dpi=600, bbox_inches="tight")
print("written")
