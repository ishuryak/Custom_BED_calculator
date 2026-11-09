"""Jackson 1997 accelerated arm: how much does the dose encoding matter at the fitted DI parameters?

Review of 2026-11-09, item 2.1 step 3 (first half). No Stan fit is run. The fitting data give all 41
twice-daily patients 66 Gy in 33 x 2 Gy, but the trial paper (Jackson et al. 1997, p.41) states that
the last ten received 62.8 Gy because the second daily fraction was reduced to 1.8 Gy (17 x 2 + 16 x 1.8).
Three encodings of that arm are compared, with the once-daily arm unchanged:
  (a) uniform: 66 Gy, sum d^2 = 132 Gy^2 (what the fit used);
  (b) weighted exposure: prescription-weighted means 31/41 and 10/41, 65.219512 Gy and 129.034146 Gy^2
      (an approximation: a mean exposure inside a nonlinear TCP model);
  (c) prescription mixture: 31/41 x TCP(66 Gy) + 10/41 x TCP(62.8 Gy), on the probability scale.

The fit's DI link (fit/stan/di.stan, fit/fit_refit.py lines 98-112):
  log(-log TCP_i) = c_j - alpha (B_i - Bbar_j) + lambda s(T_i - Tk) + g (T_i - Tbar_j) + sigma z_i,
  B = sum d + sum d^2 / (alpha/beta), s(x) = log(1 + exp(3x)) / 3.
The unit intercepts c_j are not exported with the draws, so the shared Jackson baseline is set by holding
the once-daily arm at the fitted count rate 18/41 (round(41 x 0.443) / 41) with its arm deviation at zero.
The accelerated arm then follows from the within-unit difference, in which the centering cancels:
  eta_acc = eta_od - alpha (B_acc - B_od) + lambda [s(T_acc - Tk) - s(T_od - Tk)] + g (T_acc - T_od).
Both arms keep their delivered mean times, 24.25 and 47.6 days (paper, Table 3, p.41), and the single
pooled N = 41 and outcome per arm.

Run from the repository root:  python3 fit/jackson_exposure_check.py
"""
import csv, json, math

AB_ROWS = "draws/tumor_di_primary.csv"
POINT = "draws/refit_point_values.json"
D_UNI, D2_UNI = 66.0, 132.0                       # 33 x 2 Gy
D_AMD, D2_AMD = 17 * 2 + 16 * 1.8, 17 * 4 + 16 * 1.8 ** 2   # 62.8 Gy, 119.84 Gy^2
W_AMD = 10 / 41
D_W, D2_W = (1 - W_AMD) * D_UNI + W_AMD * D_AMD, (1 - W_AMD) * D2_UNI + W_AMD * D2_AMD
T_OD, T_ACC, N_OD, P_OD_REPORTED = 47.6, 24.25, 41, 0.443
P_OD = round(N_OD * P_OD_REPORTED) / N_OD


def s(x):
    z = 3 * x
    return x if z > 30 else math.log1p(math.exp(z)) / 3


def arm(alpha, ab, lam, Tk, g, dsum, d2sum, p_od=P_OD):
    """Predicted accelerated-arm TCP and the acc-minus-od contrast on the log(-log) scale."""
    b_od, b_acc = D_UNI + D2_UNI / ab, dsum + d2sum / ab
    d_eta = -alpha * (b_acc - b_od) + lam * (s(T_ACC - Tk) - s(T_OD - Tk)) + g * (T_ACC - T_OD)
    eta_od = math.log(-math.log(p_od))
    return math.exp(-math.exp(eta_od + d_eta)), d_eta


def encodings(alpha, ab, lam, Tk, g):
    pa, ea = arm(alpha, ab, lam, Tk, g, D_UNI, D2_UNI)
    pb, eb = arm(alpha, ab, lam, Tk, g, D_W, D2_W)
    p62, _ = arm(alpha, ab, lam, Tk, g, D_AMD, D2_AMD)
    pc = (1 - W_AMD) * pa + W_AMD * p62
    return dict(uniform=pa, weighted=pb, mixture=pc, eta_uniform=ea, eta_weighted=eb,
                bed_uniform=D_UNI + D2_UNI / ab, bed_weighted=D_W + D2_W / ab)


def q(v, p):
    v = sorted(v); k = (len(v) - 1) * p; f = math.floor(k)
    return v[f] + (v[min(f + 1, len(v) - 1)] - v[f]) * (k - f)


def main():
    print(f"amended prescription: {D_AMD:.1f} Gy, sum d^2 {D2_AMD:.2f} Gy^2; weighted {D_W:.6f} Gy, {D2_W:.6f} Gy^2")
    print(f"BED at alpha/beta 10: uniform {D_UNI + D2_UNI / 10:.6f}, weighted {D_W + D2_W / 10:.6f} Gy")
    print(f"once-daily baseline TCP {P_OD:.4f} (count {round(N_OD * P_OD_REPORTED)}/{N_OD})")
    pv = json.load(open(POINT))["values"]["di"]
    e = encodings(pv["alpha"], pv["ab"], pv["lambda"], pv["Tk"], pv["g"])
    print("\nPoint values (draws/refit_point_values.json):")
    for k in ("uniform", "weighted", "mixture"):
        print(f"  {k:9s} accelerated TCP {100 * e[k]:.2f}%   acc - od {100 * (e[k] - P_OD):+.2f} points")
    print(f"  contrast on log(-log) scale: uniform {e['eta_uniform']:+.4f}, weighted {e['eta_weighted']:+.4f}")
    rows = list(csv.DictReader(open(AB_ROWS)))
    seen, res = set(), []
    for r in rows:
        if r["di_draw"] in seen:
            continue
        seen.add(r["di_draw"])
        res.append(encodings(*(float(r[k]) for k in ("alpha", "ab", "lambda", "Tk", "g"))))
    print(f"\nAcross {len(res)} distinct DI draws (draws/tumor_di_primary.csv), median (2.5%, 97.5%):")
    for k in ("uniform", "weighted", "mixture"):
        v = [100 * x[k] for x in res]
        print(f"  {k:9s} accelerated TCP {q(v, .5):.2f}% ({q(v, .025):.2f}, {q(v, .975):.2f})")
    for k in ("weighted", "mixture"):
        v = [100 * (x[k] - x["uniform"]) for x in res]
        print(f"  {k:9s} minus uniform {q(v, .5):+.2f} points ({q(v, .025):+.2f}, {q(v, .975):+.2f})")
    v = [x["bed_weighted"] - x["bed_uniform"] for x in res]
    print(f"  LQ BED of the arm, weighted minus uniform: {q(v, .5):+.3f} Gy ({q(v, .025):+.3f}, {q(v, .975):+.3f})")


if __name__ == "__main__":
    main()
