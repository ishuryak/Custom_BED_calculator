"""The one diagnostic gate (plan decision 9, review R01 of 2026-10-26).

Used by the DI fit, the calibrated-DD fits, their repair, the export and its tests, so no stage
can apply a narrower gate than the one declared. Every check is required: a missing or non-finite
diagnostic fails (fail closed), never passes by omission.
"""
import math

POLICY = dict(version="gate-2026-10-27", divergences=0, treedepth_hits=0, max_rhat=1.01,
              min_ess_bulk=400, min_ess_tail=400, min_ebfmi=0.3,
              quantities="every sampled parameter (vectors included) and every regimen contrast; each one must have "
                         "a finite R-hat and ESS of its own, checked before any minimum or maximum is taken, and "
                         "none may be constant across draws (no model in this project samples a fixed quantity)")


def failures(dg):
    """List of reasons a fit fails the gate (empty when it passes). dg: dict of diagnostics."""
    out = []
    def num(k):
        v = dg.get(k)
        try:
            v = float(v)
        except (TypeError, ValueError):
            return None
        return v if math.isfinite(v) else None
    # Counts must be the integer 0 exactly (a negative count is invalid, not "better", review F02).
    count = lambda v: v == int(v) and v == 0
    checks = [("divergences", count),
              ("treedepth_hits", count),
              ("max_rhat", lambda v: 0 < v <= POLICY["max_rhat"]),
              ("min_ess_bulk", lambda v: v >= POLICY["min_ess_bulk"]),
              ("min_ess_tail", lambda v: v >= POLICY["min_ess_tail"]),
              ("min_ebfmi", lambda v: v >= POLICY["min_ebfmi"])]
    for k, ok in checks:
        v = num(k)
        if v is None:
            out.append(f"{k} missing or not finite")
        elif not ok(v):
            out.append(f"{k} {v:g}")
    # Per-quantity problems found before aggregation (fit_refit.diagnostics). Absent in diagnostics
    # rebuilt from stored per-fit rows, which carry only the aggregates.
    bad = dg.get("invalid_quantities")
    if bad:
        out.append(f"{len(bad)} quantities with missing, non-finite or constant diagnostics ({', '.join(bad[:3])})")
    return out


def invalid_quantities(table, required):
    """Names failing the per-quantity check, before any min/max (review F02 of 2026-10-26).

    table: DataFrame indexed by quantity name with columns rhat, ess_bulk, ess_tail, sd.
    required: base names (before "[") that must be present. A missing, non-finite or non-scalar
    value, or a quantity constant across draws, is a failure: pandas min/max would skip NaN.
    """
    import numpy as np
    present = {i.split("[")[0] for i in table.index}
    out = [f"{r} absent" for r in required if r not in present]
    for name, row in table.iterrows():
        vals = []
        for c in ("rhat", "ess_bulk", "ess_tail", "sd"):
            try:
                vals.append(float(row[c]))
            except (TypeError, ValueError):
                vals.append(float("nan"))
        if not all(np.isfinite(vals)):
            out.append(f"{name} not finite")
        elif vals[3] <= 0:
            out.append(f"{name} constant")
    return out


def passes(dg):
    return not failures(dg)
