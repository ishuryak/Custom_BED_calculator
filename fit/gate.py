"""The one diagnostic gate (plan decision 9, review R01 of 2026-10-26).

Used by the DI fit, the calibrated-DD fits, their repair, the export and its tests, so no stage
can apply a narrower gate than the one declared. Every check is required: a missing or non-finite
diagnostic fails (fail closed), never passes by omission.
"""
import math

POLICY = dict(version="gate-2026-10-26", divergences=0, treedepth_hits=0, max_rhat=1.01,
              min_ess_bulk=400, min_ess_tail=400, min_ebfmi=0.3,
              quantities="every sampled parameter (vectors included) and every regimen contrast and onset; "
                         "quantities constant across draws (fixed parameters) have no R-hat or ESS and are excluded")


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
    checks = [("divergences", lambda v: v <= POLICY["divergences"]),
              ("treedepth_hits", lambda v: v <= POLICY["treedepth_hits"]),
              ("max_rhat", lambda v: v <= POLICY["max_rhat"]),
              ("min_ess_bulk", lambda v: v >= POLICY["min_ess_bulk"]),
              ("min_ess_tail", lambda v: v >= POLICY["min_ess_tail"]),
              ("min_ebfmi", lambda v: v >= POLICY["min_ebfmi"])]
    for k, ok in checks:
        v = num(k)
        if v is None:
            out.append(f"{k} missing or not finite")
        elif not ok(v):
            out.append(f"{k} {v:g}")
    return out


def passes(dg):
    return not failures(dg)
