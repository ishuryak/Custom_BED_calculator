"""The one diagnostic gate (plan decision 9, review R01 of 2026-10-26).

Used by the DI fit, the calibrated-DD fits, their repair, the export and its tests, so no stage
can apply a narrower gate than the one declared. Every check is required: a missing or non-finite
diagnostic fails (fail closed), never passes by omission.

Constituents are validated before any reduction (review F06 of 2026-10-31): each chain's energy
series and E-BFMI, the sampler arrays, and the full expected list of quantity names. Python's
min() over [finite, nan] can return the finite value, so a bad chain must be found per chain.
"""
import math
import numbers

POLICY = dict(version="gate-2026-10-27", divergences=0, treedepth_hits=0, max_rhat=1.01,
              min_ess_bulk=400, min_ess_tail=400, min_ebfmi=0.3,
              quantities="every sampled parameter (vectors included) and every regimen contrast; each one must have "
                         "a finite R-hat and ESS of its own, checked before any minimum or maximum is taken, and "
                         "none may be constant across draws (no model in this project samples a fixed quantity)")


def scalar(v):
    """One finite number, or None. Booleans are rejected, as R's draws_gate_failures rejects
    logicals (review F06 of 2026-10-31: float(False) is 0, which passed as a zero count). Numeric
    strings are accepted, as R's as.numeric accepts them (the shared fixture table is text)."""
    import numpy as np
    if isinstance(v, (bool, np.bool_, bytes, bytearray)) or v is None:
        return None
    if not isinstance(v, (numbers.Real, str)):
        return None
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def failures(dg):
    """List of reasons a fit fails the gate (empty when it passes). dg: dict of diagnostics."""
    out = []
    # Counts must be the integer 0 exactly (a negative count is invalid, not "better", review F02).
    count = lambda v: v == int(v) and v == 0
    checks = [("divergences", count),
              ("treedepth_hits", count),
              ("max_rhat", lambda v: 0 < v <= POLICY["max_rhat"]),
              ("min_ess_bulk", lambda v: v >= POLICY["min_ess_bulk"]),
              ("min_ess_tail", lambda v: v >= POLICY["min_ess_tail"]),
              ("min_ebfmi", lambda v: v >= POLICY["min_ebfmi"])]
    for k, ok in checks:
        v = scalar(dg.get(k))
        if v is None:
            out.append(f"{k} missing or not finite")
        elif not ok(v):
            out.append(f"{k} {v:g}")
    # Per-quantity and per-chain problems found before aggregation (fit_refit.diagnostics). Absent
    # in diagnostics rebuilt from stored per-fit rows, which carry only the aggregates.
    bad = dg.get("invalid_quantities")
    if bad:
        out.append(f"{len(bad)} quantities with missing, non-finite or constant diagnostics ({', '.join(bad[:3])})")
    bad = dg.get("invalid_sampler")
    if bad:
        out.append(f"{len(bad)} sampler problems ({', '.join(bad[:3])})")
    return out


def ebfmi_per_chain(energy):
    """E-BFMI of each chain, validated before any minimum is taken.

    energy: array (draws, chains). Returns (values, problems): values[k] is a finite positive
    float, or None when chain k's energy is non-finite, constant or too short, or its E-BFMI is
    not finite. problems lists each bad chain (1-based); empty only when every chain is valid.
    """
    import numpy as np
    try:
        e = np.asarray(energy, dtype=float)
    except (TypeError, ValueError):
        return [], ["energy is not numeric"]
    if e.ndim == 3 and e.shape[2] == 1:
        e = e[:, :, 0]
    if e.ndim != 2 or e.shape[0] < 3 or e.shape[1] < 1:
        return [], [f"energy has shape {e.shape}, not (draws >= 3, chains >= 1)"]
    vals, problems = [], []
    for k in range(e.shape[1]):
        x = e[:, k]
        if not np.all(np.isfinite(x)):
            vals.append(None); problems.append(f"chain {k + 1} energy not finite"); continue
        den = float(np.sum((x - x.mean()) ** 2))
        if not den > 0:
            vals.append(None); problems.append(f"chain {k + 1} energy constant"); continue
        v = float(np.sum(np.diff(x) ** 2) / den)
        if not (math.isfinite(v) and v > 0):
            vals.append(None); problems.append(f"chain {k + 1} E-BFMI {v}"); continue
        vals.append(v)
    return vals, problems


def sampler_counts(divergent, treedepth, max_treedepth, shape=None):
    """Divergence and tree-depth-hit counts from per-draw sampler arrays, validated first.

    Both arrays must be numeric, finite, of the same (draws, chains) shape (equal to shape when
    given, normally the energy array's); divergent__ must be 0 or 1 and treedepth__ a
    non-negative integer. Returns (divergences, treedepth_hits, problems); the counts are None
    when any problem is found, so the gate fails closed.
    """
    import numpy as np
    problems = []
    arrs = {}
    for name, a in (("divergent__", divergent), ("treedepth__", treedepth)):
        try:
            a = np.asarray(a, dtype=float)
        except (TypeError, ValueError):
            problems.append(f"{name} not numeric"); continue
        if a.ndim == 3 and a.shape[2] == 1:
            a = a[:, :, 0]
        arrs[name] = a
    if problems:
        return None, None, problems
    d, t = arrs["divergent__"], arrs["treedepth__"]
    if d.shape != t.shape or d.size == 0:
        problems.append(f"divergent__ {d.shape} and treedepth__ {t.shape} disagree or are empty")
    if shape is not None and (d.shape != tuple(shape) or t.shape != tuple(shape)):
        problems.append(f"sampler arrays {d.shape} do not match draws x chains {tuple(shape)}")
    if not (np.all(np.isfinite(d)) and np.all(np.isin(d, (0.0, 1.0)))):
        problems.append("divergent__ not all 0 or 1")
    if not (np.all(np.isfinite(t)) and np.all(t >= 0) and np.all(t == np.round(t))):
        problems.append("treedepth__ not all non-negative integers")
    if problems:
        return None, None, problems
    return int(d.sum()), int((t >= max_treedepth).sum()), []


def expected_quantities(sampled, N, J, R, **_):
    """Every quantity name a fit must report: the scalar sampled parameters, c[1..J], z[1..N]
    and delta_eval[1..R], from the actual data dimensions (review F06 of 2026-10-31).

    Every sampled parameter in fit/stan (di, dd, dd_cal) is a scalar except c and z.
    """
    return (list(sampled) + [f"c[{j}]" for j in range(1, int(J) + 1)] + [f"z[{i}]" for i in range(1, int(N) + 1)]
            + [f"delta_eval[{r}]" for r in range(1, int(R) + 1)])


def invalid_quantities(table, required=(), expected=None):
    """Names failing the per-quantity check, before any min/max (review F02 of 2026-10-26).

    table: DataFrame indexed by quantity name with columns rhat, ess_bulk, ess_tail, sd.
    expected: the full list of names (expected_quantities). Every one must appear exactly once,
    and no other name of the same families may appear (a wrong index or shape). This is what the
    fitting call sites use (review F06 of 2026-10-31).
    required: base names (before "[") that must be present; kept for callers without dimensions.
    It cannot see a missing vector element or a missing c[]/z[] family.
    A missing, non-finite or non-scalar value, or a quantity constant across draws, is a failure:
    pandas min/max would skip NaN.
    """
    import numpy as np
    names = [str(i) for i in table.index]
    present = {i.split("[")[0] for i in names}
    out = [f"{r} absent" for r in required if r not in present]
    if expected is not None:
        exp = list(expected)
        have = set(names)
        out += [f"{e} absent" for e in exp if e not in have]
        seen = set()
        for n in names:
            if n in seen:
                out.append(f"{n} duplicated")
            seen.add(n)
        families = {e.split("[")[0] for e in exp}
        exp_set = set(exp)
        out += [f"{n} unexpected" for n in dict.fromkeys(names) if n not in exp_set and n.split("[")[0] in families]
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
