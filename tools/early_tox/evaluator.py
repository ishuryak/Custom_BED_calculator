"""Early (acute mucosal) BED for an explicit fraction schedule (plan R8, step 4 groundwork).

One evaluator for fitting and prediction:

    E = D + Q / ab - r * max(0, T - Tk)

    D  total dose (Gy)
    Q  sum_i d_i^2 + 2 * sum_{i<j} d_i d_j exp(-mu (t_j - t_i)),  mu = ln 2 / h   (incomplete repair)
       or sum_i d_i^2 when h is None (complete repair between fractions)
    T  elapsed days from the first to the last fraction ((t_last - t_first) / 24)
    r  BED units per day = ln 2 / (alpha * Tp) in the Fowler parameterization

Times are hours from the first fraction, as in data/early_tox/*fractions.csv. The clock is
elapsed days, which is the app's convention for Fowler's acute parameters (R/bed_models.R,
acute_model_time "elapsed"). The hinge is exact here. The app's smooth max is a separate
legacy evaluator (step 7) and is not reproduced by this module.
"""
import math

FOWLER_2003 = dict(
    ab=10.0,
    alpha=0.35,
    Tk=7.0,
    Tp=2.5,
)


def fowler_rate(alpha=FOWLER_2003["alpha"], Tp=FOWLER_2003["Tp"]):
    """BED units lost to repopulation per day: ln 2 / (alpha * Tp)."""
    return math.log(2) / (alpha * Tp)


def repair_q(doses, times_h, halftime_h=None):
    """Quadratic dose term Q with incomplete repair of half-time halftime_h (None = complete)."""
    if len(doses) != len(times_h):
        raise ValueError("doses and times differ in length")
    if any(b < a for a, b in zip(times_h, times_h[1:])):
        raise ValueError("fraction times are not ordered")
    q = sum(d * d for d in doses)
    if halftime_h is None:
        return q
    if halftime_h <= 0:
        raise ValueError("repair half-time must be positive")
    mu = math.log(2) / halftime_h
    for i in range(len(doses)):
        for j in range(i + 1, len(doses)):
            q += 2 * doses[i] * doses[j] * math.exp(-mu * (times_h[j] - times_h[i]))
    return q


def elapsed_days(times_h):
    return (times_h[-1] - times_h[0]) / 24.0


def early_bed(doses, times_h, ab=FOWLER_2003["ab"], r=None, Tk=FOWLER_2003["Tk"], halftime_h=None):
    """Early BED of one schedule. r defaults to the Fowler 2003 rate."""
    if not doses:
        raise ValueError("empty schedule")
    if ab <= 0:
        raise ValueError("alpha/beta must be positive")
    r = fowler_rate() if r is None else r
    D = sum(doses)
    return D + repair_q(doses, times_h, halftime_h) / ab - r * max(0.0, elapsed_days(times_h) - Tk)


def regular_schedule(n, d, days_per_week=5, per_day=1, interval_h=6.0):
    """n fractions of d Gy from a Monday, for fixtures. Returns (doses, times_h)."""
    allowed = set(range(days_per_week))
    doses, times, day = [], [], 0
    while len(doses) < n:
        if day % 7 in allowed:
            for k in range(per_day):
                if len(doses) < n:
                    doses.append(d)
                    times.append(day * 24 + k * interval_h)
        day += 1
    return doses, times
