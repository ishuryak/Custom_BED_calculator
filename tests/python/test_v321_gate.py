"""Version 3.2.1, review F02 (2026-10-26 v3.2 review): the diagnostic gate fails closed on every
malformed value, per quantity before aggregation, and agrees with the R gate on one fixture table."""
import csv, os, sys
import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "fit"))
import gate as G

CASES = list(csv.DictReader(open(os.path.join(ROOT, "tests", "fixtures", "gate_cases.csv"))))
KEYS = ["divergences", "treedepth_hits", "max_rhat", "min_ess_bulk", "min_ess_tail", "min_ebfmi"]


@pytest.mark.parametrize("case", CASES, ids=[c["case"] for c in CASES])
def test_shared_fixture_table(case):
    dg = {k: (case[k] if case[k] != "" else None) for k in KEYS}
    dg = {k: v for k, v in dg.items() if v is not None}
    assert G.passes(dg) == (case["expect_pass"] == "TRUE")


@pytest.mark.parametrize("bad", [[1, 2], None, True])
def test_vector_null_and_bool_counts_fail(bad):
    good = {k: float(next(c for c in CASES if c["case"] == "good")[k]) for k in KEYS}
    assert G.passes(good)
    assert not G.passes(dict(good, divergences=bad))
    assert not G.passes(dict(good, min_ess_bulk=bad))


class FakeFit:
    """Minimal CmdStanMCMC stand-in for fit_refit.diagnostics."""
    def __init__(self, table):
        self.table = table
        self.column_names = ["lp__", "energy__"]
    def summary(self):
        return self.table
    def method_variables(self):
        return dict(divergent__=np.zeros(8), treedepth__=np.zeros(8))
    def draws(self, concat_chains=False):
        rng = np.random.default_rng(1)
        return rng.normal(size=(200, 4, 2))


def summary(rows):
    return pd.DataFrame(rows, columns=["name", "StdDev", "R_hat", "ESS_bulk", "ESS_tail"]).set_index("name")


GOOD_ROWS = [("lambda_dd", 0.1, 1.001, 900, 800), ("sigma_raw", 0.2, 1.002, 700, 700),
             ("c[1]", 1.0, 1.0, 1200, 900), ("z[1]", 1.0, 1.0, 1500, 1000), ("delta_eval[1]", 0.5, 1.0, 1000, 900)]


def test_diagnostics_passes_a_clean_summary():
    import fit_refit as F
    dg = F.diagnostics(FakeFit(summary(GOOD_ROWS)), ["lambda_dd", "sigma_raw"])
    assert G.passes(dg), G.failures(dg)


@pytest.mark.parametrize("i, col, val", [(0, "ESS_bulk", np.nan), (0, "R_hat", np.nan), (1, "ESS_tail", np.inf),
                                         (0, "StdDev", 0.0)])
def test_one_bad_quantity_among_good_ones_fails_before_aggregation(i, col, val):
    """Review counterexample: lambda_dd bulk ESS NaN beside sigma_raw ESS 700 gave min 700 and a pass."""
    import fit_refit as F
    rows = [list(r) for r in GOOD_ROWS]
    rows[i][["name", "StdDev", "R_hat", "ESS_bulk", "ESS_tail"].index(col)] = val
    dg = F.diagnostics(FakeFit(summary(rows)), ["lambda_dd", "sigma_raw"])
    assert not G.passes(dg), dg


def test_a_required_quantity_absent_from_the_summary_fails():
    import fit_refit as F
    dg = F.diagnostics(FakeFit(summary(GOOD_ROWS[1:])), ["lambda_dd", "sigma_raw"])
    assert not G.passes(dg)
