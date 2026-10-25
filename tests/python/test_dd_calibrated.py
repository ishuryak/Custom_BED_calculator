"""Calibrated DD (fit/fit_dd_calibrated.py): DD must equal DI at the anchor, and the
onset must match the app's own calibration in R/bed_models.R."""
import os, subprocess, sys
import numpy as np
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "fit"))
import fit_dd_calibrated as C

CASES = [  # alpha, ab, lambda_di, Tk, lambda_dd
    (0.039, 10.3, 0.0285, 38.6, 0.43),
    (0.044, 9.26, 0.0245, 34.2, 0.05),
    (0.025, 15.0, 0.06, 25.0, 3.0),
]


@pytest.mark.parametrize("a,ab,ldi,tk,ldd", CASES)
def test_dd_equals_di_at_anchor(a, ab, ldi, tk, ldd):
    n0, d0, T0 = C.ANCHOR
    _, co = C.calibrate(a, ab, ldi, tk, ldd)
    assert C.rep_dd(a, ab, ldd, co, n0, d0, T0) == pytest.approx(C.rep_di(ldi, tk, T0), rel=1e-9)


@pytest.mark.parametrize("a,ab,ldi,tk,ldd", CASES)
def test_dd_differs_from_di_away_from_anchor(a, ab, ldi, tk, ldd):
    # HYPNO has a higher BED per day than the anchor, so DD must not reduce to DI there.
    _, co = C.calibrate(a, ab, ldi, tk, ldd)
    assert abs(C.rep_dd(a, ab, ldd, co, 20, 2.75, 26) - C.rep_di(ldi, tk, 26)) > 1e-4


def test_onset_matches_app_R_calibration():
    import shutil
    win = "/mnt/c/Program Files/R/R-4.5.1/bin/Rscript.exe"
    rs = shutil.which("Rscript") or (win if os.path.exists(win) else None)
    if rs is None:
        pytest.skip("Rscript not available")
    for a, ab, ldi, tk, ldd in CASES:
        on_py, co_py = C.calibrate(a, ab, ldi, tk, ldd)
        expr = (f"source('R/bed_models.R'); r <- dd_calibration(list(m=35,d=2,T=47), {a}, {ab}, {ldi}, {ldd}, {tk});"
                f"cat(sprintf('%.12f %.12f %s', r$onset_ref, r$c_over_alpha, r$status))")
        out = subprocess.run([rs, "-e", expr], cwd=REPO, capture_output=True, text=True, check=True).stdout.split()
        assert float(out[0]) == pytest.approx(on_py, rel=1e-9, abs=1e-9)
        if out[2] == "ok":
            assert float(out[1]) == pytest.approx(co_py, rel=1e-9)


def test_inv_smax_inverts_smax():
    x = np.array([-5.0, -0.5, 0.0, 0.3, 4.0, 20.0])
    assert np.allclose(C.inv_smax(C.smax(x)), x, atol=1e-9)
