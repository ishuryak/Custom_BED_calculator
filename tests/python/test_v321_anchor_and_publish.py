"""Version 3.2.1, review F03 (anchor residual in Gy of BED) on the shipped draws."""
import os, sys
import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "fit"))
import export_draws as E
import fit_dd_calibrated as C

DRAWS = os.path.join(ROOT, "draws")


def pair(variant):
    return (pd.read_csv(os.path.join(DRAWS, f"tumor_di_{variant}.csv")),
            pd.read_csv(os.path.join(DRAWS, f"tumor_dd_{variant}.csv")))


@pytest.mark.parametrize("variant, expected", [("primary", 9.83997e-7), ("paper_restrictions", 7.75184e-7)])
def test_shipped_pairs_pass_in_gy_and_match_an_independent_bed_evaluation(variant, expected):
    di, dd = pair(variant)
    res = E.anchor_residual(di, dd)
    assert res.max() <= E.ANCHOR_TOL
    assert res.max() == pytest.approx(expected, rel=1e-4)
    # independent: evaluate both anchor BEDs and subtract
    n0, d0, T0 = C.ANCHOR
    rd = C.rep_di(di["lambda"].to_numpy(), di.Tk.to_numpy(), T0)
    rr = C.rep_dd(dd.alpha.to_numpy(), dd.ab.to_numpy(), dd["lambda"].to_numpy(), dd.c_over_alpha.to_numpy(), n0, d0, T0)
    b_di = C.eff_bed(di.alpha.to_numpy(), di.ab.to_numpy(), di.g.to_numpy(), rd, n0, d0, T0)
    b_dd = C.eff_bed(dd.alpha.to_numpy(), dd.ab.to_numpy(), dd.g.to_numpy(), rr, n0, d0, T0)
    assert np.allclose(np.abs(b_dd - b_di), res, rtol=1e-6, atol=1e-12)


def test_review_boundary_pair_fails():
    """Log-cell residual 5e-5 at alpha 0.044858313 is 0.0011146 Gy, eleven times the tolerance."""
    di, dd = pair("primary")
    di, dd = di.iloc[[0]].copy(), dd.iloc[[0]].copy()
    assert dd.alpha.iloc[0] == pytest.approx(0.044858313)
    n0, d0, T0 = C.ANCHOR
    r_di = C.rep_di(di["lambda"].to_numpy(), di.Tk.to_numpy(), T0)[0]
    r_dd = C.rep_dd(dd.alpha.to_numpy(), dd.ab.to_numpy(), dd["lambda"].to_numpy(), dd.c_over_alpha.to_numpy(), n0, d0, T0)[0]
    dd["lambda"] = dd["lambda"] * (r_di + 5e-5) / r_dd          # rep_dd is linear in lambda_dd
    res = E.anchor_residual(di, dd)[0]
    assert res == pytest.approx(5e-5 / 0.044858313, rel=1e-6)
    assert res == pytest.approx(0.0011146206, rel=1e-6)
    assert res > E.ANCHOR_TOL


def test_nonpositive_or_mismatched_alpha_is_refused():
    di, dd = pair("primary")
    di, dd = di.iloc[:3].copy(), dd.iloc[:3].copy()
    for bad in (0.0, -0.1, np.nan):
        d2 = dd.copy(); d2.loc[d2.index[0], "alpha"] = bad
        with pytest.raises(SystemExit):
            E.anchor_residual(di, d2)


# --- Version 3.2.1, review F05: publication restores the previous release on failure ---------------

def release(path, tag):
    path.mkdir()
    (path / "draws_manifest.json").write_text(tag)
    return path


def failing_rename(when):
    """os.rename that raises on the call whose (src, dst) satisfies when."""
    def r(src, dst):
        if when(str(src), str(dst)):
            raise OSError("injected failure")
        os.rename(src, dst)
    return r


def test_failure_on_the_second_rename_restores_the_previous_release(tmp_path):
    dest = release(tmp_path / "draws", "old"); stage = release(tmp_path / "stage", "new")
    with pytest.raises(OSError):
        E.publish(str(stage), str(dest), rename=failing_rename(lambda s, d: s == str(stage)))
    assert (dest / "draws_manifest.json").read_text() == "old"
    assert not (tmp_path / "draws.previous").exists()


def test_failure_on_the_first_rename_leaves_the_release_in_place(tmp_path):
    dest = release(tmp_path / "draws", "old"); stage = release(tmp_path / "stage", "new")
    with pytest.raises(OSError):
        E.publish(str(stage), str(dest), rename=failing_rename(lambda s, d: s == str(dest)))
    assert (dest / "draws_manifest.json").read_text() == "old"


def test_failure_during_rollback_keeps_the_backup_and_the_next_run_recovers_it(tmp_path):
    dest = release(tmp_path / "draws", "old"); stage = release(tmp_path / "stage", "new")
    with pytest.raises(OSError):
        E.publish(str(stage), str(dest), rename=failing_rename(lambda s, d: s in (str(stage), str(dest) + ".previous")))
    assert not dest.exists() and (tmp_path / "draws.previous" / "draws_manifest.json").read_text() == "old"
    # a later attempt that fails again still leaves a complete prior release
    with pytest.raises(OSError):
        E.publish(str(stage), str(dest), rename=failing_rename(lambda s, d: s == str(stage)))
    assert (dest / "draws_manifest.json").read_text() == "old"
    # and a clean attempt publishes the new one
    E.publish(str(stage), str(dest))
    assert (dest / "draws_manifest.json").read_text() == "new" and not (tmp_path / "draws.previous").exists()


def test_successful_publication_exposes_only_the_new_release(tmp_path):
    dest = release(tmp_path / "draws", "old"); stage = release(tmp_path / "stage", "new")
    E.publish(str(stage), str(dest))
    assert (dest / "draws_manifest.json").read_text() == "new"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["draws"]


def test_an_incomplete_stage_is_not_accepted_and_the_previous_release_survives(tmp_path):
    dest = release(tmp_path / "draws", "old"); stage = tmp_path / "stage"; stage.mkdir()
    with pytest.raises(RuntimeError):
        E.publish(str(stage), str(dest))
    assert (tmp_path / "draws.previous" / "draws_manifest.json").read_text() == "old"
