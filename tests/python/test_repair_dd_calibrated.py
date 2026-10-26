"""fit/repair_dd_calibrated.py: re-fits only gate-failing calibrated-DD fits, records each repair,
keeps the replaced file, and refuses a folder whose identity no longer matches the code or data."""
import importlib, json, os, sys
import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "fit"))
DATA = os.path.join(ROOT, "data", "refit_trials_2026_10_25c_main.csv")


def setup(monkeypatch, tmp_path, fail=(2,), tamper=None):
    out = tmp_path / "fitout"
    for k, v in dict(BED_FIT_DATA=DATA, BED_FIT_OUT=str(out), DDCAL_NDI="3", DDCAL_DI_VARIANT="primary").items():
        monkeypatch.setenv(k, v)
    import fit_refit, fit_dd_calibrated, repair_dd_calibrated
    importlib.reload(fit_refit); importlib.reload(fit_dd_calibrated)
    Rp = importlib.reload(repair_dd_calibrated)
    F, C = Rp.F, Rp.C
    di_dir = out / "di_primary"; di_dir.mkdir(parents=True)
    di = pd.DataFrame(dict(draw=range(5), alpha=0.04, ab=10.0, g=0.01, **{"lambda": 0.03}, Tk=30.0))
    di.to_csv(di_dir / "draws_params.csv", index=False)
    json.dump(dict(identity_sha256="di-x"), open(di_dir / "identity.json", "w"))
    dd = out / "dd_calibrated_onset7_n3"; dd.mkdir()
    pick = np.linspace(0, 4, 3).round().astype(int)
    ident = dict(parent_di_identity="di-x", parent_di_draws_sha256=F.sha_file(di_dir / "draws_params.csv"),
                 data_file_sha256=F.sha_file(DATA),
                 stan_sha256=F.sha_file(os.path.join(ROOT, "fit", "stan", "dd_cal.stan")),
                 common_sha256=F.sha_file(os.path.join(ROOT, "fit", "stan", "common_functions.stanfunctions")),
                 runner_sha256=F.sha_file(C.__file__), fit_refit_sha256=F.sha_file(F.__file__),
                 picks_sha256=__import__("hashlib").sha256(pick.tobytes()).hexdigest())
    if tamper:
        ident[tamper] = "other"
    json.dump(ident, open(dd / "identity.json", "w"))
    for idx in pick:
        pd.DataFrame(dict(di_draw=[idx], attempt=3, gate_pass=idx not in fail, max_rhat=1.012, tag="old")
                     ).to_csv(dd / f"di_draw_{idx:05d}.csv", index=False)
    calls = []

    def fake_fit(m, r, base, onset_min, idx, out_):
        calls.append(int(idx))
        return (pd.DataFrame(dict(di_draw=[idx], attempt=4, gate_pass=True, max_rhat=1.004, tag="new")),
                dict(max_rhat=1.004, worst_rhat="sigma_raw"), 4, True)
    monkeypatch.setattr(Rp, "fit_draw", fake_fit)
    monkeypatch.setattr(F, "model", lambda name: object())
    monkeypatch.setattr(C, "summarize", lambda *a: None)
    return Rp, dd, calls


def test_only_failing_fits_are_refitted_and_recorded(monkeypatch, tmp_path):
    Rp, dd, calls = setup(monkeypatch, tmp_path, fail=(2,))
    Rp.main()
    assert calls == [2]
    assert pd.read_csv(dd / "di_draw_00002.csv").tag.iloc[0] == "new"
    assert pd.read_csv(dd / "di_draw_00000.csv").tag.iloc[0] == "old"
    assert pd.read_csv(dd / "repaired_originals" / "di_draw_00002.csv").tag.iloc[0] == "old"
    rep = json.load(open(dd / "repair.json"))
    assert [r["di_draw"] for r in rep["repairs"]] == [2]
    assert rep["repairs"][0]["rhat_before"] == pytest.approx(1.012) and rep["repairs"][0]["gate_pass"]
    assert rep["attempts"] == [4, 5]
    calls.clear(); Rp.main()                      # a second run finds nothing left to repair
    assert calls == []


@pytest.mark.parametrize("field", ["data_file_sha256", "runner_sha256", "stan_sha256", "parent_di_identity"])
def test_refuses_a_folder_whose_identity_does_not_match(monkeypatch, tmp_path, field):
    Rp, dd, calls = setup(monkeypatch, tmp_path, tamper=field)
    with pytest.raises(SystemExit) as e:
        Rp.main()
    assert "REFUSED" in str(e.value) and field in str(e.value) and calls == []
