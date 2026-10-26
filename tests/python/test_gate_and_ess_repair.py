"""Review 2026-10-26, R01 and R06: one fail-closed gate including ESS, ESS-aware repair, equal outer weights."""
import hashlib, importlib, json, os, sys
import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "fit"))
import gate as G

GOOD = dict(divergences=0, treedepth_hits=0, max_rhat=1.001, min_ess_bulk=900, min_ess_tail=800, min_ebfmi=0.9)


def test_gate_fails_low_ess_even_when_everything_else_is_clean():
    assert G.passes(GOOD)
    assert not G.passes(dict(GOOD, min_ess_bulk=150))
    assert not G.passes(dict(GOOD, min_ess_tail=399.9))
    assert G.failures(dict(GOOD, min_ess_bulk=150)) == ["min_ess_bulk 150"]


@pytest.mark.parametrize("key", list(GOOD))
def test_gate_fails_closed_on_missing_or_nonfinite(key):
    assert not G.passes({k: v for k, v in GOOD.items() if k != key})
    assert not G.passes(dict(GOOD, **{key: float("nan")}))


def test_fit_diagnostics_counts_divergences_from_rows():
    import fit_dd_calibrated as C
    rows = pd.DataFrame(dict(divergent=[0, 1, 0, 1], treedepth_hits=0, max_rhat=1.001, min_ess_bulk=900,
                             min_ess_tail=800, min_ebfmi=0.9))
    d = C.fit_diagnostics(rows)
    assert d["divergences"] == 2 and not G.passes(d)
    assert G.passes(C.fit_diagnostics(rows.assign(divergent=0)))


def test_equal_outer_weights_ignore_chain_length():
    import fit_dd_calibrated as C
    rng = np.random.default_rng(3)
    grp = np.repeat([0, 1, 2], [100, 100, 100]); x = rng.normal(grp, 1.0)
    w = C.outer_weights(grp)
    assert np.allclose([w[grp == g].sum() for g in range(3)], 1 / 3)
    # outer draw 2 rerun with 8x the rows (same distribution): summaries must not move
    grp8 = np.concatenate([grp, np.repeat(2, 700)]); x8 = np.concatenate([x, np.tile(x[grp == 2], 7)])
    w8 = C.outer_weights(grp8)
    for p in (.1, .5, .9):
        assert C.wquantile(x8, w8, p) == pytest.approx(C.wquantile(x, w, p), abs=1e-9)
    assert np.sum(w8 * x8) == pytest.approx(np.sum(w * x))
    # the unweighted quantile does move: this is the defect being fixed
    assert abs(np.quantile(x8, .5) - np.quantile(x, .5)) > 0.1


DATA = os.path.join(ROOT, "data", "refit_trials_2026_10_25c_main.csv")


def setup(monkeypatch, tmp_path, ess=(900, 150, 900), tamper=None):
    out = tmp_path / "fitout"
    for k, v in dict(BED_FIT_DATA=DATA, BED_FIT_OUT=str(out), DDCAL_NDI="3", DDCAL_DI_VARIANT="primary").items():
        monkeypatch.setenv(k, v)
    import fit_refit, fit_dd_calibrated, repair_dd_ess
    importlib.reload(fit_refit); importlib.reload(fit_dd_calibrated)
    Rp = importlib.reload(repair_dd_ess)
    F, C = Rp.F, Rp.C
    di_dir = out / "di_primary"; di_dir.mkdir(parents=True)
    pd.DataFrame(dict(draw=range(5), alpha=0.04, ab=10.0, g=0.01, **{"lambda": 0.03}, Tk=30.0)).to_csv(di_dir / "draws_params.csv", index=False)
    json.dump(dict(identity_sha256="di-x", data_file_sha256=F.sha_file(DATA)), open(di_dir / "identity.json", "w"))
    dd = out / "dd_calibrated_onset7_n3"; dd.mkdir()
    pick = np.linspace(0, 4, 3).round().astype(int)
    ident = dict(identity_sha256="dd-x", parent_di_identity="di-x", parent_di_draws_sha256=F.sha_file(di_dir / "draws_params.csv"),
                 data_file_sha256=F.sha_file(DATA), stan_sha256=F.sha_file(os.path.join(ROOT, "fit", "stan", "dd_cal.stan")),
                 common_sha256=F.sha_file(os.path.join(ROOT, "fit", "stan", "common_functions.stanfunctions")),
                 runner_sha256="an older runner", picks_sha256=hashlib.sha256(pick.tobytes()).hexdigest())
    if tamper:
        ident[tamper] = "other"
    json.dump(ident, open(dd / "identity.json", "w"))
    for idx, e in zip(pick, ess):
        pd.DataFrame(dict(di_draw=[idx] * 2, divergent=0, treedepth_hits=0, max_rhat=1.005, min_ess_bulk=e,
                          min_ess_tail=500, min_ebfmi=0.8, gate_pass=True, tag="old")).to_csv(dd / f"di_draw_{idx:05d}.csv", index=False)
    calls = []

    def fake_fit(m, r, base, onset_min, idx, out_):
        calls.append(int(idx))
        return (pd.DataFrame(dict(di_draw=[idx], divergent=0, treedepth_hits=0, max_rhat=1.002, min_ess_bulk=700,
                                  min_ess_tail=600, min_ebfmi=0.8, gate_pass=True, tag="new")),
                dict(min_ess_bulk=700, min_ess_tail=600, max_rhat=1.002), 4, True)
    monkeypatch.setattr(Rp, "fit_draw", fake_fit)
    monkeypatch.setattr(F, "model", lambda name: object())
    monkeypatch.setattr(C, "summarize", lambda *a: None)
    return Rp, dd, calls


def test_ess_repair_refits_only_fits_failing_the_full_gate(monkeypatch, tmp_path):
    # the middle fit has gate_pass True stored but bulk ESS 150: it must be refitted
    Rp, dd, calls = setup(monkeypatch, tmp_path)
    Rp.main()
    assert calls == [2]
    assert pd.read_csv(dd / "di_draw_00002.csv").tag.iloc[0] == "new"
    assert pd.read_csv(dd / "repaired_originals_ess" / "di_draw_00002.csv").tag.iloc[0] == "old"
    rep = json.load(open(dd / "repair_ess.json"))
    assert rep["repairs"][0]["reasons_before"] == ["min_ess_bulk 150"] and rep["gate_policy"]["min_ess_bulk"] == 400
    assert rep["original_runner_sha256"] == "an older runner"
    calls.clear(); Rp.main()
    assert calls == []


@pytest.mark.parametrize("field", ["data_file_sha256", "parent_di_draws_sha256", "stan_sha256", "parent_di_identity", "picks_sha256"])
def test_ess_repair_refuses_changed_inputs(monkeypatch, tmp_path, field):
    Rp, dd, calls = setup(monkeypatch, tmp_path, tamper=field)
    with pytest.raises(SystemExit) as e:
        Rp.main()
    assert "REFUSED" in str(e.value) and field in str(e.value) and calls == []
