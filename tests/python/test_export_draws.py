"""fit/export_draws.py: refuses fits that fail the gate, and pairs DD row i with its DI parent row i."""
import importlib, json, os, sys
import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "fit"))


def fake_fit(src, n_di=6, n_outer=2, inner=10, di_pass=True, dd_fail=0):
    for v in ("primary", "paper_restrictions"):
        d = src / f"di_{v}"; d.mkdir(parents=True)
        rng = np.random.default_rng(1)
        di = pd.DataFrame(dict(draw=range(n_di), alpha=0.04 + 0.001 * np.arange(n_di), ab=10 + np.arange(n_di),
                               **{"lambda": 0.03 + 0.001 * np.arange(n_di)}, Tk=37.0, g=0.011))
        di.to_csv(d / "draws_params.csv", index=False)
        json.dump(dict(gate_pass=di_pass, divergences=0, treedepth_hits=0, max_rhat=1.001, worst_rhat="alpha",
                       min_ess_bulk=900, min_ess_tail=800, min_ebfmi=0.9, n_quantities=80, ppc_inside_90=45),
                  open(d / "summary.json", "w"))
        json.dump(dict(identity_sha256=f"di-{v}", data_file="refit_trials_2026_10_25c_main.csv"), open(d / "identity.json", "w"))
        dd = src / (f"dd_calibrated_onset7_n{n_outer}" + ("" if v == "primary" else f"_{v}")); dd.mkdir()
        for idx in np.linspace(0, n_di - 1, n_outer).round().astype(int):
            r = di.iloc[idx]
            pd.DataFrame(dict(di_draw=idx, alpha=r.alpha, ab=r.ab, g=r.g, **{"lambda": r["lambda"]}, Tk=r.Tk,
                              chain=np.repeat([1, 2], inner // 2), iter=np.tile(range(1, inner // 2 + 1), 2),
                              lambda_dd=0.4 + idx + 0.01 * np.arange(inner), c_over_alpha=300.0 + idx)
                         ).to_csv(dd / f"di_draw_{idx:05d}.csv", index=False)
        json.dump(dict(outer_fits_failing_gate=dd_fail, divergences=0, treedepth_hits=0, max_rhat=1.003,
                       min_ess_bulk_per_fit=500, min_ess_tail_per_fit=500, min_ebfmi=0.8, outer_fits_retried=0,
                       fitted_onset_below_min=dict(n_draws=0)), open(dd / "summary.json", "w"))
        json.dump(dict(identity_sha256=f"dd-{v}", parent_di_identity=f"di-{v}"), open(dd / "identity.json", "w"))


def run_export(monkeypatch, tmp_path, **kw):
    src, to = tmp_path / "fit", tmp_path / "draws"
    to.mkdir()
    (tmp_path / "data").mkdir(exist_ok=True)
    fake_fit(src, **kw)
    monkeypatch.setenv("BED_EXPORT_FROM", str(src)); monkeypatch.setenv("BED_EXPORT_DD_N", "2")
    monkeypatch.setenv("BED_EXPORT_TO", str(to))
    import export_draws
    E = importlib.reload(export_draws)
    monkeypatch.setattr(E, "ROOT", str(tmp_path))
    (tmp_path / "data" / "refit_trials_2026_10_25c_main.csv").write_text("x\n")
    E.main()
    return to


def test_dd_rows_are_paired_with_their_di_parent(monkeypatch, tmp_path):
    to = run_export(monkeypatch, tmp_path)
    di = pd.read_csv(to / "tumor_di_primary.csv"); dd = pd.read_csv(to / "tumor_dd_primary.csv")
    assert len(di) == len(dd) == 2 * 4
    assert (di.di_draw == dd.di_draw).all() and np.allclose(di.alpha, dd.alpha) and np.allclose(di.ab, dd.ab)
    # inner DD draws differ within an outer draw; DI repeats it
    assert dd.groupby("di_draw")["lambda"].nunique().min() == 4 and di.groupby("di_draw")["lambda"].nunique().max() == 1
    man = json.load(open(to / "draws_manifest.json"))
    assert man["files"]["tumor_dd_primary.csv"]["parent_di_identity"] == "di-primary"
    # the app refuses a file whose diagnostics do not say gate_pass TRUE, so DD must carry it too
    assert all(f["diagnostics"]["gate_pass"] is True for f in man["files"].values())


@pytest.mark.parametrize("kw", [dict(di_pass=False), dict(dd_fail=1)])
def test_export_refuses_a_failing_gate(monkeypatch, tmp_path, kw):
    with pytest.raises(SystemExit) as e:
        run_export(monkeypatch, tmp_path, **kw)
    assert "REFUSED" in str(e.value)
