"""fit/export_draws.py: pairs DD row i with its DI parent row i, refuses fits that fail the gate, and
refuses any break in the identity chain before replacing draws/ (review R02, 2026-10-26)."""
import hashlib, importlib, json, os, shutil, sys
import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "fit"))
import fit_dd_calibrated as C

DATA_NAME = "refit_trials_2026_10_25c_main.csv"


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def fake_fit(root, n_di=6, n_outer=2, inner=10, di_pass=True, low_ess_draw=None):
    src = root / "fit_out"
    stan = root / "fit" / "stan"; stan.mkdir(parents=True)
    for f in ("di.stan", "dd_cal.stan", "common_functions.stanfunctions"):
        shutil.copy(os.path.join(ROOT, "fit", "stan", f), stan / f)
    (root / "data").mkdir()
    data = root / "data" / DATA_NAME; data.write_text("x\n")
    picks = np.linspace(0, n_di - 1, n_outer).round().astype(int)
    for v in ("primary", "paper_restrictions"):
        d = src / f"di_{v}"; d.mkdir(parents=True)
        di = pd.DataFrame(dict(draw=range(n_di), alpha=0.04 + 0.001 * np.arange(n_di), ab=10 + np.arange(n_di),
                               **{"lambda": 0.03 + 0.001 * np.arange(n_di)}, Tk=37.0, g=0.011))
        di.to_csv(d / "draws_params.csv", index=False)
        json.dump(dict(gate_pass=di_pass, divergences=0, treedepth_hits=0, max_rhat=1.001, worst_rhat="alpha",
                       min_ess_bulk=900 if di_pass else 150, min_ess_tail=800, min_ebfmi=0.9, n_quantities=80,
                       ppc_inside_90=45), open(d / "summary.json", "w"))
        json.dump(dict(identity_sha256=f"di-{v}", variant=v, data_file=DATA_NAME, data_file_sha256=sha(data),
                       stan_sha256={f: sha(stan / f) for f in ("di.stan", "common_functions.stanfunctions")}),
                  open(d / "identity.json", "w"))
        dd = src / (f"dd_calibrated_onset7_n{n_outer}" + ("" if v == "primary" else f"_{v}")); dd.mkdir()
        for idx in picks:
            r = di.iloc[idx]
            lam = 0.4 + idx + 0.01 * np.arange(inner)
            co = np.array([C.calibrate(r.alpha, r.ab, r["lambda"], r.Tk, l)[1] for l in lam])
            ess = 150 if idx == low_ess_draw else 600
            pd.DataFrame(dict(di_draw=idx, alpha=r.alpha, ab=r.ab, g=r.g, **{"lambda": r["lambda"]}, Tk=r.Tk,
                              chain=np.repeat([1, 2], inner // 2), iter=np.tile(range(1, inner // 2 + 1), 2),
                              lambda_dd=lam, c_over_alpha=co, divergent=0, treedepth_hits=0, max_rhat=1.003,
                              min_ess_bulk=ess, min_ess_tail=600, min_ebfmi=0.8, gate_pass=True)
                         ).to_csv(dd / f"di_draw_{idx:05d}.csv", index=False)
        json.dump(dict(outer_fits_failing_gate=0, divergences=0, treedepth_hits=0,
                       fitted_onset_below_min=dict(n_draws=0)), open(dd / "summary.json", "w"))
        json.dump(dict(identity_sha256=f"dd-{v}", parent_di_identity=f"di-{v}", di_variant=v,
                       parent_di_draws_sha256=sha(d / "draws_params.csv"), data_file_sha256=sha(data),
                       stan_sha256=sha(stan / "dd_cal.stan"), common_sha256=sha(stan / "common_functions.stanfunctions"),
                       picks_sha256=hashlib.sha256(picks.tobytes()).hexdigest()), open(dd / "identity.json", "w"))
    return src


def run_export(monkeypatch, tmp_path, mutate=None, **kw):
    src = fake_fit(tmp_path, **kw)
    to = tmp_path / "draws"; to.mkdir()
    (to / "sentinel.txt").write_text("previous generation")
    if mutate:
        mutate(tmp_path, src)
    monkeypatch.setenv("BED_EXPORT_FROM", str(src)); monkeypatch.setenv("BED_EXPORT_DD_N", "2")
    monkeypatch.setenv("BED_EXPORT_TO", str(to))
    import export_draws
    E = importlib.reload(export_draws)
    monkeypatch.setattr(E, "ROOT", str(tmp_path))
    E.main()
    return to


def test_dd_rows_are_paired_with_their_di_parent(monkeypatch, tmp_path):
    to = run_export(monkeypatch, tmp_path)
    di = pd.read_csv(to / "tumor_di_primary.csv"); dd = pd.read_csv(to / "tumor_dd_primary.csv")
    assert len(di) == len(dd) == 2 * 4
    assert (di.di_draw == dd.di_draw).all() and np.allclose(di.alpha, dd.alpha) and np.allclose(di.ab, dd.ab)
    assert dd.groupby("di_draw")["lambda"].nunique().min() == 4 and di.groupby("di_draw")["lambda"].nunique().max() == 1
    man = json.load(open(to / "draws_manifest.json"))
    assert man["files"]["tumor_dd_primary.csv"]["parent_di_identity"] == "di-primary"
    assert all(f["diagnostics"]["gate_pass"] is True for f in man["files"].values())
    assert man["gate_policy"]["min_ess_bulk"] == 400 and man["draws_version"] == "3.2"
    assert man["files"]["tumor_dd_primary.csv"]["diagnostics"]["anchor_max_abs_residual_gy"] < 1e-6
    assert not (to / "sentinel.txt").exists()          # the whole folder was replaced
    pv = json.load(open(to / "refit_point_values.json"))["values"]["dd"]
    assert "lambda_over_alpha" in pv and "loss_per_day_after_onset" not in pv


def change_parent_lambda(root, src):
    p = src / "di_primary" / "draws_params.csv"; d = pd.read_csv(p); d.loc[5, "lambda"] += 0.01; d.to_csv(p, index=False)


def change_data_bytes(root, src):
    (root / "data" / DATA_NAME).write_text("y\n")


def change_stan(root, src):
    with open(root / "fit" / "stan" / "dd_cal.stan", "a") as f:
        f.write("\n// edited\n")


def change_variant(root, src):
    p = src / "dd_calibrated_onset7_n2" / "identity.json"; i = json.load(open(p)); i["di_variant"] = "paper_restrictions"
    json.dump(i, open(p, "w"))


def drop_variant_key(root, src):
    p = src / "dd_calibrated_onset7_n2" / "identity.json"; i = json.load(open(p)); del i["di_variant"]; json.dump(i, open(p, "w"))


def negative_g(root, src):
    p = src / "di_primary" / "draws_params.csv"; d = pd.read_csv(p); d["g"] = -0.01; d.to_csv(p, index=False)
    q = src / "dd_calibrated_onset7_n2"
    for f in q.glob("di_draw_*.csv"):
        x = pd.read_csv(f); x["g"] = -0.01; x.to_csv(f, index=False)
    i = json.load(open(q / "identity.json")); i["parent_di_draws_sha256"] = sha(p); json.dump(i, open(q / "identity.json", "w"))


def change_pair_tk(root, src):
    p = src / "dd_calibrated_onset7_n2" / "di_draw_00005.csv"; d = pd.read_csv(p); d["Tk"] = 30.0; d.to_csv(p, index=False)


def break_anchor(root, src):
    p = src / "dd_calibrated_onset7_n2" / "di_draw_00005.csv"; d = pd.read_csv(p); d["c_over_alpha"] *= 1.5; d.to_csv(p, index=False)


@pytest.mark.parametrize("mutate,msg", [
    (change_parent_lambda, "parent DI draws differ"),
    (change_data_bytes, "data file"),
    (change_stan, "Stan source"),
    (change_variant, "variant"),
    (change_pair_tk, "Tk differs from the parent"),
    (break_anchor, "at the anchor"),
    (drop_variant_key, "variant"),
    (negative_g, "negative g"),
])
def test_export_refuses_a_broken_identity_chain_and_leaves_draws_untouched(monkeypatch, tmp_path, mutate, msg):
    with pytest.raises(SystemExit) as e:
        run_export(monkeypatch, tmp_path, mutate=mutate)
    assert "REFUSED" in str(e.value) and msg in str(e.value)
    assert sorted(os.listdir(tmp_path / "draws")) == ["sentinel.txt"]


@pytest.mark.parametrize("kw", [dict(di_pass=False), dict(low_ess_draw=5)])
def test_export_refuses_a_failing_gate(monkeypatch, tmp_path, kw):
    # low_ess_draw: a conditional fit with 0 divergences and R-hat 1.003 but bulk ESS 150 (review R01)
    with pytest.raises(SystemExit) as e:
        run_export(monkeypatch, tmp_path, **kw)
    assert "REFUSED" in str(e.value) and "gate" in str(e.value)


def test_v321_a_held_export_lock_refuses_and_leaves_draws_untouched(monkeypatch, tmp_path):
    """Review F05: two exports must not publish into the same folder at once."""
    (tmp_path / "draws.export.lock").write_text("held")
    with pytest.raises(SystemExit) as e:
        run_export(monkeypatch, tmp_path)
    assert "export.lock" in str(e.value)
    assert (tmp_path / "draws" / "sentinel.txt").read_text() == "previous generation"
    assert not [p for p in tmp_path.iterdir() if ".staging" in p.name]
