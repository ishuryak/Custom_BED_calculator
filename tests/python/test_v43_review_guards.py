"""Review of 2026-10-31 (v4), F06 and F07: diagnostics are validated per chain and per expected
quantity before aggregation, and an ESS repair cannot take its anchor from mutable current code.

No test compiles or samples Stan: fits are replaced by small synthetic stand-ins.
"""
import hashlib, importlib, itertools, json, os, sys
import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "fit"))
import gate as G

DATA = os.path.join(ROOT, "data", "refit_trials_2026_10_25c_main.csv")
GOOD = dict(divergences=0, treedepth_hits=0, max_rhat=1.001, min_ess_bulk=900, min_ess_tail=800, min_ebfmi=0.9)
N_DRAWS, N_CHAINS = 200, 4


# --- synthetic CmdStanMCMC stand-ins -----------------------------------------------------------

def energy(seed=11):
    """Clean four-chain energy series (E-BFMI near 2 for white noise)."""
    return np.random.default_rng(seed).normal(size=(N_DRAWS, N_CHAINS))


def full_names(sampled, N, J, R):
    return list(sampled) + [f"c[{j}]" for j in range(1, J + 1)] + [f"z[{i}]" for i in range(1, N + 1)] \
        + [f"delta_eval[{r}]" for r in range(1, R + 1)]


def summary_table(names):
    return pd.DataFrame([(n, 0.5, 1.001, 900.0, 800.0) for n in names],
                        columns=["name", "StdDev", "R_hat", "ESS_bulk", "ESS_tail"]).set_index("name")


class FakeFit:
    def __init__(self, table, e=None, divergent=None, treedepth=None, data=None):
        self.table = table
        self.e = energy() if e is None else e
        self.divergent = np.zeros(self.e.shape) if divergent is None else divergent
        self.treedepth = np.full(self.e.shape, 5.0) if treedepth is None else treedepth
        self.data = data
        self.column_names = ["lp__", "energy__"]

    def summary(self):
        return self.table

    def method_variables(self):
        return dict(divergent__=self.divergent, treedepth__=self.treedepth)

    def draws(self, concat_chains=False):
        out = np.zeros(self.e.shape + (2,))
        out[:, :, 1] = self.e
        return out

    def draws_pd(self):
        n = self.e.size
        d = dict(chain__=np.repeat(np.arange(1, self.e.shape[1] + 1), self.e.shape[0]),
                 iter__=np.tile(np.arange(1, self.e.shape[0] + 1), self.e.shape[1]),
                 lambda_dd=np.full(n, 0.5), onset_ref=np.full(n, 20.0), c_over_alpha=np.full(n, 30.0),
                 onset_min_fitted=np.full(n, 15.0), divergent__=np.zeros(n), rep_anchor=np.full(n, self.data["R0"]))
        for i in range(1, self.data["N"] + 1):
            d[f"log_lik[{i}]"] = np.full(n, -1.0)
        for j in range(1, self.data["R"] + 1):
            d[f"delta_eval[{j}]"] = np.zeros(n)
            d[f"onset_eval[{j}]"] = np.full(n, 10.0)
        return pd.DataFrame(d)


SAMPLED_DD = ["lambda_dd", "sigma_raw"]
DIMS = dict(N=5, J=3, R=4)


def diag(table=None, **kw):
    import fit_refit as F
    t = summary_table(full_names(SAMPLED_DD, **DIMS)) if table is None else table
    return F.diagnostics(FakeFit(t, **kw), SAMPLED_DD, dims=DIMS)


# --- F06: E-BFMI per chain ---------------------------------------------------------------------

def test_clean_four_chain_fixture_passes_and_records_constituents():
    dg = diag()
    assert G.passes(dg), G.failures(dg)
    assert len(dg["ebfmi_per_chain"]) == 4 and all(np.isfinite(dg["ebfmi_per_chain"]))
    assert dg["min_ebfmi"] == pytest.approx(min(dg["ebfmi_per_chain"]))
    assert dg["quantities_checked"] == full_names(SAMPLED_DD, **DIMS)
    assert dg["n_quantities"] == 2 + 3 + 5 + 4


@pytest.mark.parametrize("kind", ["constant", "nan", "inf"])
@pytest.mark.parametrize("pos", range(N_CHAINS))
def test_one_bad_energy_chain_fails_in_every_position(kind, pos):
    """Review counterexample: a constant second chain gave E-BFMI 0/0 = NaN and min() kept the finite first value."""
    import fit_refit as F
    e = energy()
    if kind == "constant":
        e[:, pos] = 3.0
    else:
        e[7, pos] = np.nan if kind == "nan" else np.inf
    t = summary_table(full_names(SAMPLED_DD, **DIMS))
    dg = F.diagnostics(FakeFit(t, e=e), SAMPLED_DD)     # also through the base-name call form
    assert not G.passes(dg), (kind, pos, dg["min_ebfmi"])
    # both guards, each on its own: the aggregate is not finite and the chain is named
    assert not np.isfinite(dg["min_ebfmi"]) and any(f"chain {pos + 1}" in p for p in dg["invalid_sampler"])


def test_bad_chain_classification_does_not_depend_on_chain_order():
    import fit_refit as F
    e = energy(); e[:, 1] = 3.0
    t = summary_table(full_names(SAMPLED_DD, **DIMS))
    verdicts = {G.passes(F.diagnostics(FakeFit(t, e=e[:, list(p)]), SAMPLED_DD)) for p in itertools.permutations(range(4))}
    assert verdicts == {False}


def test_ebfmi_helper_reports_the_bad_chain():
    e = energy(); e[:, 2] = 1.0
    vals, problems = G.ebfmi_per_chain(e)
    assert vals[2] is None and any("chain 3" in p for p in problems)
    vals, problems = G.ebfmi_per_chain(energy())
    assert problems == [] and len(vals) == 4


@pytest.mark.parametrize("what", ["divergent_fraction", "divergent_negative", "treedepth_nan", "shape"])
def test_malformed_sampler_arrays_fail(what):
    e = energy()
    div, td = np.zeros(e.shape), np.full(e.shape, 5.0)
    if what == "divergent_fraction":
        div[0, 0], div[1, 0] = 0.5, -0.5            # sums to 0
    elif what == "divergent_negative":
        div[0, 0], div[0, 1] = -1, 1                # sums to 0
    elif what == "treedepth_nan":
        td[3, 2] = np.nan
    else:
        div = np.zeros((e.shape[0] - 10, e.shape[1]))
    assert not G.passes(diag(e=e, divergent=div, treedepth=td))


# --- F06: expected quantity identities ---------------------------------------------------------

def drop(names, pred):
    return [n for n in names if not pred(n)]


@pytest.mark.parametrize("case", ["no_c", "no_z", "one_delta", "sampled", "duplicate", "unexpected_index", "wrong_shape"])
def test_incomplete_or_malformed_quantity_set_fails(case):
    names = full_names(SAMPLED_DD, **DIMS)
    if case == "no_c":
        names = drop(names, lambda n: n.startswith("c["))
    elif case == "no_z":
        names = drop(names, lambda n: n.startswith("z["))
    elif case == "one_delta":
        names = drop(names, lambda n: n == "delta_eval[3]")
    elif case == "sampled":
        names = drop(names, lambda n: n == "sigma_raw")
    elif case == "duplicate":
        names = names + ["z[2]"]
    elif case == "unexpected_index":
        names = names + ["c[4]"]
    else:
        names = names + ["c[1,1]"]
    dg = diag(summary_table(names))
    assert not G.passes(dg), case


def test_invalid_quantities_requires_every_expected_index():
    names = full_names(SAMPLED_DD, **DIMS)
    expected = G.expected_quantities(SAMPLED_DD, **DIMS)
    assert expected == names
    assert G.invalid_quantities(summary_table(names).rename(columns=str.lower).rename(
        columns={"r_hat": "rhat", "stddev": "sd"}), expected=expected) == []
    t = summary_table(drop(names, lambda n: n.startswith("c["))).rename(columns={"R_hat": "rhat", "ESS_bulk": "ess_bulk",
                                                                               "ESS_tail": "ess_tail", "StdDev": "sd"})
    assert G.invalid_quantities(t, expected=expected)[:1] == ["c[1] absent"]


# --- F06: scalar types agree with the R gate ---------------------------------------------------

@pytest.mark.parametrize("key", list(GOOD))
@pytest.mark.parametrize("bad", [False, np.bool_(False), True, "abc", "", b"0"])
def test_bool_and_nonnumeric_diagnostics_are_rejected(key, bad):
    """R's draws_gate_failures rejects logicals; Python used to take False as the count 0."""
    assert not G.passes(dict(GOOD, **{key: bad}))


def test_numeric_strings_still_accepted_as_in_the_shared_fixtures():
    assert G.passes({k: str(v) for k, v in GOOD.items()})
    assert G.passes(dict(GOOD, divergences=np.int64(0), max_rhat=np.float64(1.001)))


# --- F06 at a real call site: the ESS repair's fit loop ----------------------------------------

def load_modules(monkeypatch, out):
    for k, v in dict(BED_FIT_DATA=DATA, BED_FIT_OUT=str(out), DDCAL_NDI="3", DDCAL_DI_VARIANT="primary",
                     BED_FIT_PRIOR_C_SD="3").items():
        monkeypatch.setenv(k, v)
    import fit_refit, fit_dd_calibrated, repair_dd_ess
    importlib.reload(fit_refit); importlib.reload(fit_dd_calibrated)
    Rp = importlib.reload(repair_dd_ess)
    return Rp, Rp.F, Rp.C


class FakeModel:
    def __init__(self, names):
        self.names, self.calls = names, []

    def sample(self, data, **kw):
        self.calls.append(data)
        return FakeFit(summary_table(self.names(data)), data=data)


def di_row():
    return pd.Series(dict(draw=0, alpha=0.04, ab=10.0, g=0.01, **{"lambda": 0.03}, Tk=30.0))


def test_repair_fit_loop_rejects_a_fit_missing_its_nuisance_vectors(monkeypatch, tmp_path):
    Rp, F, C = load_modules(monkeypatch, tmp_path)
    base = F.load()[2]
    full = lambda d: full_names(SAMPLED_DD, d["N"], d["J"], d["R"])
    m = FakeModel(full)
    *_, ok = Rp.fit_draw(m, di_row(), base, 7.0, 0, str(tmp_path / "x"))
    assert ok and len(m.calls) == 1                     # complete fit: accepted at 4x
    m = FakeModel(lambda d: drop(full(d), lambda n: n.startswith(("c[", "z["))))
    *_, ok = Rp.fit_draw(m, di_row(), base, 7.0, 0, str(tmp_path / "y"))
    assert not ok and len(m.calls) == 3                 # every multiplier refused


# --- F07: the repair anchor is the historical one, never the current constant -------------------

def make_folder(monkeypatch, tmp_path, extra=None, drop_keys=()):
    """A legacy (anchor-less) DD folder whose parent DI identity reproduces the real stan-data hash."""
    out = tmp_path / "fitout"
    Rp, F, C = load_modules(monkeypatch, out)
    di_dir = out / "di_primary"; di_dir.mkdir(parents=True)
    pd.DataFrame(dict(draw=range(5), alpha=0.04, ab=10.0, g=0.01, **{"lambda": 0.03}, Tk=30.0)).to_csv(di_dir / "draws_params.csv", index=False)
    real_di = json.load(open(os.path.join(ROOT, "fit", "output_refit_2026_10_25c_main", "di_primary", "identity.json")))
    json.dump(dict(identity_sha256="di-x", data_file_sha256=F.sha_file(DATA), runner_sha256=real_di["runner_sha256"],
                   stan_data_sha256=real_di["stan_data_sha256"]), open(di_dir / "identity.json", "w"))
    dd = out / "dd_calibrated_onset7_n3"; dd.mkdir()
    pick = np.linspace(0, 4, 3).round().astype(int)
    ident = dict(identity_sha256="dd-x", parent_di_identity="di-x", parent_di_draws_sha256=F.sha_file(di_dir / "draws_params.csv"),
                 data_file_sha256=F.sha_file(DATA), stan_sha256=F.sha_file(os.path.join(ROOT, "fit", "stan", "dd_cal.stan")),
                 common_sha256=F.sha_file(os.path.join(ROOT, "fit", "stan", "common_functions.stanfunctions")),
                 runner_sha256="an older runner", fit_refit_sha256=real_di["runner_sha256"],
                 picks_sha256=hashlib.sha256(pick.tobytes()).hexdigest(),
                 prior_c_sd=3.0, onset_min=7.0, onset_scale=1.0, prior_log_lambda=[-0.6931471805599453, 1.5],
                 di_variant="primary", n_di=3)
    ident.update(extra or {})
    for k in drop_keys:
        ident.pop(k)
    json.dump(ident, open(dd / "identity.json", "w"))
    for idx, e in zip(pick, (900, 150, 900)):
        pd.DataFrame(dict(di_draw=[idx] * 2, divergent=0, treedepth_hits=0, max_rhat=1.005, min_ess_bulk=e,
                          min_ess_tail=500, min_ebfmi=0.8)).to_csv(dd / f"di_draw_{idx:05d}.csv", index=False)
    calls = []

    def no_compile(name):
        raise AssertionError("compile attempted after the target check")
    monkeypatch.setattr(F, "model", no_compile)
    monkeypatch.setattr(Rp, "fit_draw", lambda *a: calls.append(a) or (_ for _ in ()).throw(AssertionError("sampled")))
    monkeypatch.setattr(C, "summarize", lambda *a: None)
    return Rp, F, C, dd, calls


def test_changed_current_anchor_is_refused_against_a_legacy_identity(monkeypatch, tmp_path):
    """Review probe: current anchor 30 x 2 Gy / 40 d with matching Stan data [60, 120, 40] was accepted."""
    Rp, F, C, dd, calls = make_folder(monkeypatch, tmp_path)
    monkeypatch.setattr(F, "ANCHOR", (30, 2.0, 40.0))
    monkeypatch.setattr(C, "ANCHOR", (30, 2.0, 40.0))
    assert F.load()[2]["anchor_dsum"] == 60.0          # the current Stan data really moved
    with pytest.raises(SystemExit) as e:
        Rp.main()
    assert "REFUSED" in str(e.value) and "anchor" in str(e.value) and calls == []


def test_changed_current_evaluation_regimen_is_refused_against_a_legacy_identity(monkeypatch, tmp_path):
    Rp, F, C, dd, calls = make_folder(monkeypatch, tmp_path)
    regs = list(F.REGIMENS); regs[0] = ("HYPNO 55 Gy/20 (28 d)", 20, 2.75, 28)
    monkeypatch.setattr(F, "REGIMENS", regs)
    with pytest.raises(SystemExit) as e:
        Rp.main()
    assert "REFUSED" in str(e.value) and "eval_regimens" in str(e.value) and calls == []


@pytest.mark.parametrize("mutate", ["prior_c_sd", "prior_log_lambda", "onset_scale", "onset_min", "variant", "parent"])
def test_each_changed_current_setting_is_refused_individually(monkeypatch, tmp_path, mutate):
    Rp, F, C, dd, calls = make_folder(monkeypatch, tmp_path,
                                      extra=dict(di_variant="paper_restrictions") if mutate == "variant" else
                                      dict(parent_di_identity="another parent") if mutate == "parent" else None)
    if mutate == "prior_c_sd":
        monkeypatch.setattr(F, "PRIOR_C_SD", 10.0)
    elif mutate == "prior_log_lambda":
        monkeypatch.setattr(C, "PRIOR_LOG_LAMBDA", (np.log(0.5), 1.0))
    elif mutate == "onset_scale":
        monkeypatch.setattr(C, "ONSET_SCALE", 2.0)
    elif mutate == "onset_min":
        monkeypatch.setattr(Rp, "folder", lambda: (str(dd), 3, 5.0, "primary"))
    with pytest.raises(SystemExit) as e:
        Rp.main()
    assert "REFUSED" in str(e.value) and calls == []


def test_legacy_identity_without_verifiable_parent_provenance_is_refused(monkeypatch, tmp_path):
    Rp, F, C, dd, calls = make_folder(monkeypatch, tmp_path, drop_keys=("fit_refit_sha256",))
    with pytest.raises(SystemExit) as e:
        Rp.main()
    assert "REFUSED" in str(e.value) and "legacy" in str(e.value) and calls == []


def test_new_schema_identity_must_state_its_anchor(monkeypatch, tmp_path):
    Rp, F, C, dd, calls = make_folder(monkeypatch, tmp_path, extra=dict(target_schema="dd-target-2026-10-31"))
    assert C.TARGET_SCHEMA == "dd-target-2026-10-31"
    with pytest.raises(SystemExit) as e:
        Rp.main()
    assert "REFUSED" in str(e.value) and "anchor" in str(e.value) and calls == []


def test_unchanged_legacy_target_proceeds_to_repair(monkeypatch, tmp_path):
    Rp, F, C, dd, calls = make_folder(monkeypatch, tmp_path)
    with pytest.raises(AssertionError, match="compile attempted"):
        Rp.main()                                       # target accepted; the run reaches the compile step


def test_fit_draw_calibrates_at_the_checked_anchor_not_the_current_constant(monkeypatch, tmp_path):
    Rp, F, C = load_modules(monkeypatch, tmp_path)
    base = F.load()[2]                                  # checked target: 70 / 140 / 47
    monkeypatch.setattr(C, "ANCHOR", (30, 2.0, 40.0))
    m = FakeModel(lambda d: full_names(SAMPLED_DD, d["N"], d["J"], d["R"]))
    Rp.fit_draw(m, di_row(), base, 7.0, 0, str(tmp_path / "x"))
    data = m.calls[0]
    assert data["T0"] == 47.0 and data["K0"] == pytest.approx(70.0 + 140.0 / 10.0)
    assert data["R0"] == pytest.approx(float(C.rep_di(0.03, 30.0, 47.0)))


def test_new_dd_identity_records_anchor_and_evaluation_regimens(monkeypatch, tmp_path):
    out = tmp_path / "fitout"
    Rp, F, C = load_modules(monkeypatch, out)
    di_dir = out / "di_primary"; di_dir.mkdir(parents=True)
    pd.DataFrame(dict(draw=range(5), alpha=0.04, ab=10.0, g=0.01, **{"lambda": 0.03}, Tk=30.0)).to_csv(di_dir / "draws_params.csv", index=False)
    json.dump(dict(identity_sha256="di-x", data_file_sha256=F.sha_file(DATA)), open(di_dir / "identity.json", "w"))
    monkeypatch.setattr(F, "versions", lambda: {})

    class Stop(Exception):
        pass

    def stop(name):
        raise Stop()
    monkeypatch.setattr(F, "model", stop)
    with pytest.raises(Stop):
        C.main()
    ident = json.load(open(out / "dd_calibrated_onset7_n3" / "identity.json"))
    assert ident["target_schema"] == C.TARGET_SCHEMA
    assert ident["anchor"] == [70.0, 140.0, 47.0] and ident["anchor_regimen"] == [35, 2.0, 47.0]
    assert len(ident["eval_regimens"]) == len(F.REGIMENS) and ident["eval_regimens"][0] == [55.0, 151.25, 26]


# --- F07: every real intended repair target still passes ----------------------------------------

JOBS = [("main", "main", 1000, "primary", 3), ("main", "main", 1000, "paper_restrictions", 3),
        ("effn_half", "effn_half", 200, "primary", 3), ("km_only", "km_only", 200, "primary", 3),
        ("rtog9003_beitler", "rtog9003_beitler", 200, "primary", 3), ("artscan", "artscan", 200, "primary", 3),
        ("chemo", "chemo", 200, "primary", 3), ("new_trials", "new_trials", 200, "primary", 3),
        ("nominal_clock", "nominal_clock", 200, "primary", 3), ("prior_c10", "main", 200, "primary", 10)]


@pytest.mark.parametrize("job", JOBS, ids=[f"{j[0]}_{j[3]}" for j in JOBS])
def test_real_intended_targets_pass_the_target_check(monkeypatch, job):
    name, data, ndi, variant, prior = job
    monkeypatch.setenv("BED_FIT_DATA", os.path.join(ROOT, "data", f"refit_trials_2026_10_25c_{data}.csv"))
    monkeypatch.setenv("BED_FIT_PRIOR_C_SD", str(prior))
    import fit_refit, fit_dd_calibrated, repair_dd_ess
    importlib.reload(fit_refit); importlib.reload(fit_dd_calibrated); Rp = importlib.reload(repair_dd_ess)
    folder = os.path.join(ROOT, "fit", f"output_refit_2026_10_25c_{name}")
    sub = f"dd_calibrated_onset7_n{ndi}" + ("" if variant == "primary" else f"_{variant}")
    ident = json.load(open(os.path.join(folder, sub, "identity.json")))
    di_ident = json.load(open(os.path.join(folder, f"di_{variant}", "identity.json")))
    assert "anchor" not in ident                        # legacy: resolved, never filled from current code
    base = Rp.F.load()[2]
    Rp.check_target(ident, Rp.posterior_target(base, 7.0, variant, ndi), di_ident, base)
