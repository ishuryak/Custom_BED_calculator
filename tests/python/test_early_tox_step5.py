"""Early-toxicity step 5: data assembly for the refit, the C1 chemotherapy coding and the Stan model."""
import csv
import importlib.util
import math
import os

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, "tools", "early_tox", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


s5 = load("step5_fit")
ev = s5.ev


def test_chemo_coding_matches_the_planned_therapy_of_every_arm():
    # A nonzero code exactly where an arm's planned systemic therapy differs from its reference arm's.
    for m4 in (False, True):
        rows, fams = s5.load_rows(m4=m4)
        ref = {f: next(x for x in rows if x["arm_id"] == s5.REF_ARM[f]) for f in fams}
        for x in rows:
            differs = x["chemo_text"] != ref[x["family"]]["chemo_text"] and x["arm_id"] not in s5.CHEMO_SAME_DOSE
            assert (s5.CHEMO_CUMULATIVE.get(x["arm_id"], 0) != 0) == differs, x["arm_id"]
    assert {a.split("_")[0] for a in s5.CHEMO_CUMULATIVE} <= {
        "GORTEC9902", "RTOG0129", "Vienna", "Olmi", "Bartelink", "GhoshLaskar", "GuptaP", "Sanchiz"}
    assert not set(s5.CHEMO_SAME_DOSE) & set(s5.CHEMO_CUMULATIVE)
    assert set(s5.CHEMO_ONOFF) <= set(s5.CHEMO_CUMULATIVE)
    assert all(s5.CHEMO_ONOFF[k] == s5.CHEMO_CUMULATIVE[k] for k in s5.CHEMO_ONOFF)


def test_analysis_sets_have_the_approved_sizes():
    rows, fams = s5.load_rows()
    assert (len(rows), len(fams)) == (51, 22)
    rows, fams = s5.load_rows(w1_only=True)
    assert (len(rows), len(fams)) == (32, 13)
    rows, fams = s5.load_rows(drop_chemo=True)
    assert len(fams) == 15 and not set(fams) & s5.CHEMO_FAMILIES
    rows, fams = s5.load_rows(m4=True)
    assert (len(rows), len(fams)) == (64, 28)


def test_stan_contrasts_reproduce_the_step4_evaluator():
    rows, fams = s5.load_rows()
    d = s5.stan_data(rows, fams, 7.0, {})
    ab, r = 10.0, ev.fowler_rate()
    dE = np.array(d["dD"]) + np.array(d["dQ"]) / ab - r * np.array(d["dTe"])
    rows4, fams4 = s5.s4.load_design()
    dE4, _, _ = s5.s4.contrasts(rows4, fams4, ab, r, 7.0)
    assert [x["arm_id"] for x in rows] == [x["arm_id"] for x in rows4]
    assert np.allclose(dE, dE4, atol=1e-9)
    refs = [i for i, x in enumerate(rows) if x["arm_id"] == s5.REF_ARM[x["family"]]]
    assert len(refs) == 22
    assert all(d["dD"][i] == d["dQ"][i] == d["dTe"][i] == d["chemo"][i] == 0 for i in refs)
    # ORO 93-01: identical schedules, chemotherapy differs
    oro = [i for i, x in enumerate(rows) if x["arm_id"] == "Olmi_CFCT"][0]
    assert dE[oro] == 0 and d["chemo"][oro] == 1


def test_scale_groups_and_decision_schedules():
    rows, fams = s5.load_rows(m4=True)
    d = s5.stan_data(rows, fams, 7.0, {})
    grp = np.array(d["grp"])
    J = len(fams)
    assert d["G"] == 2 and grp.shape == (J, 2)
    assert {fams[j] for j in range(J) if grp[j, 0]} == {"Tandon", "Sikdar", "Shivananjappa", "Poddar", "Grover"}
    assert {fams[j] for j in range(J) if grp[j, 1]} == {"Pinto", "Sanchiz", "Dwivedi"}
    # decision contrasts equal early BED differences from the evaluator (3.2 h repair)
    for Tk in (7.0, 21.0):
        d = s5.stan_data(*s5.load_rows(), Tk, {})
        ref = ev.early_bed(*ev.regular_schedule(35, 2.0), ab=10, Tk=Tk, halftime_h=3.2)
        for k, (n, dd) in enumerate(s5.DECISION):
            want = ev.early_bed(*ev.regular_schedule(n, dd), ab=10, Tk=Tk, halftime_h=3.2) - ref
            got = d["kD"][k] + d["kQ"][k] / 10 - ev.fowler_rate() * d["kTe"][k]
            assert got == pytest.approx(want, abs=1e-9)


@pytest.fixture(scope="module")
def model():
    pytest.importorskip("cmdstanpy")
    return s5.compile_model()


def test_stan_model_recovers_known_parameters_from_large_noise_free_counts(model):
    # Expected counts at 200 times the real arm sizes, from known values: kappa 0.1, gamma 1.0,
    # alpha/beta 10, Fowler r. The fit must land near them, so a sign or scale error fails.
    rows, fams = s5.load_rows()
    d = s5.stan_data(rows, fams, 7.0, {})
    eta = np.linspace(-1.0, 0.5, len(fams))
    dE = np.array(d["dD"]) + np.array(d["dQ"]) / 10 - ev.fowler_rate() * np.array(d["dTe"])
    lin = eta[np.array(d["fam"]) - 1] + 0.1 * dE + 1.0 * np.array(d["chemo"])
    N = np.array(d["N"]) * 200
    d.update(N=N.tolist(), y=np.round(N / (1 + np.exp(-lin))).astype(int).tolist())
    fit = model.sample(data=d, chains=2, parallel_chains=1, iter_warmup=300, iter_sampling=300,
                       seed=1, show_progress=False, show_console=False)
    med = {k: float(np.median(fit.stan_variable(k))) for k in ("kappa", "gamma", "r", "ab")}
    assert med["kappa"] == pytest.approx(0.1, rel=0.1)
    assert med["gamma"] == pytest.approx(1.0, abs=0.1)
    assert med["r"] == pytest.approx(ev.fowler_rate(), rel=0.15)
    assert med["ab"] == pytest.approx(10.0, rel=0.3)


def test_overdispersion_switch_reaches_the_fit_and_finished_units_are_skipped(model, tmp_path, monkeypatch):
    monkeypatch.setattr(s5, "wait_for_memory", lambda: 99.0)
    import json
    for od, out in ((True, tmp_path / "od"), (False, tmp_path / "bin")):
        s5.run_unit(model, "U02", str(out), smoke=True, od=od)
        diag = json.loads((out / "U02" / "diagnostics.json").read_text())
        summ = (out / "U02" / "summary.csv").read_text()
        assert diag["settings"].get("use_od", 0) == int(od)
        assert ("sigma_arm[1]" in summ) == od
    stamp = os.path.getmtime(tmp_path / "od" / "U02" / "diagnostics.json")
    s5.run_unit(model, "U02", str(tmp_path / "od"), smoke=True, od=True)
    assert os.path.getmtime(tmp_path / "od" / "U02" / "diagnostics.json") == stamp


def test_leave_one_family_out_removes_exactly_that_family():
    assert len(s5.LOFO_UNITS) == 22
    base, _ = s5.load_rows()
    for name, (_, st) in s5.LOFO_UNITS.items():
        rows, fams = s5.load_rows(drop_family=st["drop_family"])
        assert st["drop_family"] not in fams and len(fams) == 21
        assert len(rows) == len(base) - sum(1 for x in base if x["family"] == st["drop_family"])


def test_sbc_prior_draws_follow_the_registered_priors():
    sbc = load("step5_sbc")
    d = s5.stan_data(*s5.load_rows(), 7.0, dict(use_od=1))
    rng = np.random.default_rng(0)
    th = [sbc.prior_draw(rng, d)[0] for _ in range(4000)]
    col = lambda k: np.array([t[k] for t in th])
    assert np.std(np.log(col("ab"))) == pytest.approx(0.354, rel=0.06)
    assert np.mean(np.log(col("ab"))) == pytest.approx(math.log(10), abs=0.03)
    assert np.mean(np.log(col("r"))) == pytest.approx(math.log(0.8), abs=0.04)
    assert np.mean(col("Lstar")) == pytest.approx(math.sqrt(2 / math.pi), rel=0.06)     # HalfNormal(1)
    assert np.mean(col("sigma_arm")) == pytest.approx(0.5 * math.sqrt(2 / math.pi), rel=0.06)
    y = sbc.prior_draw(np.random.default_rng(1), d)[1]
    assert len(y) == d["I"] and all(0 <= a <= n for a, n in zip(y, d["N"]))


def test_an_arm_entering_the_analysis_twice_is_refused(monkeypatch):
    real = s5.read

    def doubled(name):
        rows = real(name)
        if name == "endpoint_mapping.csv":
            rows = rows + [dict(r, outcome_id="O055") for r in rows if r["outcome_id"] == "O055"]
        return rows
    monkeypatch.setattr(s5, "read", doubled)
    with pytest.raises(ValueError, match=r"entered twice: \[.GuptaP_CFCCT.\]"):
        s5.load_rows()


def test_adequate_ascertainment_unit_loads_only_the_strict_set():
    rows, fams = s5.load_rows(adequate_only=True)
    assert fams == ["CAIR1", "ORO9301"] and len(rows) == 5
    assert s5.UNITS["U15"][1] == dict(adequate_only=True)


def _cum(row):
    v = (row or {}).get("planned_cumulative_mg_m2", "").split(" ")[0]
    return float(v) if v else None


def test_planned_therapy_differences_are_coded_or_declared_as_scenarios():
    # Exposure completeness (review F1): a known planned-dose difference against the reference arm
    # must carry a matching chemotherapy code (fitted arms) or declared scenarios (external arms).
    pt = s5.planned_therapy()
    arms = {a["arm_id"]: a for a in s5.read("arms.csv")}
    ref_of = dict(s5.REF_ARM, HYPNO="HYPNO_NFX")
    checked = 0
    for aid, a in arms.items():
        ref = ref_of.get(a["family_id"])
        if ref is None or aid == ref or aid not in pt or ref not in pt:
            continue
        if pt[aid]["status"] == "assumed_equal":
            assert s5.CHEMO_CUMULATIVE.get(aid, 0) == 0, aid
            continue
        x, y = _cum(pt[aid]), _cum(pt[ref])
        if x is None or y is None or x == y:
            continue
        checked += 1
        if aid in s5.CHEMO_CUMULATIVE:
            assert (s5.CHEMO_CUMULATIVE[aid] > 0) == (x > y), aid
        else:
            assert aid in s5.EXTERNAL_CHEMO_SCENARIOS, f"{aid}: planned {x} vs {y} but no code or scenario"
    assert checked >= 5


def test_hypno_scenarios_and_tandon_sensitivity_are_declared():
    sc = s5.EXTERNAL_CHEMO_SCENARIOS["HYPNO_HFX"]
    assert sc["code_0"] == 0 and sc["code_minus1"] == -1 and sc["mixture_75pct"] == ("mixture", 0.75, -1.0)
    pt = s5.planned_therapy()
    assert (pt["HYPNO_HFX"]["planned_courses"], pt["HYPNO_NFX"]["planned_courses"]) == ("4", "5")
    assert pt["Tandon_SMART"]["status"] == "assumed_equal" and s5.CHEMO_SENSITIVITY == {"Tandon_SMART": -1}


def test_reported_totals_come_from_the_loaded_rows():
    rows, _ = s5.load_rows()
    assert s5.analysis_totals(rows) == dict(arms=51, families=22, evaluable=6042, events=2945)
    with pytest.raises(ValueError, match="counted twice"):
        s5.analysis_totals(rows + rows[:1])
