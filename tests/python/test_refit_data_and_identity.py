"""Data corrections and fit-runner contracts from the 2026-10-25 review (plan stages B and C).

Run: python3 -m pytest tests/python -q  (from the repository root)
The builder tests need the extraction folder (BED_EXTRACT_DIR or ../References_refit_trials/...);
they are skipped when it is absent. BED_TEST_DATA points the dataset tests at another CSV, which
is how they were shown to fail on the 2026-10-25 file before the corrections."""
import csv, importlib, json, os, sys
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "fit"))
import build_refit_dataset as B

HAVE_EXTRACT = os.path.exists(os.path.join(B.EXTRACT, "all_extract.json"))
needs_extract = pytest.mark.skipif(not HAVE_EXTRACT, reason="extraction folder not available")


@pytest.fixture(scope="module")
def main_rows(tmp_path_factory):
    path = os.environ.get("BED_TEST_DATA")
    if not path:
        if not HAVE_EXTRACT:
            pytest.skip("extraction folder not available")
        path = str(tmp_path_factory.mktemp("d") / "main.csv")
        B.main(["x", path])
    return list(csv.DictReader(open(path, encoding="utf-8")))


# ---- data corrections ---------------------------------------------------------------

def test_krog_hypo_counts_nine_t2_patients(main_rows):
    # Moon 2014 Table 1: HYPO T1a 45 + T1b 20 = 65, T2a 7 + T2b 2 = 9, total 74.
    r = [x for x in main_rows if x["ref_no"] == "KROG0201" and x["arm"].startswith("HYPO")][0]
    assert float(r["n_fractions"]) == pytest.approx((65 * 28 + 9 * 30) / 74, abs=5e-4)
    assert int(r["n_patients"]) == 74


def test_rtog9003_rate_and_n_from_the_same_report(main_rows):
    r = [x for x in main_rows if x["ref_no"] == "RTOG9003"]
    assert {x["pmid"] for x in r} == {"10924966"}
    assert sorted(float(x["tcp_actual_pct"]) for x in r) == [46.0, 47.5, 54.4, 54.5]
    assert sum(int(x["n_patients"]) for x in r) == 1073


def test_mixed_fraction_sizes_kept_as_components(main_rows):
    g = [x for x in main_rows if x["ref_no"] == "Ghoshal" and x["arm"].startswith("Concomitant")][0]
    assert float(g["sum_d2_gy2"]) == pytest.approx(25 * 1.8 ** 2 + 15 * 1.5 ** 2)
    assert float(g["total_dose_gy"]) == pytest.approx(67.5)
    # mean-fraction BED would understate the quadratic term
    assert float(g["sum_d2_gy2"]) > 40 * (67.5 / 40) ** 2


def test_source_rate_estimator_and_both_clocks_are_carried(main_rows):
    for x in main_rows:
        assert x["estimate_type"] in ("actuarial/KM", "cumulative incidence")
        assert x["T_nominal_days"] != "" and x["T_source_days"] != ""
        v = float(x["source_value_pct"])
        tcp = 100 - v if x["source_value_is_failure"] == "True" else v
        assert float(x["tcp_actual_pct"]) == pytest.approx(tcp)


# ---- builder mechanics --------------------------------------------------------------

@pytest.mark.parametrize("F,week,gap,expected", [
    (35, B.MON_FRI, None, 47), (33, B.MON_SAT, None, 38), (20, B.MON_FRI, None, 26),
    (35, B.MON_SAT, None, 40), (20, B.MON_FRI, (10, 14), 40),        # ELAN split course, planned 40 d
    (42, B.BID_MON_FRI, (24, 14), 43),                                 # RTOG 90-03 / 88-09 split arm
    (36, (3,) * 7, None, 12), (33, (2, 2, 2, 2, 1, 0, 0), None, 24),   # CHART, TROG ART
])
def test_calendar(F, week, gap, expected):
    assert B.last_day(F, week, gap) == expected


@needs_extract
def test_equal_component_split_changes_nothing(monkeypatch):
    base = {o["ref_no"] + o["arm"]: o for o in B.build(["x"])[0]}
    monkeypatch.setitem(B.COMPONENTS, ("TROG9101", "CRT"), [(20, 2.0), (15, 2.0)])
    split = {o["ref_no"] + o["arm"]: o for o in B.build(["x"])[0]}
    k = next(k for k in base if k.startswith("TROG9101CRT"))
    assert split[k]["sum_d2_gy2"] == pytest.approx(base[k]["sum_d2_gy2"]) == pytest.approx(35 * 4.0)


@needs_extract
def test_builder_refuses_unexplained_patient_number(monkeypatch):
    bad = [a for a in B.A]
    i = next(k for k, a in enumerate(bad) if a[0] == "TROG9101")
    bad[i] = bad[i][:7] + (999,) + bad[i][8:]
    monkeypatch.setattr(B, "A", bad)
    _, problems = B.build(["x"])
    assert any("TROG9101" in p and "999" in p for p in problems)


@needs_extract
def test_builder_refuses_components_that_do_not_sum(monkeypatch):
    monkeypatch.setitem(B.COMPONENTS, ("Ghoshal", "Concomitant boost"), [(25, 1.8), (14, 1.5)])
    _, problems = B.build(["x"])
    assert any("Ghoshal" in p and "components" in p for p in problems)


# ---- fit runner ---------------------------------------------------------------------

def _fit_refit(monkeypatch, **env):
    for k in ("BED_FIT_DATA", "BED_FIT_OUT", "BED_FIT_PRIOR_C_SD"):
        monkeypatch.delenv(k, raising=False)
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    import fit_refit
    return importlib.reload(fit_refit)


def test_load_uses_fraction_components(monkeypatch, tmp_path):
    p = tmp_path / "d.csv"
    p.write_text("row,first_author,ref_no,arm,dose_per_fraction_gy,n_fractions,total_time_days,n_patients,"
                 "tcp_actual_pct,total_dose_gy,sum_d2_gy2\n"
                 "1,G,G,boost,1.6875,40,36,143,73.6,67.5,114.75\n"
                 "2,G,G,conv,2.0,33,46,142,54.5,66.0,132.0\n")
    F = _fit_refit(monkeypatch, BED_FIT_DATA=str(p))
    _, _, base = F.load()
    assert base["d2sum"][0] == pytest.approx(114.75)          # not 40 x 1.6875^2 = 113.9
    assert base["dsum"][0] == pytest.approx(67.5)


def test_load_without_component_columns_matches_n_times_d(monkeypatch):
    F = _fit_refit(monkeypatch)                                # 2018 Table 1, no component columns
    rows, _, base = F.load()
    r = rows[0]
    n, d = float(r["n_fractions"]), float(r["dose_per_fraction_gy"])
    assert base["dsum"][0] == pytest.approx(n * d) and base["d2sum"][0] == pytest.approx(n * d * d)


def test_finished_fit_with_unknown_identity_is_not_reused(monkeypatch, tmp_path):
    # Review M5: sentinel results from an earlier run must not be skipped and relabeled.
    out = tmp_path / "o"
    for name in ("di", "dd"):
        for v in ("primary", "binomial", "paper_restrictions", "wide_ab"):
            (out / f"{name}_{v}").mkdir(parents=True)
            (out / f"{name}_{v}" / "summary.json").write_text('{"sentinel": true}')
    F = _fit_refit(monkeypatch, BED_FIT_OUT=str(out))
    monkeypatch.setattr(F, "model", lambda name: object())     # no compilation, no sampling
    monkeypatch.setattr(sys, "argv", ["fit_refit.py"])
    with pytest.raises(SystemExit) as e:
        F.main()
    assert "REFUSED" in str(e.value)
    assert not (out / "provenance.json").exists()


def test_identity_changes_with_every_component(monkeypatch):
    F = _fit_refit(monkeypatch)
    _, _, base = F.load()
    data = F.stan_data(base, "primary")
    ref = F.identity("di", data, 1, 4, 10, 10)["identity_sha256"]
    assert F.identity("di", data, 2, 4, 10, 10)["identity_sha256"] != ref
    assert F.identity("di", dict(data, prior_c_sd=10.0), 1, 4, 10, 10)["identity_sha256"] != ref
    assert F.identity("di", data, 1, 4, 11, 10)["identity_sha256"] != ref
    assert F.identity("dd", data, 1, 4, 10, 10)["identity_sha256"] != ref
    assert F.identity("di", data, 1, 4, 10, 10)["identity_sha256"] == ref


def test_matching_identity_is_reused(monkeypatch, tmp_path):
    F = _fit_refit(monkeypatch)
    ident = {"identity_sha256": "abc"}
    (tmp_path / "summary.json").write_text("{}")
    (tmp_path / "identity.json").write_text(json.dumps(ident))
    assert F.check_identity(str(tmp_path), ident) is True
    assert F.check_identity(str(tmp_path / "missing"), ident) is False
