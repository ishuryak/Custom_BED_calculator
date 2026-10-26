"""Review R08 (2026-10-26): source intervals and prescription-mixture provenance are kept in a
table beside the fitting file, and the builder refuses ambiguous matches and absent quote checks."""
import csv, importlib, json, os, shutil, sys
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "fit"))
EXTRACT = os.path.join(ROOT, "..", "References_refit_trials", "search", "extraction")
needs_extraction = pytest.mark.skipif(not os.path.isdir(EXTRACT), reason="extraction folder not present")


def mods(monkeypatch, extract=EXTRACT):
    # patch the attribute (restored after the test) rather than reloading the shared module
    import build_refit_dataset as B, build_provenance as P
    monkeypatch.setattr(B, "EXTRACT", str(extract))
    return B, P


@needs_extraction
def test_known_intervals_and_mixtures_survive(monkeypatch):
    B, P = mods(monkeypatch)
    rows = P.provenance(["x", "out.csv"])
    ci = {(r["first_author"], r["arm"]): r["source_ci"] for r in rows if r["source_ci"]}
    assert sorted(ci.values()) == sorted(["53-93", "12-48", "63-93", "15-55", "40-56", "32-47", "39-56", "43-61"])
    for r in rows:
        if not r["source_ci"]:
            assert r["ci_missing_reason"].startswith("not extracted")
    dah = [r for r in rows if r["first_author"] == "DAHANCA67"]
    assert len(dah) == 2 and all("62 Gy in 31 fractions" in r["prescription_categories"] and r["prescription_weights"] == "unknown" for r in dah)
    assert all("nasopharyngeal" in r["population_exceptions"] for r in rows if r["first_author"] == "RTOG8313")
    # the committed provenance file is what the builder produces now
    committed = list(csv.DictReader(open(os.path.join(ROOT, "data", "refit_trials_2026_10_25c_main_provenance.csv"), encoding="utf-8")))
    assert [r["source_ci"] for r in committed] == [r["source_ci"] for r in rows]


def copy_extraction(tmp_path):
    d = tmp_path / "extraction"; d.mkdir()
    for f in ("all_extract.json", "supp_extract_2026_10_24.json"):
        shutil.copy(os.path.join(EXTRACT, f), d / f)
    return d


@needs_extraction
def test_builder_refuses_an_ambiguous_match(monkeypatch, tmp_path):
    d = copy_extraction(tmp_path)
    R = json.load(open(d / "all_extract.json", encoding="utf-8"))
    r = next(x for x in R if x["trial"] == "PMH_HF" and x.get("value_ci"))
    R.append(dict(r, n_analyzed=999))          # a second row with the same value, different N
    json.dump(R, open(d / "all_extract.json", "w", encoding="utf-8"))
    B, P = mods(monkeypatch, d)
    out, problems = B.build(["x", "out.csv"])
    assert any("match with different content" in p for p in problems)


@needs_extraction
def test_builder_refuses_a_missing_quote_check(monkeypatch, tmp_path):
    d = copy_extraction(tmp_path)
    R = json.load(open(d / "all_extract.json", encoding="utf-8"))
    for x in R:
        if x["trial"] == "PMH_HF":
            x["_quote_ok"] = {}
    json.dump(R, open(d / "all_extract.json", "w", encoding="utf-8"))
    B, P = mods(monkeypatch, d)
    out, problems = B.build(["x", "out.csv"])
    assert any("no recorded quote check" in p for p in problems)


@pytest.mark.parametrize("name", ["main", "effn_half", "nominal_clock", "km_only", "rtog9003_beitler", "artscan", "new_trials", "chemo"])
def test_each_provenance_table_lines_up_with_its_data_file(name):
    d = list(csv.DictReader(open(os.path.join(ROOT, "data", f"refit_trials_2026_10_25c_{name}.csv"), encoding="utf-8")))
    p = list(csv.DictReader(open(os.path.join(ROOT, "data", f"refit_trials_2026_10_25c_{name}_provenance.csv"), encoding="utf-8")))
    assert [(r["row"], r["ref_no"], r["arm"]) for r in d] == [(r["row"], r["ref_no"], r["arm"]) for r in p]
