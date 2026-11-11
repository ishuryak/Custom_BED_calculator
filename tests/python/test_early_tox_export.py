"""Export of the fitted early-toxicity model for the app's research mode."""
import hashlib
import importlib.util
import json
import os
import shutil

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_spec = importlib.util.spec_from_file_location("x", os.path.join(ROOT, "tools", "early_tox", "export_fit_for_app.py"))
x = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(x)


def _unit(tmp_path, **change):
    u = tmp_path / "U02"
    shutil.copytree(x.UNIT_DEFAULT, u)
    m = json.loads((u / "manifest.json").read_text())
    m.update(change)
    (u / "manifest.json").write_text(json.dumps(m))
    return str(u)


def test_exported_draws_match_their_manifest():
    m = json.load(open(os.path.join(ROOT, "draws", "early_fit_manifest.json")))
    p = os.path.join(ROOT, "draws", "early_fit_draws.csv")
    assert hashlib.sha256(open(p, "rb").read()).hexdigest() == m["files"]["early_fit_draws.csv"]["sha256"]
    assert m["n_draws"] == 4000 and sum(1 for _ in open(p)) == 4001
    src = json.load(open(os.path.join(x.UNIT_DEFAULT, "manifest.json")))
    assert m["source_identity"] == src["identity"] and m["model_id"].endswith(src["identity"][:12])
    assert m["fixed"]["repair_halftime_h"] == 3.2 and m["fixed"]["Tk_e_elapsed_days"] == 7.0


def test_failed_smoke_or_changed_units_are_refused(tmp_path):
    with pytest.raises(RuntimeError, match="gate did not pass"):
        x.read_unit(_unit(tmp_path / "a", gate_pass=False))
    with pytest.raises(RuntimeError, match="smoke run"):
        x.read_unit(_unit(tmp_path / "b", smoke=True))
    u = _unit(tmp_path / "c")
    with open(os.path.join(u, "draws.csv.gz"), "ab") as f:
        f.write(b"x")
    with pytest.raises(RuntimeError, match="changed since the run"):
        x.read_unit(u)


def test_trial_table_is_hashed_cited_and_one_control_arm_per_trial():
    import csv
    m = json.load(open(os.path.join(ROOT, "draws", "early_fit_manifest.json")))
    p = os.path.join(ROOT, "draws", "early_fit_trials.csv")
    assert hashlib.sha256(open(p, "rb").read()).hexdigest() == m["files"]["early_fit_trials.csv"]["sha256"]
    rows = list(csv.DictReader(open(p)))
    fitted = [r for r in rows if r["role"] == "fitted"]
    assert len({r["trial"] for r in fitted}) == 22 and sum(int(r["patients"]) for r in fitted) == 6042
    assert all(r["pmid"].isdigit() for r in rows)
    trials = {r["trial"] for r in rows}
    assert sum(int(r["control_arm"]) for r in rows) == len(trials)


# Review R4 (v4.6): only the approved primary profile can be exported; failed or altered sources refuse.
@pytest.mark.parametrize("unit,field", [("U04", "settings"), ("U17", "settings"), ("U19", "settings"), ("U15", "settings")])
def test_other_scientific_profiles_cannot_be_exported_as_the_primary(unit, field):
    u = os.path.join(os.path.dirname(x.UNIT_DEFAULT), unit)
    with pytest.raises(RuntimeError, match=f"not the approved primary profile \\({field}"):
        x.read_unit(u)


def _rehash(u, fn):
    m = json.loads(open(os.path.join(u, "manifest.json")).read())
    m["outputs"][fn] = x.sha256(os.path.join(u, fn))
    open(os.path.join(u, "manifest.json"), "w").write(json.dumps(m))


def test_altered_diagnostics_or_types_refuse_even_with_matching_hashes(tmp_path):
    u = _unit(tmp_path / "d")
    d = json.loads(open(os.path.join(u, "diagnostics.json")).read())
    d.update(gate_pass=False, rhat=2.0, divergences=99)
    open(os.path.join(u, "diagnostics.json"), "w").write(json.dumps(d))
    _rehash(u, "diagnostics.json")
    with pytest.raises(RuntimeError, match="do not meet the gate"):
        x.read_unit(u)
    with pytest.raises(RuntimeError, match="gate_pass is not a boolean"):
        x.read_unit(_unit(tmp_path / "s", gate_pass="true"))
    u2 = _unit(tmp_path / "m")
    os.remove(os.path.join(u2, "ppc_arms.csv"))
    with pytest.raises(RuntimeError, match="changed since the run"):
        x.read_unit(u2)


def test_a_staged_package_that_does_not_match_its_manifest_publishes_nothing(tmp_path):
    stage = tmp_path / "stage"
    stage.mkdir()
    for fn in ("early_fit_draws.csv", "early_fit_trials.csv", "early_fit_evaluator.json"):
        (stage / fn).write_text("x")
    (stage / "early_fit_manifest.json").write_text(json.dumps(
        {"files": {"early_fit_draws.csv": {"sha256": "0" * 64}, "early_fit_trials.csv": {"sha256": "0" * 64}}}))
    root = tmp_path / "root"
    (root / "draws").mkdir(parents=True)
    (root / "tests" / "testthat" / "fixtures").mkdir(parents=True)
    (root / "draws" / "early_fit_manifest.json").write_text("previous")
    with pytest.raises(RuntimeError, match="does not match its manifest; nothing published"):
        x.publish_package(str(stage), root=str(root))
    assert (root / "draws" / "early_fit_manifest.json").read_text() == "previous"
    assert not (root / "draws" / "early_fit_draws.csv").exists()
