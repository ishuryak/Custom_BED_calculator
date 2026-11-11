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
