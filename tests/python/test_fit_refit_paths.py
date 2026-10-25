"""fit_refit.py honours BED_FIT_DATA / BED_FIT_OUT and keeps the app's defaults without them.

Run: python3 -m pytest tests/python -q  (from the repository root)"""
import importlib, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "fit"))


def load(monkeypatch, **env):
    for k in ("BED_FIT_DATA", "BED_FIT_OUT"):
        monkeypatch.delenv(k, raising=False)
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    import fit_refit
    return importlib.reload(fit_refit)


def test_defaults_unchanged(monkeypatch):
    m = load(monkeypatch)
    assert m.DATA == os.path.join(ROOT, "data", "shuryak2018_table1.csv")
    assert m.OUT == os.path.join(ROOT, "fit", "output")


def test_overrides_reach_module(monkeypatch, tmp_path):
    data = tmp_path / "d.csv"; out = tmp_path / "o"
    m = load(monkeypatch, BED_FIT_DATA=str(data), BED_FIT_OUT=str(out))
    assert m.DATA == str(data) and m.OUT == str(out)
