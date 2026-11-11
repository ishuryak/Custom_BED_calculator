"""Applicability table (review F3): schedule classes follow the stated rule on known trial arms."""
import importlib.util
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_spec = importlib.util.spec_from_file_location("app", os.path.join(ROOT, "tools", "early_tox", "applicability_table.py"))
app = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(app)


def test_schedule_classes_of_known_arms():
    rows = {r["arm_id"]: r["schedule_class"] for r in app.build("/nonexistent")}
    expected = {
        "RTOG9003_SF": "conventional", "HYPNO_HFX": "hypofractionated", "HYPNO_NFX": "accelerated",
        "ELAN_HSC": "split course", "RTOG9003_HFX": "hyperfractionated", "Vienna_VCHART": "accelerated",
        "PMH_CF": "hypofractionated", "PMH_HF": "hyperfractionated", "Tandon_SMART": "hypofractionated",
        "Ghoshal_CB": "accelerated", "GuptaM_6F": "accelerated", "Grover_SIB": "hypofractionated",
    }
    assert {k: rows[k] for k in expected} == expected


def test_every_primary_and_external_arm_is_listed_once():
    rows = app.build("/nonexistent")
    ids = [r["arm_id"] for r in rows]
    assert len(ids) == len(set(ids)) == 66
    assert sum(r["role"].startswith("primary") for r in rows) == 51
