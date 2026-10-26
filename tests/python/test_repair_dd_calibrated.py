"""fit/repair_dd_calibrated.py is retired (review R01, 2026-10-26): its gate omitted ESS, so it
must refuse to run rather than select fits by the narrower rule. fit/repair_dd_ess.py replaces it
(tests/python/test_gate_and_ess_repair.py)."""
import os, sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "fit"))


def test_retired_repair_refuses_to_run():
    import repair_dd_calibrated as R
    with pytest.raises(SystemExit) as e:
        R.main()
    assert "retired" in str(e.value) and "repair_dd_ess.py" in str(e.value)
