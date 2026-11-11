"""Early-toxicity step 5.3: historical classifier on complete-course Pardo-Montero labels."""
import importlib.util
import os

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_spec = importlib.util.spec_from_file_location("h", os.path.join(ROOT, "tools", "early_tox", "step5_historical.py"))
h = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(h)


def test_selections_and_partial_courses_never_enter(monkeypatch):
    for sel, n, e in (("primary", 103, 13), ("sensitivity", 107, 15)):
        recs = h.load_records(sel)
        assert (len(recs), sum(x["y"] for x in recs)) == (n, e)
    real = h.read

    def with_partial(name):
        rows = real(name)
        if name == "labels.csv":
            rows = [dict(r) for r in rows]
            rows[0]["course_scope"] = "partial"
        return rows
    monkeypatch.setattr(h, "read", with_partial)
    with pytest.raises(ValueError, match="partial course"):
        h.load_records("primary")


def test_fowler_exposure_reproduces_the_plan_fixture():
    recs = h.load_records("primary")
    i = next(k for k, x in enumerate(recs) if x["schedule_id"] == "S07394dd030")     # 70 Gy / 35
    assert h.exposure(recs, "fowler")[i] == pytest.approx(53.10544, abs=1e-5)


def test_folds_hold_out_each_source_once():
    recs = h.load_records("primary")
    held = np.zeros(len(recs), int)
    for s, tr, te in h.folds(recs):
        assert not (tr & te).any() and (tr | te).all()
        assert {recs[i]["source"] for i in np.flatnonzero(te)} == {s}
        held += te
    assert (held == 1).all()


@pytest.mark.parametrize("model", ["fowler", "boundary"])
def test_held_out_labels_do_not_leak_into_their_own_prediction(model):
    recs = h.load_records("primary")
    src = next(x["source"] for x in recs if x["y"] == 1)
    ls0, _ = h.loso(recs, model)
    flipped = [dict(x, y=1 - x["y"]) if x["source"] == src else x for x in recs]
    ls1, _ = h.loso(flipped, model)
    idx = [i for i, x in enumerate(recs) if x["source"] == src]
    # p is unchanged when the held-out labels flip: exp(log p) + exp(log(1 - p)) = 1
    assert np.allclose(np.exp(ls0[idx]) + np.exp(ls1[idx]), 1.0, atol=1e-6)


def test_tk_is_chosen_without_the_held_out_source(monkeypatch):
    recs = h.load_records("primary")
    seen = []
    real = h.choose_tk
    monkeypatch.setattr(h, "choose_tk", lambda rs: seen.append({x["source"] for x in rs}) or real(rs))
    h.loso(recs, "boundary")
    srcs = [s for s, _, _ in h.folds(recs)]
    assert len(seen) == len(srcs)
    assert all(s not in used for s, used in zip(srcs, seen))


def test_matched_repair_variants_differ_only_in_repair():
    recs = h.load_records("primary")
    s5pt = (12.0, 1.5)
    bid = next(i for i, x in enumerate(recs) if len(set(int(t // 24) for t in x["times"])) < len(x["times"]))
    one = next(i for i, x in enumerate(recs) if len(set(int(t // 24) for t in x["times"])) == len(x["times"]))
    f0, f3 = h.exposure(recs, "fowler"), h.exposure(recs, "fowler_h3.2")
    c3, c0 = h.exposure(recs, "step5", step5=s5pt), h.exposure(recs, "step5_complete", step5=s5pt)
    assert f3[bid] > f0[bid] and c3[bid] > c0[bid]          # incomplete repair raises BID exposure
    assert f3[one] == pytest.approx(f0[one], abs=0.2)       # once daily: 24 h apart, about 0.1 Gy10
