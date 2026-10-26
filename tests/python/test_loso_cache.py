"""fit/fit_refit.py leave-one-study-out cache (review R02, 2026-10-26): saved scores are reused only
for the same identity; a cache from other data or another runner is refused before any fit."""
import importlib, json, os, sys
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "fit"))
DATA = os.path.join(ROOT, "data", "refit_trials_2026_10_25c_main.csv")


def load(monkeypatch, tmp_path):
    monkeypatch.setenv("BED_FIT_DATA", DATA); monkeypatch.setenv("BED_FIT_OUT", str(tmp_path))
    import fit_refit
    F = importlib.reload(fit_refit)
    rows, trials, base = F.load()
    return F, rows, trials, base


def test_stale_holdout_cache_is_refused(monkeypatch, tmp_path):
    F, rows, trials, base = load(monkeypatch, tmp_path)
    # a cache written before identities existed, or by another data set: scores only
    json.dump({"CAIR": dict(elpd=-1.0)}, open(tmp_path / "loso_di.json", "w"))
    with pytest.raises(SystemExit) as e:
        F.loto("di", None, base, rows, trials)
    assert "REFUSED" in str(e.value)


def test_matching_holdout_cache_is_reused_without_fitting(monkeypatch, tmp_path):
    F, rows, trials, base = load(monkeypatch, tmp_path)
    studies = sorted({r["first_author"] for r in rows}, key=[r["first_author"] for r in rows].index)
    ident = F.identity("di", F.stan_data(base, "primary"), F.SEED, 4, min(800, F.WARMUP), min(500, F.SAMPLES),
                       extra=dict(analysis="leave one study out", studies=studies, units=list(trials),
                                  prior_c_sd=F.PRIOR_C_SD, score="conditional plug-in, 20 normal draws per posterior draw"))
    cache = {"_identity": ident, **{s: dict(elpd=-1.0) for s in studies}}
    json.dump(cache, open(tmp_path / "loso_di.json", "w"))
    res = F.loto("di", None, base, rows, trials)        # m=None: any attempt to fit would fail
    assert set(res) - {"_identity"} == set(studies)
    # the same cache with one patient count changed in the data is refused
    base2 = dict(base, npat=[n + (i == 0) for i, n in enumerate(base["npat"])])
    with pytest.raises(SystemExit):
        F.loto("di", None, base2, rows, trials)
