"""draws/model_methods.json (version 4.6.1, review 2026-11-12 item M3): the methods page's priors,
fixed settings and diagnostics must be proven against the retained fits, and the build must refuse
when any of them changes without a refit. Each mutation test edits a temporary copy of the tree."""
import json, os, shutil, subprocess, sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TOOL = os.path.join(ROOT, "tools", "build_model_methods.py")
UNIT = "docs/development/2026_11_10_early_toxicity/step5_review/overdispersed/U02"
FIT = "fit/output_refit_2026_10_25c_main"


def run(root, *args):
    # A fresh interpreter per build, so modules are imported from `root`, not cached from ROOT.
    return subprocess.run([sys.executable, TOOL, "--root", root, *args], capture_output=True, text=True)


@pytest.fixture
def tree(tmp_path):
    root = str(tmp_path / "app")
    for d in ("data", "tools/early_tox", "fit/stan", "draws", "R"):
        shutil.copytree(os.path.join(ROOT, d), os.path.join(root, d),
                        ignore=shutil.ignore_patterns("__pycache__", "*.csv.gz"))
    for f in ("fit/fit_refit.py", "fit/gate.py", "tools/build_model_methods.py",
              f"{FIT}/di_primary/identity.json", f"{FIT}/dd_calibrated_onset7_n1000/identity.json",
              f"{UNIT}/manifest.json", f"{UNIT}/diagnostics.json"):
        os.makedirs(os.path.dirname(os.path.join(root, f)), exist_ok=True)
        shutil.copy(os.path.join(ROOT, f), os.path.join(root, f))
    return root


def edit(root, rel, old, new):
    p = os.path.join(root, rel)
    txt = open(p, encoding="utf-8").read()
    assert old in txt, f"{old!r} not in {rel}"
    open(p, "w", encoding="utf-8", newline="\n").write(txt.replace(old, new, 1))


def test_committed_record_is_current():
    r = run(ROOT, "--check")
    assert r.returncode == 0, r.stdout + r.stderr


def test_unchanged_copy_builds(tree):
    r = run(tree)
    assert r.returncode == 0, r.stderr
    m = json.load(open(os.path.join(tree, "draws", "model_methods.json"), encoding="utf-8"))
    assert m["tumor"]["di"]["stan_data_reproduced"] and m["early"]["identity_reproduced"]
    assert m["late"]["status"] == "not fitted"


def test_record_states_the_active_settings():
    m = json.load(open(os.path.join(ROOT, "draws", "model_methods.json"), encoding="utf-8"))
    di = {p["parameter"]: p["prior"] for p in m["tumor"]["di"]["priors"]}
    assert di["Tumor α/β"] == "lognormal(log 10, 0.354)" and di["DI regrowth onset Tk"] == "normal(28, 10), limited to 0–60"
    assert m["tumor"]["di"]["ab_prior_95"] == [5.0, 20.0]
    dd = {p["parameter"]: p["prior"] for p in m["tumor"]["dd"]["priors"]}
    assert dd["DD regrowth rate λ"] == "lognormal(log 0.5, 1.5)"
    assert dd["DD onset"].startswith("soft penalty favoring onset at or after day 7") and "not a bound" in dd["DD onset"]
    early = {p["parameter"]: p["prior"] for p in m["early"]["priors"]}
    assert early["Mucosal α/β"] == "lognormal(log 10, 0.354)" and early["Time-loss rate r"] == "lognormal(log 0.8, 0.5)"
    assert m["early"]["fixed"] == {"onset_elapsed_days": 7.0, "repair_halftime_h": 3.2, "hinge": "exact"}
    assert m["tumor"]["draws"]["distinct_di_draws"] == 1000 and m["tumor"]["draws"]["dd_per_di_draw"] == 4
    assert (m["tumor"]["n_arms"], m["tumor"]["n_fitting_units"]) == (45, 21)
    assert m["late"]["uncertainty"]["central_95"] == [2.0, 4.5]


@pytest.mark.parametrize("rel, old, new, why", [
    ("fit/stan/di.stan", "real<lower=0, upper=60> Tk;", "real<lower=0, upper=70> Tk;", "di.stan differs"),
    ("fit/stan/di.stan", "Tk ~ normal(28, 10);", "Tk ~ normal(28, 12);", "di.stan differs"),
    ("fit/stan/dd_cal.stan", "exp((onset_min - onset) / onset_scale)", "exp((onset_min - onset) / (2 * onset_scale))",
     "dd_cal.stan differs"),
    (f"{FIT}/dd_calibrated_onset7_n1000/identity.json", '"onset_scale": 1.0', '"onset_scale": 2.0',
     "does not hash to its recorded identity_sha256"),
    ("fit/fit_refit.py", 'dict(prior_log_ab_sd=0.354, fix_ab=0, fix_g=0, overdispersion=1),\n    "binomial"',
     'dict(prior_log_ab_sd=0.5, fix_ab=0, fix_g=0, overdispersion=1),\n    "binomial"', "stan_data_sha256"),
    ("tools/early_tox/step5_fit.py", "ab_sd=0.354", "ab_sd=0.5", "early fit's identity does not reproduce"),
    ("tools/early_tox/step5_fit.py", "gamma_sd=1.0", "gamma_sd=2.0", "early fit's identity does not reproduce"),
    ("tools/early_tox/step5_model.stan", "sigma_arm ~ normal(0, od_sd);", "sigma_arm ~ normal(0, 2 * od_sd);",
     "early fit's identity does not reproduce"),
    (f"{UNIT}/diagnostics.json", '"halftime_h": 3.2', '"halftime_h": 4.0', "diagnostics differ"),
])
def test_a_changed_prior_bound_or_setting_is_refused(tree, rel, old, new, why):
    edit(tree, rel, old, new)
    r = run(tree)
    assert r.returncode != 0 and "REFUSED" in (r.stdout + r.stderr) and why in (r.stdout + r.stderr), r.stdout + r.stderr


def test_a_changed_late_assumption_makes_the_record_stale(tree):
    assert run(tree).returncode == 0
    edit(tree, "R/uncertainty.R", "median = 3, sd_log = 0.207", "median = 3, sd_log = 0.3")
    r = run(tree, "--check")
    assert r.returncode != 0 and "out of date" in (r.stdout + r.stderr)


def test_a_described_prior_missing_from_the_model_is_refused(tree):
    # The page text and the model disagree: the generator's statement list names a prior the
    # Stan file does not contain. Simulated by editing the generator copy's expected statement.
    edit(tree, "tools/build_model_methods.py", 'stan=["alpha ~ lognormal(log(0.2), 1);"]',
         'stan=["alpha ~ lognormal(log(0.3), 1);"]')
    r = subprocess.run([sys.executable, os.path.join(tree, "tools", "build_model_methods.py"), "--root", tree],
                       capture_output=True, text=True)
    assert r.returncode != 0 and "is not in the Stan model" in (r.stdout + r.stderr)
