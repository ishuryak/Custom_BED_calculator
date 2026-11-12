"""Build draws/model_methods.json: the priors, fixed settings, sampler settings and diagnostics of
the fits the app actually uses, proven against their retained fit identities (version 4.6.1,
review of 2026-11-12, item M3).

Nothing here is typed in from memory. Each prior is a statement that must occur verbatim in the
Stan file whose sha256 the fit identity recorded, and each value passed to Stan as data is read
from Stan data that is rebuilt from the data file and must reproduce the identity's recorded hash.
Any mismatch stops the build with the reason; the app's page then keeps the last valid record.

Usage:  python3 tools/build_model_methods.py [--root DIR] [--check]
  --check  rebuild in memory and exit 1 if draws/model_methods.json differs (no write).
"""
import argparse, hashlib, importlib, json, math, os, re, sys, types

SCHEMA = "model-methods-1"
TUMOR_FIT = "fit/output_refit_2026_10_25c_main"
TUMOR_DATA = "data/refit_trials_2026_10_25c_main.csv"
EARLY_UNIT = "docs/development/2026_11_10_early_toxicity/step5_review/overdispersed/U02"
# Versions of the environment that ran the early fit. The early manifest records the CmdStan
# version only; these two are what make its recorded identity reproduce (a sha256 match).
EARLY_ENV = dict(cmdstan="cmdstan-2.40.0", cmdstanpy="1.3.0", numpy="2.4.3")

# Priors and bounds of each fitted model. `stan`: statements that must occur in the model file;
# `data`: names of values passed as data, filled in from the reproduced Stan data.
DI_PRIORS = [
    dict(parameter="Tumor α/β", units="Gy", prior="lognormal(log {prior_log_ab_mean_exp:g}, {prior_log_ab_sd:g})",
         stan=["log_ab_raw ~ std_normal();",
               "exp(prior_log_ab_mean + prior_log_ab_sd * log_ab_raw)"],
         data=["prior_log_ab_mean", "prior_log_ab_sd"], status="assumption: the plausible range of the 2018 paper"),
    dict(parameter="Tumor α", units="per Gy", prior="lognormal(log 0.2, 1)",
         stan=["alpha ~ lognormal(log(0.2), 1);"], status="weakly informative"),
    dict(parameter="DI regrowth rate λ", units="per day", prior="lognormal(log 0.1, 1.5)",
         stan=["lambda ~ lognormal(log(0.1), 1.5);"], status="weakly informative"),
    dict(parameter="DI regrowth onset Tk", units="days", prior="normal(28, 10), limited to 0–60",
         stan=["Tk ~ normal(28, 10);", "real<lower=0, upper=60> Tk;"], status="centered on the classic onset near 4 weeks"),
    dict(parameter="Background growth g", units="per day", prior="uniform 0.0046–0.021",
         stan=["real<lower=0.0046, upper=0.021> g_free;"], absent=["g_free ~"],
         status="the plausible range of the 2018 paper"),
    dict(parameter="Baseline of each fitting unit", units="log(−log) scale", prior="normal(0, {prior_c_sd:g})",
         stan=["c ~ normal(0, prior_c_sd);"], data=["prior_c_sd"], status="weak, the same fixed prior for every unit"),
    dict(parameter="Spread of arm-level effects σ", units="log(−log) scale", prior="half-normal(0.5)",
         stan=["sigma_raw ~ normal(0, 0.5);", "real<lower=0> sigma_raw;"], status="weak; estimated common spread"),
]
DD_PRIORS = [
    dict(parameter="DD regrowth rate λ", units="per day", prior="lognormal(log {lam_med:g}, {lam_sd:g})",
         stan=["lambda_dd ~ lognormal(prior_log_lambda_mean, prior_log_lambda_sd);"],
         status="weakly informative; fitted given each DI draw"),
    dict(parameter="DD onset", units="days", prior="soft penalty favoring onset at or after day {onset_min:g} in every arm (scale {onset_scale:g} day), not a bound",
         stan=["target += -sum(exp((onset_min - onset) / onset_scale));"], status="assumption"),
    dict(parameter="Tumor α, α/β and g", units="", prior="taken from the DI draw (not refitted)",
         stan=["real<lower=0> alpha;        // from the DI draw", "real<lower=0> ab;           // from the DI draw"],
         status="conditional construction"),
    dict(parameter="DD onset constant C/α", units="Gy·day", prior="solved so that DD equals DI at 35 × 2 Gy in 47 days (not fitted)",
         stan=["real c_over_alpha = K0 * onset_ref / T0;"], status="calibration"),
    dict(parameter="Baseline of each fitting unit", units="log(−log) scale", prior="normal(0, {prior_c_sd:g})",
         stan=["c ~ normal(0, prior_c_sd);"], status="as in the DI fit"),
    dict(parameter="Spread of arm-level effects σ", units="log(−log) scale", prior="half-normal(0.5)",
         stan=["sigma_raw ~ normal(0, 0.5);", "real<lower=0> sigma_raw;"], status="as in the DI fit"),
]
EARLY_PRIORS = [
    dict(parameter="Mean trial baseline", units="log odds", prior="normal(0, 2.5)", stan=["m0 ~ normal(0, 2.5);"],
         status="weak"),
    dict(parameter="Spread of trial baselines", units="log odds", prior="half-normal(1)",
         stan=["tau0 ~ normal(0, 1);", "real<lower=0> tau0;"], status="weak; estimated"),
    dict(parameter="CTCAE version 4 grading offset", units="log odds", prior="normal(0, 1)", stan=["delta ~ normal(0, 1);"],
         status="weak"),
    dict(parameter="Slope per 10 Gy of early BED", units="log odds", prior="half-normal({Lstar_sd:g})",
         stan=["Lstar ~ normal(0, Lstar_sd);", "real<lower=0> Lstar;", "real kappa = Lstar / 10;"], data=["Lstar_sd"],
         status="weak; risk cannot fall with dose"),
    dict(parameter="Mucosal α/β", units="Gy", prior="lognormal(log {ab_mu_exp:g}, {ab_sd:g})",
         stan=["log_ab ~ normal(ab_mu, ab_sd);"], data=["ab_mu", "ab_sd"], status="assumption centered on Fowler's 10 Gy"),
    dict(parameter="Time-loss rate r", units="Gy of BED per day", prior="lognormal(log {r_mu_exp:g}, {r_sd:g})",
         stan=["log_r ~ normal(r_mu, r_sd);"], data=["r_mu", "r_sd"], status="centered near Fowler's 0.79"),
    dict(parameter="Planned-chemotherapy coefficient", units="log odds", prior="normal(0, {gamma_sd:g})",
         stan=["gamma ~ normal(0, gamma_sd);"], data=["gamma_sd"], status="weak"),
    dict(parameter="Spread of arm-level effects", units="log odds", prior="half-normal({od_sd:g})",
         stan=["sigma_arm ~ normal(0, od_sd);", "array[use_od] real<lower=0> sigma_arm;"], data=["od_sd"], status="weak; estimated"),
]


def sha_file(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def fail(msg):
    raise SystemExit(f"REFUSED: {msg}")


def check_statements(spec, stan_text, model):
    for p in spec:
        for s in p.get("stan", []):
            if s not in stan_text:
                fail(f"{model}: '{s}' ({p['parameter']}) is not in the Stan model the fit was run with")
        for s in p.get("absent", []):
            if s in stan_text:
                fail(f"{model}: '{s}' should not occur in the Stan model ({p['parameter']})")


def table(spec, values):
    return [dict(parameter=p["parameter"], units=p["units"], prior=p["prior"].format(**values), status=p["status"])
            for p in spec]


def tumor(root):
    sys.path.insert(0, os.path.join(root, "fit"))
    F = importlib.import_module("fit_refit")
    out = os.path.join(root, TUMOR_FIT)
    di = json.load(open(os.path.join(out, "di_primary", "identity.json")))
    dd = json.load(open(os.path.join(out, "dd_calibrated_onset7_n1000", "identity.json")))
    for name, ident in (("DI", di), ("calibrated DD", dd)):
        body = {k: v for k, v in ident.items() if k != "identity_sha256"}
        if hashlib.sha256(json.dumps(body, sort_keys=True).encode()).hexdigest() != ident["identity_sha256"]:
            fail(f"the {name} identity does not hash to its recorded identity_sha256 (edited after the fit)")
    stan_dir = os.path.join(root, "fit", "stan")
    if sha_file(os.path.join(stan_dir, "di.stan")) != di["stan_sha256"]["di.stan"]:
        fail("fit/stan/di.stan differs from the model the DI fit was run with")
    if sha_file(os.path.join(stan_dir, "dd_cal.stan")) != dd["stan_sha256"]:
        fail("fit/stan/dd_cal.stan differs from the model the calibrated DD fit was run with")
    if dd["parent_di_identity"] != di["identity_sha256"]:
        fail("the calibrated DD fit was not conditioned on this DI fit")
    data_path = os.path.join(root, TUMOR_DATA)
    if sha_file(data_path) != di["data_file_sha256"]:
        fail(f"{TUMOR_DATA} differs from the data the DI fit used")
    # Rebuild the DI Stan data with the historical anchor and evaluation regimens of every
    # 2026-10-25c identity (fit/repair_dd_ess.py LEGACY_*); it must reproduce the recorded hash.
    legacy_anchor = (35, 2.0, 47.0)
    legacy_regs = ((20, 2.75, 26), (15, 3.1, 19), (17, 3.05, 23), (25, 2.4, 33), (30, 65 / 30, 40),
                   (33, 2.0, 45), (30, 2.0, 40), (33, 2.0, 38), (35, 2.0, 40), (68, 1.2, 47))
    form = lambda n, d, T: [float(n * d), float(n * d * d), T]
    rows, trials, base = F.load(data_path)
    a = form(*legacy_anchor); regs = [form(*r) for r in legacy_regs]
    hist = dict(base, prior_c_sd=float(dd["prior_c_sd"]), R=len(regs), anchor_dsum=a[0], anchor_d2sum=a[1], anchor_T=a[2],
                eval_dsum=[r[0] for r in regs], eval_d2sum=[r[1] for r in regs], eval_T=[r[2] for r in regs])
    data = F.stan_data(hist, di["variant"])
    if hashlib.sha256(json.dumps(data, sort_keys=True, default=float).encode()).hexdigest() != di["stan_data_sha256"]:
        fail("the DI Stan data do not reproduce the identity's stan_data_sha256, so the data-passed priors are unproven")
    check_statements(DI_PRIORS, open(os.path.join(stan_dir, "di.stan"), encoding="utf-8").read(), "di.stan")
    check_statements(DD_PRIORS, open(os.path.join(stan_dir, "dd_cal.stan"), encoding="utf-8").read(), "dd_cal.stan")
    dval = dict(data, prior_log_ab_mean_exp=round(math.exp(data["prior_log_ab_mean"]), 6))
    lam = dd["prior_log_lambda"]
    ddval = dict(lam_med=round(math.exp(lam[0]), 6), lam_sd=lam[1], onset_min=dd["onset_min"],
                 onset_scale=dd["onset_scale"], prior_c_sd=dd["prior_c_sd"])
    man = json.load(open(os.path.join(root, "draws", "draws_manifest.json")))
    di_f, dd_f = man["files"]["tumor_di_primary.csv"], man["files"]["tumor_dd_primary.csv"]
    if di_f["fit_identity"] != di["identity_sha256"] or dd_f["fit_identity"] != dd["identity_sha256"]:
        fail("the exported tumor draws were not produced by these fits")
    ab_q = [round(math.exp(data["prior_log_ab_mean"] + z * data["prior_log_ab_sd"]), 1) for z in (-1.959964, 1.959964)]
    return dict(
        data_file=os.path.basename(TUMOR_DATA), data_sha256=di["data_file_sha256"],
        n_arms=base["N"], n_fitting_units=base["J"],
        likelihood="binomial: patients with tumor control in each arm, a rounded count from the reported control rate times the number analyzed",
        response="Poisson tumor control, TCP = exp(−exp(η))",
        di=dict(identity=di["identity_sha256"], stan_sha256=di["stan_sha256"]["di.stan"],
                stan_data_sha256=di["stan_data_sha256"], stan_data_reproduced=True,
                sampler=dict(chains=di["chains"], warmup=di["warmup"], samples=di["samples"],
                             adapt_delta=di["adapt_delta"], max_treedepth=di["max_treedepth"], seed=di["seed"]),
                software=di["software"], priors=table(DI_PRIORS, dval), ab_prior_95=ab_q,
                diagnostics=di_f["diagnostics"], draws_sha256=di_f["sha256"]),
        dd=dict(identity=dd["identity_sha256"], parent_di_identity=dd["parent_di_identity"], stan_sha256=dd["stan_sha256"],
                n_di_draws=dd["n_di"], onset_min=dd["onset_min"], onset_scale=dd["onset_scale"], sampler=dict(dd["inner"], max_attempts=dd["max_attempts"], seed=dd["seed"],
                                                    adapt_delta=0.99, max_treedepth=12),
                priors=table(DD_PRIORS, ddval), diagnostics=dd_f["diagnostics"], draws_sha256=dd_f["sha256"]),
        draws=dict(rows=di_f["n_draws"], distinct_di_draws=di_f["n_distinct_di_draws"],
                   dd_per_di_draw=di_f["n_draws"] // di_f["n_distinct_di_draws"],
                   pairing=man["pairing"], gate_policy=man["gate_policy"]["version"]),
    )


def early(root):
    if "cmdstanpy" not in sys.modules:
        try:
            import cmdstanpy  # noqa: F401
        except ImportError:  # the fit is not run here; only its data builder is used
            stub = types.ModuleType("cmdstanpy")
            stub.__version__, stub.CmdStanModel, stub.cmdstan_path = EARLY_ENV["cmdstanpy"], object, lambda: EARLY_ENV["cmdstan"]
            sys.modules["cmdstanpy"] = stub
    sys.path.insert(0, os.path.join(root, "tools", "early_tox"))
    S = importlib.import_module("step5_fit")
    unit = os.path.join(root, EARLY_UNIT)
    man = json.load(open(os.path.join(unit, "manifest.json")))
    diag = json.load(open(os.path.join(unit, "diagnostics.json")))
    if sha_file(os.path.join(unit, "diagnostics.json")) != man["outputs"]["diagnostics.json"]:
        fail("the early-fit diagnostics differ from those its manifest recorded")
    desc, settings = S.UNITS["U02"]
    settings = dict(settings, **diag["settings"])
    rows, fams = S.load_rows(False, False, False, None, False)
    data = S.stan_data(rows, fams, diag["Tk"], settings)
    payload = dict(data=S._normalize(data), settings=S._normalize(settings), sampler=S._normalize(S.SAMPLER),
                   stan_sha256=sha_file(S.STAN_FILE), cmdstan=EARLY_ENV["cmdstan"], cmdstanpy=EARLY_ENV["cmdstanpy"],
                   numpy=EARLY_ENV["numpy"], extra=S._normalize(dict(unit="U02", families=fams)))
    ident = hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()
    if ident != man["identity"]:
        fail("the early fit's identity does not reproduce from the current data, settings, priors and Stan model")
    check_statements(EARLY_PRIORS, open(S.STAN_FILE, encoding="utf-8").read(), "step5_model.stan")
    val = dict(data, ab_mu_exp=round(math.exp(data["ab_mu"]), 6), r_mu_exp=round(math.exp(data["r_mu"]), 6))
    pkg = json.load(open(os.path.join(root, "draws", "early_fit_manifest.json")))
    if pkg["source_identity"] != man["identity"]:
        fail("the exported early draws were not produced by this fit")
    s = diag["sampler"]
    return dict(
        identity=man["identity"], model_id=pkg["model_id"], stan_sha256=sha_file(S.STAN_FILE), identity_reproduced=True,
        n_families=diag["n_families"], n_arms=diag["n_arms"], n_evaluable=diag["totals"]["evaluable"],
        likelihood="binomial: patients with grade 3 or worse (or confluent) acute mucositis in each arm, logit link",
        fixed=dict(onset_elapsed_days=diag["Tk"], repair_halftime_h=diag["halftime_h"], hinge="exact"),
        sampler=dict(chains=s["chains"], warmup=s["iter_warmup"], samples=s["iter_sampling"],
                     adapt_delta=s["adapt_delta"], seed=s["seed"]),
        priors=table(EARLY_PRIORS, val),
        diagnostics=dict(divergences=diag["divergences"], max_rhat=diag["rhat"], min_ess_bulk=diag["ess_bulk"],
                         min_ess_tail=diag["ess_tail"], min_ebfmi=diag["ebfmi_min"],
                         treedepth_hits=diag["treedepth_saturated"], n_quantities=diag["n_required"],
                         gate_pass=diag["gate_pass"],
                         gate_note="tree-depth hits are recorded but are not part of this fit's pass rule"),
        draws=dict(rows=pkg["n_draws"], sha256=pkg["files"]["early_fit_draws.csv"]["sha256"]),
    )


def late(root):
    txt = open(os.path.join(root, "R", "uncertainty.R"), encoding="utf-8").read()
    m = re.search(r"ab_l = list\(dist = \"lognormal\", median = ([0-9.]+), sd_log = ([0-9.]+)", txt)
    b = re.search(r"ab_l = 3\)", open(os.path.join(root, "R", "bed_models.R"), encoding="utf-8").read())
    if not m or not b:
        fail("the late-toxicity α/β assumption could not be read from R/uncertainty.R and R/bed_models.R")
    med, sd = float(m.group(1)), float(m.group(2))
    q = [round(med * math.exp(z * sd), 2) for z in (-1.959964, 1.959964)]
    return dict(status="not fitted", point_ab=3.0, uncertainty=dict(dist="lognormal", median=med, sd_log=sd, central_95=q),
                note="an analyst-assumed range, not an estimate from trial data; complete repair, no time factor")


def build(root):
    t, e, l = tumor(root), early(root), late(root)
    sources = [TUMOR_DATA, "fit/stan/di.stan", "fit/stan/dd_cal.stan", "fit/fit_refit.py",
               f"{TUMOR_FIT}/di_primary/identity.json", f"{TUMOR_FIT}/dd_calibrated_onset7_n1000/identity.json",
               "draws/draws_manifest.json", "tools/early_tox/step5_model.stan", "tools/early_tox/step5_fit.py",
               f"{EARLY_UNIT}/manifest.json", f"{EARLY_UNIT}/diagnostics.json", "draws/early_fit_manifest.json",
               "R/uncertainty.R", "tools/build_model_methods.py"]
    return dict(schema=SCHEMA, created_by="tools/build_model_methods.py", tumor=t, early=e, late=l,
                toxicity_pairing="tumor and toxicity draws are paired independently (an assumption, not a fitted joint distribution)",
                sources={p: sha_file(os.path.join(root, p)) for p in sources})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    text = json.dumps(build(a.root), indent=1, ensure_ascii=False, sort_keys=True) + "\n"
    path = os.path.join(a.root, "draws", "model_methods.json")
    if a.check:
        cur = open(path, encoding="utf-8").read() if os.path.exists(path) else ""
        if cur != text:
            sys.exit("draws/model_methods.json is out of date: run python3 tools/build_model_methods.py")
        print("draws/model_methods.json is current")
        return
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    print("wrote", path)


if __name__ == "__main__":
    main()
