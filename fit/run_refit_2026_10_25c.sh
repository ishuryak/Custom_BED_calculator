#!/usr/bin/env bash
# Corrected refit of 2026-10-25 (plan PLAN_codex_review_2026_10_25.md, stages B-D).
# Builds every data set, then fits main (all variants, leave-one-study-out, calibrated DD on 1000
# DI draws for the primary and the fixed-alpha/beta variant) and the sensitivity sets (primary
# variant, calibrated DD on 200 DI draws). Each fit resumes if killed: finished fits are reused
# only when their identity matches. Run from the repository root:
#   bash fit/run_refit_2026_10_25c.sh 2>&1 | tee -a fit/run_refit_2026_10_25c.log
set -euo pipefail
cd "$(dirname "$0")/.."
TAG=2026_10_25c
echo "start $(date -Is) HEAD $(git rev-parse HEAD) (code may be uncommitted; each fit's identity.json binds its file hashes)"
sha256sum fit/fit_refit.py fit/fit_dd_calibrated.py fit/build_refit_dataset.py fit/stan/*.stan fit/stan/*.stanfunctions

declare -A FLAGS=([main]="" [effn_half]="--n-scale 0.5" [nominal_clock]="--nominal-clock" [km_only]="--km-only"
                  [rtog9003_beitler]="--rtog9003-beitler" [artscan]="--artscan" [new_trials]="--add-2026-10-24"
                  [chemo]="--chemo")
for s in "${!FLAGS[@]}"; do
  python3 -I fit/build_refit_dataset.py "data/refit_trials_${TAG}_${s}.csv" ${FLAGS[$s]}
done
sha256sum data/refit_trials_${TAG}_*.csv

fit_set () {   # name data_csv ndi variants loto prior_c_sd
  local name=$1 csv=$2 ndi=$3 variants=$4 loto=$5 csd=$6
  local out="fit/output_refit_${TAG}_${name}"
  echo "=== $name $(date -Is) data $csv out $out"
  BED_FIT_DATA="$csv" BED_FIT_OUT="$out" BED_FIT_VARIANTS="$variants" BED_FIT_PRIOR_C_SD="$csd" \
    python3 fit/fit_refit.py $loto
  for v in ${7:-primary}; do
    BED_FIT_DATA="$csv" BED_FIT_OUT="$out" BED_FIT_PRIOR_C_SD="$csd" DDCAL_NDI="$ndi" DDCAL_DI_VARIANT="$v" \
      python3 -I fit/fit_dd_calibrated.py
  done
}

fit_set main data/refit_trials_${TAG}_main.csv 1000 primary,binomial,paper_restrictions,wide_ab --loto 3 "primary paper_restrictions"
for s in effn_half nominal_clock km_only rtog9003_beitler artscan new_trials chemo; do
  fit_set "$s" "data/refit_trials_${TAG}_${s}.csv" 200 primary "" 3
done
fit_set prior_c10 data/refit_trials_${TAG}_main.csv 200 primary "" 10
echo "ALLDONE $(date -Is)"
