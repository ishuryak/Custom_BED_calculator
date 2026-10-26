#!/usr/bin/env bash
# Re-fit the calibrated-DD fits that failed the gate (fit/repair_dd_calibrated.py), main sets first.
set -euo pipefail
cd "$(dirname "$0")/.."
TAG=2026_10_25c
echo "revision $(git rev-parse HEAD) $(date -Is)"
sha256sum fit/repair_dd_calibrated.py fit/fit_dd_calibrated.py fit/fit_refit.py
rep () {   # name data_csv ndi variant prior_c_sd
  echo "=== repair $1 $4 $(date -Is)"
  BED_FIT_DATA="$2" BED_FIT_OUT="fit/output_refit_${TAG}_$1" BED_FIT_PRIOR_C_SD="$5" DDCAL_NDI="$3" DDCAL_DI_VARIANT="$4" \
    python3 -I fit/repair_dd_calibrated.py
}
rep main data/refit_trials_${TAG}_main.csv 1000 primary 3
rep main data/refit_trials_${TAG}_main.csv 1000 paper_restrictions 3
for s in effn_half km_only rtog9003_beitler artscan chemo; do
  rep "$s" "data/refit_trials_${TAG}_${s}.csv" 200 primary 3
done
rep prior_c10 data/refit_trials_${TAG}_main.csv 200 primary 10
echo "ALLDONE $(date -Is)"
