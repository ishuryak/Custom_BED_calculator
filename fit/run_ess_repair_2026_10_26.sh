#!/usr/bin/env bash
# Review R01 (2026-10-26): refit every calibrated-DD fit that fails the full gate (fit/gate.py,
# including ESS >= 400) with fit/repair_dd_ess.py. One folder per process, five in parallel
# (4 chains each, 20 of 22 cores). Each fit is saved as it finishes, so a killed run resumes.
set -euo pipefail
cd "$(dirname "$0")/.."
TAG=2026_10_25c
LOGDIR=fit/logs_ess_repair_2026_10_26; mkdir -p "$LOGDIR"
echo "revision $(git rev-parse HEAD) $(date -Is)"
[ -z "$(git status --porcelain --untracked-files=no -- fit)" ] || { echo "REFUSED: fit/ has uncommitted changes"; exit 1; }
sha256sum fit/repair_dd_ess.py fit/fit_dd_calibrated.py fit/fit_refit.py fit/gate.py
cat > "$LOGDIR/jobs.txt" <<JOBS
main data/refit_trials_${TAG}_main.csv 1000 primary 3
main data/refit_trials_${TAG}_main.csv 1000 paper_restrictions 3
effn_half data/refit_trials_${TAG}_effn_half.csv 200 primary 3
km_only data/refit_trials_${TAG}_km_only.csv 200 primary 3
rtog9003_beitler data/refit_trials_${TAG}_rtog9003_beitler.csv 200 primary 3
artscan data/refit_trials_${TAG}_artscan.csv 200 primary 3
chemo data/refit_trials_${TAG}_chemo.csv 200 primary 3
new_trials data/refit_trials_${TAG}_new_trials.csv 200 primary 3
nominal_clock data/refit_trials_${TAG}_nominal_clock.csv 200 primary 3
prior_c10 data/refit_trials_${TAG}_main.csv 200 primary 10
JOBS
run1 () {   # name data_csv ndi variant prior_c_sd
  BED_FIT_DATA="$2" BED_FIT_OUT="fit/output_refit_${TAG}_$1" BED_FIT_PRIOR_C_SD="$5" DDCAL_NDI="$3" DDCAL_DI_VARIANT="$4" \
    python3 -I fit/repair_dd_ess.py > "$LOGDIR/$1_$4.log" 2>&1 && echo "DONE $1 $4 $(date -Is)" || echo "FAILED $1 $4 $(date -Is)"
}
export -f run1; export TAG LOGDIR
xargs -P 5 -L 1 bash -c 'run1 "$@"' _ < "$LOGDIR/jobs.txt"
echo "ALLDONE $(date -Is)"
