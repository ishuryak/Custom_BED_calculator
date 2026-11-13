#!/bin/bash
# Full run of plan revision 3 (after the T1-T8 gate and the B = 20 pilot). Every step resumes:
# completed, verified records of the same identity are kept; nothing completed is rerun.
# Run from the repository root:  bash fit/frequentist/run_full.sh   (log: fit/output_frequentist/run_full.log)
set -euo pipefail
FM="python3 -W ignore fit/frequentist/fit_ml.py"
echo "revision $(git rev-parse HEAD) $(date -Is)"
$FM validate
$FM profiles --model tumor --support matched
$FM profiles --model tumor --support expanded
$FM profiles --model early
$FM sigma_sensitivity --model tumor --support matched
$FM bootstrap --model tumor --support matched --B 1000
$FM bootstrap --model early --B 500
python3 -W ignore fit/frequentist/r_check.py
echo "done $(date -Is)"
