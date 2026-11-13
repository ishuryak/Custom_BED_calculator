# November 12 software review evidence

All original project reads were read-only. The canonical snapshot is `source/`; its revision is `79e9c812f91482fa57ffa1dddabe89dc329a6303`. Two untracked current plans are separately captured in `latest_plans/`. No fits, deployment, dependency installation, or original file modifications were performed.

Run the following in this evidence directory in PowerShell. The scripts accept source/cache paths so a later reviewer can substitute a corresponding checkout. `artifact_audit.py` writes JSONs and `rebuilt/` only in the current directory; mutation probes create temporary directories under the current directory. They do not modify the source/cache paths supplied.

```powershell
$python = 'C:/Users/ishur/AppData/Local/Programs/Python/Python312/python.exe'
$source = (Resolve-Path ./source).Path
$cache = 'C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v45/docs/development/2026_11_10_early_toxicity/step6_review/cache'
& $python -B ./source/tools/build_model_methods.py --root $source --check
& $python -B ./guard_probes.py $source
& $python -B ./cache_probe.py $source $cache
& $python -B ./methods_probe.py $source
& $python -B ./artifact_audit.py $source $cache
& $python -B ./retained_checks.py $source $cache
```

`guard_probes.py` and `cache_probe.py` explicitly stub CmdStan imports and prohibit sampling. The cache probe uses the original recorded NumPy version string for deterministic identity reconstruction only; actual array calculations run with installed NumPy 1.26.4. `methods_probe.py` runs the actual metadata builder in fresh subprocesses against temporary copies. No stub alters its generated scientific metadata.

Native Python 3.12.10 with NumPy 1.26.4, SciPy 1.15.3 and pandas 2.2.3 is available; pytest and cmdstanpy are not. The attempted suite command failed before collection (`No module named pytest`), so this evidence does not claim a full Python suite pass. The checked raw-chain recheck CSV is existing maintainer evidence. Its 114 identities match current artifacts, but the Linux raw-chain archive was not re-read or recomputed by this review.

The report `software_review.md` distinguishes passing controls, deliberately exposed residuals and future acceptance tests. A zero exit status from a reproduction script means the script ran; it does not mean every malformed input was rejected. In particular, the exporter and methods probes intentionally print accepted cases that constitute findings.

`evidence_inventory.csv` lists file SHA-256 digests for this contribution, excluding the full `source/` tree and `source.zip`. `source_manifest.json` independently records every tracked snapshot file's Git blob/SHA-256 and original-worktree equality. Do not package copied runtime libraries or the 52 NPZ caches into the review evidence bundle.
