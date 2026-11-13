from pathlib import Path
import hashlib,json
here=Path(__file__).resolve().parent
project=Path('C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v45')
sha=lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
rows=json.loads((here/'software/source_manifest.json').read_text())
print('source_manifest_keys',list(rows[0]))
changed=[]
for row in rows:
    # Snapshot bytes establish the comparison regardless of the manifest's naming convention.
    p=Path(row['path'])
    if (project/p).read_bytes() != (here/'software/source'/p).read_bytes(): changed.append(p.as_posix())
plans={}
for name in ['2026_11_12_frequentist_comparison_plan.md','2026_11_12_frequentist_comparison_plan_audit.md']:
    plans[name]={'sha256':sha(project/'docs/development'/name),'matches_capture':(project/'docs/development'/name).read_bytes()==(here/'software/latest_plans'/name).read_bytes()}
orig=Path('C:/Users/ishur/AppData/Local/R/cache/R/HNSCC_BED_calculator/library-2026-10-15-R4.5.1/DT')
local=Path('C:/Users/ishur/Documents/Codex/2026-10-18/rev-2/work/review_2026_11_09/software/local_render_library/DT')
original_files={p.relative_to(orig).as_posix():sha(p) for p in orig.rglob('*') if p.is_file()}
local_files={p.relative_to(local).as_posix():sha(p) for p in local.rglob('*') if p.is_file()}
ans={'tracked_files':len(rows),'changed_original_vs_snapshot':changed,'plans':plans,
     'DT_original_files':len(original_files),'DT_local_files':len(local_files),'DT_byte_identical':original_files==local_files}
(here/'final_checks.json').write_text(json.dumps(ans,indent=2))
print(json.dumps(ans,indent=2))
assert not changed and all(x['matches_capture'] for x in plans.values()) and original_files==local_files
