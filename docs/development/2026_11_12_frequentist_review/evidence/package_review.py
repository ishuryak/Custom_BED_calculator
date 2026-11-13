from pathlib import Path
import csv,hashlib,io,json,re,zipfile

HERE=Path(__file__).resolve().parent
WORK=HERE.parent.parent
REPORT=WORK/'HNSCC_BED_frequentist_plan_and_Bayesian_refit_review_2026-11-12.md'
ZIP=WORK/'HNSCC_BED_frequentist_Bayesian_review_evidence_2026-11-12.zip'
body=REPORT.read_text(encoding='utf-8')
assert all(f'## F{i:02}.' in body for i in range(1,7))
assert all(f'## B{i:02}.' in body for i in range(1,5))
assert not re.search(r'RESULT_TO_INSERT|\bTODO\b|\bTBD\b|after the active run',body)
assert '\ufffd' not in body and body.count('```')%2==0
items={REPORT.name:REPORT.read_bytes()}
for p in sorted(HERE.iterdir()):
    if p.is_file() and p.suffix in ('.R','.py','.json','.log','.csv','.md'):
        items[p.name]=p.read_bytes()
for sub in ('math','clinical'):
    for p in sorted((HERE/sub).iterdir()):
        if p.is_file() and p.suffix in ('.py','.json','.log','.csv','.md','.txt','.png'):
            items[f'{sub}/{p.name}']=p.read_bytes()
with (HERE/'software/evidence_inventory.csv').open(encoding='utf-8-sig',newline='') as f:
    for row in csv.DictReader(f):
        data=(HERE/'software'/row['path']).read_bytes()
        assert hashlib.sha256(data).hexdigest()==row['sha256'],row['path']
        items['software/'+row['path'].replace('\\','/')]=data
items['software/evidence_inventory.csv']=(HERE/'software/evidence_inventory.csv').read_bytes()
items['README.md']='''# Review evidence, 12 November 2026

The governing consolidated assessment is the root report
HNSCC_BED_frequentist_plan_and_Bayesian_refit_review_2026-11-12.md.
It contains six frequentist-plan findings and four current Bayesian
software/source/documentation findings, with one ordered execution roadmap.
Domain reports provide additional derivation and reproduction details.

Reviewed source: Custom_BED_calculator_v45, version 4.6.1,
HEAD 79e9c812f91482fa57ffa1dddabe89dc329a6303.
The current frequentist plan and audit were physical uncommitted additions;
their exact copies and digests are included under software/latest_plans
and in the mathematical contribution. They were not implemented or fitted
during this review. All 1,542 tracked files and both plans were checked again
at the end and remained unchanged.

## Evidence map

- plan_numeric_checks.py/.json: independently verified tumor binomial constant,
  historical likelihood and all ten selected-posterior contrast medians.
- math/plan_probes.py/.json/.log: quadrature counterexamples and positive
  adaptive controls; finite-cluster illustration; rounded early contrast
  tolerances; exact unchanged current primary Stan input and artifact checks.
- runtime_checks.R/.json/.log: current calendar and replay corrections,
  legitimate package positive, malformed package negatives that still pass.
- run_r_review.R, r_full.log, r_full_results.csv: complete testthat run,
  6,781 expectations, zero failures, one environment-level error. All 285
  expectations in the v4.6/v4.6.1 test files passed. Not a release-pass marker.
- readability_probe.R/.log: direct dataset read/hash succeeds although
  file.access reports unreadable; four recorded gate flags also checked.
- final_checks.py/.json: final tracked-file/plan equality and all 368 copied
  DT dependency files byte-identical to the pinned release installation.
- software/: canonical identities, cache/report/output inventories,
  exporter and metadata negative probes, builder reconstruction and limits.
- clinical/: clinical report and selected pages 391/394 of the user's local
  HYPNO PDF used for source verification. The full publication is not bundled.

The live main/methods check observed version 4.6.1, deployed revision 7ec2faf,
69 main points, correct early-section routing and the rendered tumor formula.
It was a targeted visual/accessible-UI check, not a new 73-item walkthrough.

## Reproduction and limits

Python numerical probes used 3.12.10, NumPy 1.26.4, SciPy 1.15.3, pandas 2.2.3.
pytest/CmdStanPy were unavailable in that interpreter. Stubbed imports in
specific probes explicitly prohibit sampling; methods-builder --check uses
its existing data-only fallback. No new MCMC/bootstrap was run.

R used 4.5.1, Windows UTF-8 and the project's release library, with a verified
byte-identical existing workspace DT copy to avoid resource-directory access
restrictions. Full-suite error remains recorded, not waived or marked passed.
readability_probe independently checks the blocked digest/gate assertions.

See software/README.md and math/math_review.md for probe commands. Supply
an exact source checkout and, for cache checks, the original retained caches.
The archive intentionally excludes copied application source, installed
libraries, full publications, posterior NPZ files and external raw chains.
Root R scripts and some root Python scripts record absolute review paths;
edit their source/output/library path literals for an isolated workspace
before rerunning elsewhere. Retention comparisons use the earlier reviewed
snapshot; restore that exact revision if reproducing those checks elsewhere.

Proposed acceptance tests and future runner invocations are not claimed to
have run. A zero-exit probe can demonstrate an invalid input was accepted;
read its labeled result instead of interpreting execution success as a gate
pass. Deliberately invalid fixtures were never installed in the project.

EVIDENCE_MANIFEST.csv records the SHA-256 and byte count of every other
archive member. No production edits, Git changes or deployment were made.
'''.encode('utf-8')
buf=io.StringIO(newline='');w=csv.writer(buf);w.writerow(['path','bytes','sha256'])
for name,data in sorted(items.items()):w.writerow([name,len(data),hashlib.sha256(data).hexdigest()])
items['EVIDENCE_MANIFEST.csv']=buf.getvalue().encode('utf-8')
with zipfile.ZipFile(ZIP,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
    for name,data in sorted(items.items()):z.writestr(name,data)
with zipfile.ZipFile(ZIP) as z:
    assert z.testzip() is None
    for row in csv.DictReader(io.StringIO(z.read('EVIDENCE_MANIFEST.csv').decode())):
        data=z.read(row['path'])
        assert len(data)==int(row['bytes']) and hashlib.sha256(data).hexdigest()==row['sha256']
print(json.dumps({'report':str(REPORT),'words':len(body.split()),'report_sha256':hashlib.sha256(REPORT.read_bytes()).hexdigest(),
                  'zip':str(ZIP),'members':len(items),'zip_bytes':ZIP.stat().st_size,'zip_sha256':hashlib.sha256(ZIP.read_bytes()).hexdigest(),
                  'qa':'10 findings; no draft placeholders; balanced code fences; all archive members verified'},indent=2))
