"""Read-only retained output verification; no sampler import. python artifact_audit.py SOURCE CACHE"""
import argparse,csv,hashlib,importlib.util,json,pathlib
import numpy as np
p=argparse.ArgumentParser();p.add_argument('repo',type=pathlib.Path);p.add_argument('cache',type=pathlib.Path);a=p.parse_args();root=a.repo.resolve()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
base=root/'docs/development/2026_11_10_early_toxicity/step5_review';records=[]
expected={'summary.csv','ppc_arms.csv','diagnostics.json','draws.csv.gz'}
for path in sorted(base.glob('*/*/manifest.json')):
    m=json.loads(path.read_text());d=json.loads((path.parent/'diagnostics.json').read_text())
    records.append(dict(unit=str(path.relative_to(base).parent),identity=m['identity'],gate=m['gate_pass'],declared_keys=sorted(m['outputs']),required_keys_present=expected<=set(m['outputs']),hashes_match=all((path.parent/f).exists() and sha(path.parent/f)==h for f,h in m['outputs'].items()),diagnostic_identity_match=m['identity']==d['identity'],diagnostic_gate_match=m['gate_pass']==d['gate_pass'],chain_hash_count=len(m.get('chain_csv_sha256',{})),archive=m['chain_archive'],ebfmi_min=d['ebfmi_min']))
pathlib.Path('step5_manifest_audit.json').write_text(json.dumps(records,indent=1))
print('Step5 manifests',len(records),'output hashes all match',all(r['hashes_match'] for r in records),'required files all present',all(r['required_keys_present'] for r in records),'gates passing',sum(r['gate'] for r in records),'diagnostic identities/gates match',all(r['diagnostic_identity_match'] and r['diagnostic_gate_match'] for r in records),'chain hash counts',sorted(set(r['chain_hash_count'] for r in records)))
out=[]
for p in sorted(a.cache.glob('*.npz')):
    with np.load(p,allow_pickle=False) as z:
        g=json.loads(str(z['gate_json']));archive=json.loads(str(z['chain_archive']))
        out.append(dict(file=p.name,sha256=sha(p),identity=str(z['identity']),n_draws=len(z['kappa']),gate=g['gate_pass'],scope=g['scope'],ebfmi_min=g['ebfmi_min'],chain_hash_count=len(archive),all_numeric_finite=all(np.isfinite(z[k]).all() for k in z.files if z[k].dtype.kind in 'fiu')))
pathlib.Path('step6_cache_audit.json').write_text(json.dumps(out,indent=1))
print('Step6 caches',len(out),'n_draws',sorted(set(r['n_draws'] for r in out)),'finite',all(r['all_numeric_finite'] for r in out),'gates pass',sum(r['gate'] for r in out),'scope',sorted(set(r['scope'] for r in out)),'chain hash counts',sorted(set(r['chain_hash_count'] for r in out)))
spec=importlib.util.spec_from_file_location('builder',root/'tools/early_tox/build_arms.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
od=pathlib.Path('rebuilt');od.mkdir(exist_ok=True)
for name,rr in zip(['arms','outcomes','arm_schedules','arm_fractions','evidence_families'],b.build(str(root))):
    fn=od/(name+'.csv');b.write_csv(str(fn),rr);print('builder',name,len(rr),'byte match',fn.read_bytes()==(root/'data/early_tox'/fn.name).read_bytes())
o={'report_prefix':'Synthetic','endpoint_class':'grade_3_4_mucositis'}
for ep in ['other: grade 1 mucositis','other: grade 2 mucositis','max_grade_ge3_mucositis']:
    r=dict(report_file='Synthetic_paper.pdf',events=8,denominator=30,endpoint=ep)
    print('endpoint',repr(ep),'problem',repr(b.row_problem(r,o,30,8)))
print('Raw chain bytes were not checked; manifests name a separate Linux archive.')
