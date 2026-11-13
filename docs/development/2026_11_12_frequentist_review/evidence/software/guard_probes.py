"""Current guard positives/negatives, no sampler. python -B guard_probes.py SOURCE"""
import argparse, importlib.util,json,pathlib,sys,tempfile,types,hashlib,shutil,ast,csv
import numpy as np
import pandas as pd
ap=argparse.ArgumentParser();ap.add_argument('repo',type=pathlib.Path);a=ap.parse_args();root=a.repo.resolve()
def forbidden(*a,**k): raise AssertionError('Sampling disabled')
stub=types.ModuleType('cmdstanpy');stub.__version__='review-stub';stub.CmdStanModel=forbidden;stub.cmdstan_path=lambda:'/stub/cmdstan-2.40.0';sys.modules['cmdstanpy']=stub
def load(n):
    sp=importlib.util.spec_from_file_location(n,root/'tools/early_tox'/f'{n}.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
s5=load('step5_fit');s6=load('step6_validate');b=load('build_arms');ex=load('export_fit_for_app')
rr,ff=s5.load_rows();d=s5.stan_data(rr,ff,7,{'use_od':1});req=s5.required_quantities(d)
summ=pd.DataFrame(dict(Mean=.5,StdDev=.1,R_hat=1.001,ESS_bulk=1000.,ESS_tail=1000.),index=req)
gate=lambda e:s5.gate_from_summary(summ,d,0,e,0)['gate_pass']
print('valid 4-chain gate',gate([.9]*4))
print('NaN at each chain position',[gate([np.nan if i==j else .9 for i in range(4)]) for j in range(4)])
print('Inf at each chain position',[gate([np.inf if i==j else .9 for i in range(4)]) for j in range(4)])
print('wrong chain count accepted',gate([.9]))
print('old incomplete manifest problem',s5.manifest_problem(dict(identity='a'*64,gate_pass=True)))
print('old grade1 endpoint problem',b.row_problem(dict(report_file='Example',denominator=30,events=8,endpoint='other: grade 1 mucositis'),dict(report_prefix='Example',endpoint_class='grade_3_4_mucositis'),30,8))
base=root/'docs/development/2026_11_10_early_toxicity/step5_review/overdispersed'
for unit in ['U02','U04','U15','U17','U19']:
    try:ex.read_unit(str(base/unit));print('export source',unit,'accepted')
    except RuntimeError as e:print('export source',unit,'refused',str(e).replace(str(base), 'SOURCE_UNITS'))
with tempfile.TemporaryDirectory(dir='.') as temp:
    target=pathlib.Path(temp)/'U02';shutil.copytree(base/'U02',target)
    diag=json.loads((target/'diagnostics.json').read_text());man=json.loads((target/'manifest.json').read_text())
    diag.update(rhat=-float('inf'),ess_bulk=float('inf'),ess_tail=float('inf'),ebfmi_min=float('inf'))
    (target/'diagnostics.json').write_text(json.dumps(diag));man['outputs']['diagnostics.json']=hashlib.sha256((target/'diagnostics.json').read_bytes()).hexdigest();(target/'manifest.json').write_text(json.dumps(man))
    try:ex.read_unit(str(target));print('export nonfinite gate metadata accepted True')
    except Exception as e:print('export nonfinite gate refused',type(e).__name__,str(e))
records=list(csv.DictReader((root/'docs/development/2026_11_10_early_toxicity/chain_recheck_2026_11_11.csv').open()))
print('retained recheck rows',len(records),'checked',sum(r['status']=='checked' for r in records),'passed',sum(r['gate_pass']=='True' for r in records))
bad=[r['fit'] for r in records if r['gate_pass']!='True'];print('recorded failures',bad)
energy=[ast.literal_eval(r['ebfmi_per_chain']) for r in records]
print('recheck energy complete finite',all(len(x)==4 and np.isfinite(x).all() for x in energy),'minimum',min(min(x) for x in energy))
print('No raw chains or sampling rechecked in this harness.')

