from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.special import gammaln
from scipy.optimize import minimize_scalar

HERE=Path(__file__).resolve().parent
SRC=HERE/'software/source'
r=pd.read_csv(SRC/'data/refit_trials_2026_10_25c_main.csv')
n=r.n_patients.to_numpy(int)
y=np.rint(n*r.tcp_actual_pct.to_numpy()/100).astype(int)
constant=float((gammaln(n+1)-gammaln(y+1)-gammaln(n-y+1)).sum())
record=json.loads((SRC/'docs/refit_explainer/main_free.json').read_text())
p=record['di']
t=r.total_time_days.to_numpy()
K=r.total_dose_gy.to_numpy()+r.sum_d2_gy2.to_numpy()/10
offset=-p['alpha']*K+p['lam']*np.logaddexp(0,3*(t-p['third']))/3+p['g']*t
def nll(c, idx):
    eta=c+offset[idx]
    logp=-np.exp(eta)
    log1mp=np.log(-np.expm1(logp))
    return float(-(y[idx]*logp+(n[idx]-y[idx])*log1mp).sum())
unit_fits=[minimize_scalar(nll,args=(r.ref_no.to_numpy()==unit,),bounds=(-5,10),method='bounded',options={'xatol':1e-12}) for unit in r.ref_no.unique()]
old_nll=float(sum(x.fun for x in unit_fits))
dr=pd.read_csv(SRC/'draws/tumor_di_primary.csv')
ss=json.loads((SRC/'fit/output_refit_2026_10_25c_main/di_primary/summary.json').read_text())
regimens=[(20,2.75,26),(15,3.1,19),(17,3.05,23),(25,2.4,33),(30,65/30,40),(33,2,45),(30,2,40),(33,2,38),(35,2,40),(68,1.2,47)]
def E(m,d,T):
    return m*d+m*d*d/dr.ab-(dr['lambda']*np.logaddexp(0,3*(T-dr.Tk))/3+dr.g*T)/dr.alpha
anchor=E(35,2,47)
comparisons=[]
for triple,archived in zip(regimens,ss['contrasts_vs_anchor']):
    delta=E(*triple)-anchor
    comparisons.append({'regimen':archived['regimen'],'current_export_median':float(np.median(delta)),
                        'archived_8000_median':archived['median'],'difference':float(np.median(delta)-archived['median'])})
ans={'source':'79e9c812f91482fa57ffa1dddabe89dc329a6303',
     'n_arms':len(r),'n_units':len(r.ref_no.unique()),'time_range':[float(min(t)),float(max(t))],
     'log_binomial_constant':constant,'old_nll_record':p['neg_log_lik'],'old_nll_reprofiled_intercepts':old_nll,
     'full_nll':p['neg_log_lik']-constant,'full_log_likelihood':-p['neg_log_lik']+constant,
     'incorrect_plan_addition':p['neg_log_lik']+constant,
     'export_rows':len(dr),'distinct_di_rows':len(dr.drop_duplicates(subset=['alpha','ab','lambda','Tk','g'])),
     'tumor_T4_comparisons':comparisons,
     'max_abs_median_difference':max(abs(x['difference']) for x in comparisons),
     'T4_exceeds_0_1':sum(abs(x['difference'])>0.1 for x in comparisons)}
(HERE/'plan_numeric_checks.json').write_text(json.dumps(ans,indent=2))
print(json.dumps(ans,indent=2))
