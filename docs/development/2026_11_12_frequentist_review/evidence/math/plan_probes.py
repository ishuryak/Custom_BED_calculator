"""Bounded numerical checks for the frequentist PLAN. No model refitting or original writes.
python -B plan_probes.py ORIGINAL_PROJECT
Requires numpy/scipy/pandas. Outputs beside this script.
"""
from pathlib import Path
import sys,types,importlib.util,math,json,hashlib,shutil
import numpy as np
import pandas as pd
from scipy.special import gammaln,logsumexp
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.stats import chi2,gamma,binom
repo=Path(sys.argv[1]).resolve();out=Path(__file__).resolve().parent
stub=types.ModuleType('cmdstanpy')
def refuse(*a,**k):raise RuntimeError('No MCMC or compilation allowed')
stub.CmdStanModel=refuse;stub.cmdstan_path=refuse;sys.modules['cmdstanpy']=stub
sp=importlib.util.spec_from_file_location('s5',repo/'tools/early_tox/step5_fit.py');s5=importlib.util.module_from_spec(sp);sp.loader.exec_module(s5)
early,fams=s5.load_rows();tum=pd.read_csv(repo/'data/refit_trials_2026_10_25c_main.csv')
result={'target':str(repo),'plan_files':{}}
for f in ('2026_11_12_frequentist_comparison_plan.md','2026_11_12_frequentist_comparison_plan_audit.md'):
    p=repo/'docs/development'/f
    result['plan_files'][f]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
    shutil.copyfile(p,out/f)

def lb(y,n,eta,link):
    eta=np.asarray(eta)
    if link=='logit':lp=-np.logaddexp(0,-eta);lq=-np.logaddexp(0,eta)
    else:
        t=np.exp(np.minimum(eta,700));lp=-t
        with np.errstate(divide='ignore'):lq=np.where(eta < -36,eta,np.log(-np.expm1(-t)))
    return gammaln(n+1)-gammaln(y+1)-gammaln(n-y+1)+y*lp+(n-y)*lq
def gh(logf,n):
    x,w=np.polynomial.hermite.hermgauss(n)
    return float(logsumexp(np.log(w)+logf(np.sqrt(2)*x))-.5*math.log(math.pi))
def oracle(logf):
    h=lambda z:float(np.asarray(logf(z)-.5*z*z-.5*math.log(2*math.pi)).item())
    fit=minimize_scalar(lambda z:-h(z),bounds=(-20,20),method='bounded',options={'xatol':1e-12})
    mx=h(fit.x)
    val,err=quad(lambda z:math.exp(h(z)-mx),-np.inf,np.inf,epsabs=1e-12,epsrel=1e-12,limit=300)
    return float(mx+math.log(val)),float(err/val)
def agh(logf,n):
    h=lambda z:np.asarray(logf(z))-.5*np.asarray(z)**2-.5*math.log(2*math.pi)
    mode=minimize_scalar(lambda z:-float(h(z).item()),bounds=(-20,20),method='bounded',options={'xatol':1e-12}).x
    eps=1e-4
    curv=float((h(mode+eps)-2*h(mode)+h(mode-eps)).item()/eps**2)
    scale=1/math.sqrt(-curv);x,w=np.polynomial.hermite.hermgauss(n)
    return float(math.log(math.sqrt(2)*scale)+logsumexp(np.log(w)+h(mode+math.sqrt(2)*scale*x)+x*x))

checks=[]
tr=tum.iloc[tum.n_patients.argmax()];N=int(tr.n_patients);y=round(N*tr.tcp_actual_pct/100);eta=math.log(-math.log(y/N))
er=max(early,key=lambda r:r['N']);Ne,ye=er['N'],er['y'];etae=math.log(ye/(Ne-ye))
for label,nn,yy,mu,link in (('tumor largest arm',N,y,eta,'cloglog_control'),(er['arm_id'],Ne,ye,etae,'logit')):
 for sig in (.05,.25,1.):
    f=lambda z:lb(yy,nn,mu+sig*z,link)
    exact,relerr=oracle(f)
    checks.append({'case':label,'N':nn,'y':yy,'eta':mu,'sigma':sig,'link':link,'quad_loglik':exact,'quad_estimated_relative_error':relerr,'adaptive_GH30_error':agh(f,30)-exact,
                   'gh':{str(n):{'loglik':gh(f,n),'error':gh(f,n)-exact}for n in (20,30,60,120)}})
result['single_arm_quadrature']=checks
# Outer integral at sigma_arm=0 is a required nested-model boundary and is independently tractable.
fr=[r for r in early if r['family']=='RTOG0129'];Nf=np.array([r['N']for r in fr]);yf=np.array([r['y']for r in fr]);mu=np.log(yf/(Nf-yf))
outer=[]
for tau in (.5,1.,2.):
    def f(z):return np.sum(lb(yf[:,None],Nf[:,None],mu[:,None]+tau*np.atleast_1d(z)[None,:],'logit'),axis=0)
    exact,relerr=oracle(f)
    outer.append({'family':'RTOG0129','N':Nf.tolist(),'y':yf.tolist(),'eta':mu.tolist(),'tau':tau,'sigma_arm':0,
                  'quad_loglik':exact,'quad_estimated_relative_error':relerr,'adaptive_GH20_error':agh(f,20)-exact,'gh':{str(n):{'loglik':gh(f,n),'error':gh(f,n)-exact}for n in (20,40,80,160)}})
result['outer_family_quadrature']=outer
# T4 recomputes four decision log ORs from the same rounded stored parameter rows.
u=repo/'docs/development/2026_11_10_early_toxicity/step5_review/overdispersed/U02/draws.csv.gz';dr=pd.read_csv(u)
ref=s5.schedule_terms(*s5.ev.regular_schedule(35,2),7,3.2);errs=[]
for k,(m,d) in enumerate(s5.DECISION,1):
    dt=np.subtract(s5.schedule_terms(*s5.ev.regular_schedule(m,d),7,3.2),ref)
    pred=dr.kappa*(dt[0]+dt[1]/dr.ab-dr.r*dt[2]);stored=dr[f'lor_dec[{k}]'];diff=pred-stored
    errs.append({'schedule':[m,d],'terms':dt.tolist(),'max_abs_error':float(np.max(np.abs(diff))),'n_exceed_1e_minus6':int(np.sum(np.abs(diff)>1e-6)),
                 'row_worst':int(np.argmax(np.abs(diff)))+1,'max_relative_error':float(np.max(np.abs(diff)/np.maximum(np.abs(stored),1e-12)))})
result['early_T4_rounding']=errs
# Even observing the 22 Gaussian trial effects perfectly does not yield2% variance-scale recovery.
J=22
result['T5_oracle_family_SD']={'J':J,'assumption':'22 iid Normal(mu,tau^2) effects observed exactly, unknown mean, MLE SD',
    'tau_hat_over_tau_95_interval':np.sqrt(chi2.ppf([.025,.975],J-1)/J).tolist(),
    'probability_within_2_percent':float(chi2.cdf(J*1.02**2,J-1)-chi2.cdf(J*.98**2,J-1))}
result['flat_prior_skew_counterexample']={'model':'one exponential waiting time y=1; likelihood lambda*exp(-lambda), flat prior on positive rate',
    'mle_and_posterior_mode':1.,'posterior_median':float(gamma.ppf(.5,a=2,scale=1))}
result['percentile_order_MC']={str(B):{'expected_order_at_025':B*.025,'95_binomial_count_range':binom.ppf([.025,.975],B,.025).tolist()}for B in (500,1000)}
# Identify which previous independent Bayesian calculations remain applicable.
previous=out.parents[1]/'review_2026_11_11_v46/software/source'
files=('tools/early_tox/evaluator.py','tools/early_tox/step5_model.stan','data/early_tox/arm_fractions.csv','data/early_tox/outcomes.csv','draws/early_fit_draws.csv')
result['previous_verified_math_inputs']={f:{'same_normalized_content':((repo/f).read_text().replace('\r\n','\n')==(previous/f).read_text().replace('\r\n','\n'))if(previous/f).exists()else None,'sha256_current':hashlib.sha256((repo/f).read_bytes()).hexdigest()}for f in files}
if (previous/'tools/early_tox/step5_fit.py').exists():
    oldspec=importlib.util.spec_from_file_location('old_s5',previous/'tools/early_tox/step5_fit.py');old=importlib.util.module_from_spec(oldspec);oldspec.loader.exec_module(old)
    oldrows,oldfamilies=old.load_rows()
    result['current_primary_stan_input_equals_previous_verified']=s5.stan_data(early,fams,7,{'use_od':1})==old.stan_data(oldrows,oldfamilies,7,{'use_od':1})
(out/'plan_probes.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
