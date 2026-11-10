import pandas as pd, numpy as np, sys
x = pd.read_excel(sys.argv[1], header=None)
S=[];i=0
while i<len(x):
    if str(x.iat[i,0]).strip()=="Reference":
        d=pd.to_numeric(x.iloc[i+1,1:],errors='coerce').dropna().values
        t=pd.to_numeric(x.iloc[i+2,1:],errors='coerce').dropna().values
        S.append(dict(row=i+1,ref=str(x.iat[i,1]).strip(),d=d,t=t,cls=str(x.iat[i+3,1]).strip()));i+=4
    else:i+=1
print(len(S))
for k,s in enumerate(S[104:],105): print(k,s['row'],s['ref'],len(s['d']),s['d'].sum(),s['cls'])
intol=lambda s:'intol' in s['cls'].lower()
part=set(range(108,116))
sel={'all':S,'dropDisp':[s for s in S if s['cls']!='Disp'],
 'dropAllDisp':[s for s in S if not s['cls'].startswith('Disp')],
 'complete_dropDisp':[s for k,s in enumerate(S) if k not in part and s['cls']!='Disp'],
 'complete_dropAllDisp':[s for k,s in enumerate(S) if k not in part and not s['cls'].startswith('Disp')]}
for k,v in sel.items(): print(k,len(v),sum(map(intol,v)))
v=sel['dropDisp']
c={};f={}
for s in v:
    c.setdefault((len(s['d']),round(s['d'].sum(),2),int(np.floor((s['t'].max()-s['t'].min())/24)+1)),[]).append(s)
    f.setdefault((tuple(np.round(s['d'],3)),tuple(np.round(s['t']-s['t'].min(),2))),[]).append(s)
print('coarse',len(c),'full',len(f),'event groups coarse',sum(any(map(intol,g)) for g in c.values()),'full',sum(any(map(intol,g)) for g in f.values()))
for k,g in c.items():
    sub={tuple(np.round(s['t']-s['t'].min(),2)) for s in g}
    if len(sub)>1: print('merged',k,[s['ref'] for s in g])
# BED checks
def early(n,d,Te,ab=10,r=np.log(2)/2.5/0.35,Tk=7): return n*d*(1+d/ab)-r*max(0,Te-Tk)
print(early(35,2,46),early(20,2.75,25),early(33,2,37),early(20,2.75,25)-early(33,2,37))
mu=np.log(2)/3.2
s=[s for s in S if len(s['d'])==35 and abs(s['d'].sum()-70)<.1][0]
W=np.exp(-mu*np.abs(s['t'][:,None]-s['t'][None,:])); print('repair70/35',s['d'].sum()+s['d']@W@s['d']/10, s['ref'])
t=np.arange(5)*24.;d=np.full(5,8.);W=np.exp(-mu*np.abs(t[:,None]-t[None,:]));print('stress',d@W@d/10-(d**2).sum()/10)
p=190/373;q=208/380;dd=p-q;se=np.sqrt(p*(1-p)/373+q*(1-q)/380);print(dd*100,(dd-1.96*se)*100,(dd+1.96*se)*100)
L=np.log(q/(1-q))
for sc in (50,10): print(sc,(1/(1+np.exp(-(L+3.9/sc)))-q)*100)
