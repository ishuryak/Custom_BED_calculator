import pandas as pd, numpy as np, pickle, sys
x = pd.read_excel(sys.argv[1], header=None)
S=[]; i=0
while i < len(x):
    if str(x.iat[i,0]).strip()=="Reference":
        ref=str(x.iat[i,1]).strip()
        d=pd.to_numeric(x.iloc[i+1,1:],errors='coerce').dropna().values
        t=pd.to_numeric(x.iloc[i+2,1:],errors='coerce').dropna().values
        cl=str(x.iat[i+3,1]).strip()
        assert len(d)==len(t), (i,ref)
        S.append(dict(ref=ref,d=d,t=t,cls=cl)); i+=4
    else: i+=1
print(len(S), pd.Series([s['cls'] for s in S]).value_counts().to_dict())
def bed(s, ab, th=3.2):
    d,t=s['d'],s['t']; mu=np.log(2)/th
    dt=np.abs(t[:,None]-t[None,:]); W=np.exp(-mu*dt)
    return d.sum()+ (d@W@d)/ab
for s in S:
    s['T']=np.floor((s['t'].max()-s['t'].min())/24)+1   # inclusive days
    s['n']=len(s['d']); s['D']=s['d'].sum(); s['y']=1 if 'intol' in s['cls'].lower() else 0
rows=[(s['ref'],s['n'],round(s['D'],1),round(s['D']/s['n'],2),int(s['T']),round(bed(s,10),1),s['cls']) for s in S]
df=pd.DataFrame(rows,columns=['ref','n','D','dmean','T_incl','BEDck10','cls'])
print(df.sort_values('T_incl').to_string())
pickle.dump(S,open(sys.argv[2],'wb'))
