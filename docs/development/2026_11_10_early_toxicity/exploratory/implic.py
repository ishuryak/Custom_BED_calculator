import numpy as np
def ot(m): w=(m-1)//5; return 7*w+(m-1)%5+1
def smax(x,k=3.0): return np.where(k*x>30,x,np.log1p(np.exp(np.minimum(k*x,30)))/k)
def early(m,d,T,Tk,r,ab): return m*d*(1+d/ab)-r*smax(T-1-Tk)
ds=np.arange(2.05,8.0001,0.05)
cur=dict(Tk=7.0,r=np.log(2)/2.5/0.35,ab=10.0)
for name in ['weak','informed']:
    X=np.load(f'draws_{name}.npy')[::5]; B0,Tk,r,s,ab=X.T
    std=early(35,2,47,Tk,r,ab)
    print(f"== {name}: P(intol) at 70/35: median {np.median(1/(1+np.exp(-(std-B0)/s))):.3f}")
    print(" m  T  dmax_now  dmax(+5) median [2.5-97.5%]   P(intol) at dmax_now median [95%]")
    for m in [13,15,16,20,25,30]:
        T=ot(m); E=early(m,ds[:,None],T,Tk,r,ab)-std   # grid x draws
        dm=np.array([ds[np.where(E[:,j]<=5)[0]].max() if (E[:,j]<=5).any() else np.nan for j in range(E.shape[1])])
        c=early(m,ds,T,**cur)-early(35,2,47,**cur); dnow=ds[c<=5].max()
        p=1/(1+np.exp(-(early(m,dnow,T,Tk,r,ab)-B0)/s))
        print(f"{m:2d} {T:2d}  {dnow:.2f}   {np.nanmedian(dm):.2f} [{np.nanpercentile(dm,2.5):.2f}-{np.nanpercentile(dm,97.5):.2f}]   {np.median(p):.2f} [{np.percentile(p,2.5):.2f}-{np.percentile(p,97.5):.2f}]")
