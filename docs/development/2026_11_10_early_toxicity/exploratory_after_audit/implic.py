import numpy as np, sys
def ot(m): w=(m-1)//5; return 7*w+(m-1)%5+1
def smax(x,k=3.0): return np.where(k*x>30,x,np.log1p(np.exp(np.minimum(k*x,30)))/k)
def early(m,d,Te,Tk,r,ab): return m*d*(1+d/ab)-r*smax(Te-1-Tk)
ds=np.arange(2.05,8.0001,0.05); cur=dict(Tk=7.0,r=np.log(2)/2.5/0.35,ab=10.0)
def q(dm):  # infeasible draws (no dose >= 2.05 Gy) count as below the grid (2.00)
    x=np.where(np.isnan(dm),2.0,dm); p=np.percentile(x,[2.5,50,97.5])
    f=lambda v: "none" if v<2.04 else f"{v:.2f}"
    return f"{f(p[1])} [{f(p[0])}-{f(p[2])}] infeas {np.mean(np.isnan(dm)):.0%}"
for name in sys.argv[1:]:
    X=np.load(f'draws_{name}.npy')[::5]; B0,Tk,r,s,ab=X.T
    std=early(35,2,47,Tk,r,ab); hy=early(20,2.75,26,Tk,r,ab); hyc=early(33,2,38,Tk,r,ab)
    print(f"== {name}: standard early BED<=0 in {np.mean(std<=0):.1%} of draws")
    print(f"   HYPNO 55/20 minus standard: median {np.median(hy-std):+.1f} [{np.percentile(hy-std,2.5):+.1f},{np.percentile(hy-std,97.5):+.1f}] Gy10, passes +5 in {np.mean(hy-std<=5):.0%}")
    print(f"   HYPNO 55/20 minus 66/33(6/wk): median {np.median(hy-hyc):+.1f} [{np.percentile(hy-hyc,2.5):+.1f},{np.percentile(hy-hyc,97.5):+.1f}]")
    for m in [13,15,16,20,25,30]:
        T=ot(m); E=early(m,ds[:,None],T,Tk,r,ab)-std
        dm=np.array([ds[E[:,j]<=5].max() if (E[:,j]<=5).any() else np.nan for j in range(E.shape[1])])
        c=early(m,ds,T,**cur)-early(35,2,47,**cur)
        print(f"   m={m:2d} T={T:2d} app {ds[c<=5].max():.2f}  refit dmax {q(dm)}")
