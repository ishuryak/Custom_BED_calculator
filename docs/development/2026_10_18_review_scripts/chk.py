import math
def K(n,d,ab=10): return n*d*(1+d/ab)
def hmax(x): return max(0.0,x)
def smax(x,k=3): z=k*x; return x if z>30 else math.log1p(math.exp(z))/k
def di(n,d,T,a,lam,Tk,g=0,mx=hmax): return K(n,d)-lam*mx(T-Tk)/a - g*T/a
def dd_fit(n,d,T,a,gam,C,g=0,mx=hmax):
    k=K(n,d); p=1-math.exp(-a*k/T); on=C*T/(a*k)
    return k-gam*p*mx(T-on)/a - g*T/a, on, gam*p, gam*p/a
# trial pairs from Shuryak 2018 Table 1 (avg schedules) : (exp, ctrl)
pairs={'GORTEC':((31.5,2.0,22),(35,2.0,48)),'TROG(Poulsen)':((33,1.8,24),(35,2.0,50)),
 'CHART':((36,1.5,12),(33,2.0,45)),'DAHANCA':((33.5,2.0,39),(33.5,2.0,46)),
 'CAIR':((40,1.9,40),(40,1.9,54)),'RTOG CB':((42,1.7,43),(35,2.0,50)),'RTOG HF':((68,1.2,50),(35,2.0,50)),
 'EORTC HF':((70,1.2,51),(35,2.0,51)),'HYPNO':((20,2.75,26),(33,2.0,38))}
print("== 2018 FITTED models (Table 2), with g=0.0154 and without ==")
for g in (0.0154,0.0):
  print("g=",g)
  for nm,(e,c) in pairs.items():
    dDI=di(*e,0.069,0.035,28.6,g)-di(*c,0.069,0.035,28.6,g)
    de=dd_fit(*e,0.224,1.17,14.5,g); dc=dd_fit(*c,0.224,1.17,14.5,g)
    print(f"  {nm:14s} dDI {dDI:+6.2f}  dDD {de[0]-dc[0]:+6.2f}   DD onset exp {de[1]:5.1f} ctrl {dc[1]:5.1f}")
r=dd_fit(35,2,47,0.224,1.17,14.5)
print("2018 DD at 35x2/47: onset %.1f d, rate %.3f/d (Tp %.2f d), loss %.2f Gy10/d"%(r[1],r[2],math.log(2)/r[2],r[3]))
print("2018 DI loss %.2f Gy10/d"%(0.035/0.069))
# 55/20 vs 70/35 under fitted models
for g in (0.0154,0):
  a=di(20,2.75,26,0.069,0.035,28.6,g)-di(35,2,47,0.069,0.035,28.6,g)
  b=dd_fit(20,2.75,26,0.224,1.17,14.5,g)[0]-dd_fit(35,2,47,0.224,1.17,14.5,g)[0]
  print("55/20 vs 70/35 fitted g=%g: dDI %+.2f dDD %+.2f"%(g,a,b))
