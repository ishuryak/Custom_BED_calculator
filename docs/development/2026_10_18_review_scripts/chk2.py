import math
def K(n,d,ab=10): return n*d*(1+d/ab)
def smax(x,k=3): z=k*x; return x if z>30 else math.log1p(math.exp(z))/k
def hmax(x): return max(0,x)
def T_incl(m): w=(m-1)//5; return 7*w+(m-1)%5+1
def early(n,d,T,ab=10,a=.35,Tp=2.5,Tk=7,mx=smax): return K(n,d,ab)-math.log(2)/Tp*mx(T-Tk)/a
print("Fowler 2003: 70/35 T=46 ->%.2f ; CAIR 35x2 T=34 -> %.2f"%(early(35,2,46,mx=hmax),early(35,2,34,mx=hmax)))
print("App now (incl T=47): 70/35 -> %.2f ; 55/20 T=26 -> %.2f ; 66/33 T=38 -> %.2f"%(early(35,2,47),early(20,2.75,26),early(33,2,38)))
# D18 settings
a,ldi,Tk,ldd=.35,.23,28,.40
def di(n,d,T,mx=smax): return K(n,d)-ldi*mx(T-Tk)/a
def cal(ref=(35,2,47)):
    n,d,T=ref; k=K(n,d); p=1-math.exp(-a*k/T); y=ldi*smax(T-Tk)/(ldd*p)
    x0=math.log(math.expm1(3*y))/3; on=T-x0; return k*on/T,on
ca,on=cal()
def dd(n,d,T):
    k=K(n,d); p=1-math.exp(-a*k/T); o=ca*T/k; return k-ldd*p*smax(T-o)/a,o,ldd*p
print("D18: ref DI %.2f, DD onset ref %.1f"%(di(35,2,47),on))
for nm,(n,d,T) in {'55/20':(20,2.75,26),'66/33 6/wk':(33,2,38),'HYHOPE':(15,3.1,19),'GORTEC':(31.5,2,22),'TROG':(33,1.8,24)}.items():
    r=dd(n,d,T); print(f"  {nm:10s} dDI {di(n,d,T)-di(35,2,47):+.2f} dDD {r[0]-dd(35,2,47)[0]:+.2f} onset {r[1]:.1f} Tp {math.log(2)/r[2]:.2f}")
# tumor time-convention effect for 55/20 vs 70/35 under D18 DI
print("D18 dDI 55/20 inclusive %.2f vs elapsed %.2f"%(di(20,2.75,26)-di(35,2,47), di(20,2.75,25)-di(35,2,46)))
# smooth max bias at onset
print("smooth max at 0: %.3f d"%smax(0))
# reference-dependence of DD (current defaults)
