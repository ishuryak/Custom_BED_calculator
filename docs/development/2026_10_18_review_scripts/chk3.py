import math
def K(n,d,ab=10): return n*d*(1+d/ab)
def smax(x,k=3): z=k*x; return x if z>30 else math.log1p(math.exp(z))/k
def T_incl(m): w=(m-1)//5; return 7*w+(m-1)%5+1
def early(n,d,T): return K(n,d)-math.log(2)/2.5*smax(T-1-7)/.35   # elapsed convention (D1)
def late(n,d,ab=3): return n*d*(1+d/ab)
def run(a,ldi,Tk,ldd,label):
  def di(n,d,T): return K(n,d)-ldi*smax(T-Tk)/a
  n0,d0,T0=35,2,47; k0=K(n0,d0); p0=1-math.exp(-a*k0/T0); y=ldi*smax(T0-Tk)/(ldd*p0)
  on=T0-math.log(math.expm1(3*y))/3; ca=k0*on/T0
  def dd(n,d,T):
    k=K(n,d); p=1-math.exp(-a*k/T); return k-ldd*p*smax(T-ca*T/k)/a
  H=(20,2.75,26); eH,lH=early(*H),late(20,2.75)
  print(f"\n{label}: HYPNO early {eH:.1f} late {lH:.1f}; DI {di(*H):.1f} DD {dd(*H):.1f}")
  print(" m   d     T  early  late   dDI_vsHYPNO dDD_vsHYPNO  dDI_vs70/35 dDD_vs70/35")
  for m in range(12,20):
    T=T_incl(m); best=None
    for i in range(0,200):
      d=2.75+0.05*i
      if early(m,d,T)<=eH+1e-9 and late(m,d)<=lH+1e-9: best=d
    if best is None: continue
    d=best
    print(f"{m:3d} {d:4.2f} {T:3d} {early(m,d,T):5.1f} {late(m,d):6.1f}   {di(m,d,T)-di(*H):+6.2f}     {dd(m,d,T)-dd(*H):+6.2f}      {di(m,d,T)-di(35,2,47):+6.2f}     {dd(m,d,T)-dd(35,2,47):+6.2f}")
run(.2,.2,28,.5,"v1 defaults")
run(.35,.23,28,.40,"D18 proposal")
