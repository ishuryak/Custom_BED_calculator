import math
def K(n,d,ab=10): return n*d*(1+d/ab)
def s(x): z=3*x; return x if z>30 else math.log1p(math.exp(z))/3
def inv(y): z=3*y; return y if z>30 else math.log(math.expm1(z))/3
def acute(n,d,Tincl): return K(n,d)-math.log(2)/2.5*s(Tincl-1-7)/.35
def model(a,ldi,Tk,ldd,anchor=(35,2,47),g=0):
  n,d,T=anchor; k=K(n,d); p=-math.expm1(-a*k/T); y=ldi*s(T-Tk)/(ldd*p); on=T-inv(y); ca=k*on/T
  di=lambda n,d,T: K(n,d)-ldi*s(T-Tk)/a-g*T/a
  def dd(n,d,T):
    k=K(n,d); p=-math.expm1(-a*k/T); return k-ldd*p*s(T-ca*T/k)/a-g*T/a
  return di,dd,on,ca
L=model(.2,.2,28,.5); D=model(.35,.23,28,.40)
for nm,M in (('legacy',L),('D18',D)):
  di,dd,on,ca=M; print(nm,'onset %.8f C/a %.8f'%(on,ca))
  for r in [(35,2,47),(20,2.75,26),(33,2,38),(15,3.1,19),(22,2.75,30)]:
    print('  ',r,'%.8f %.8f acute %.8f'%(di(*r),dd(*r),acute(*r)))
print('late 55/20 %.8f'%(55*(1+2.75/3)))
# legacy coupled-anchor 55/20 DD with ref 35x2/47, 33x2/38, 20x2.75/26
for ref in [(35,2,47),(33,2,38),(20,2.75,26)]:
  print('coupled ref',ref,'%.4f'%model(.2,.2,28,.5,anchor=ref)[1](20,2.75,26))
# Withers r=5
r=5; k=K(35,2,5); print('Withers r5 loss %.4f refBED %.4f'%(0.6*(1+2/r), k-0.6*(1+2/r)*s(47-28)))
print('Fowler2010 60x1.3 T39 %.8f T42 %.8f'%(K(60,1.3)-math.log(2)/2.5*max(0,39-7)/.35, K(60,1.3)-math.log(2)/2.5*max(0,42-7)/.35))
print('matched ldd %.12f'%(0.2/(-math.expm1(-.2*84/47))))
g=0.0154; print('g ref %.4f 26d %.4f'%(g*47/.2,g*26/.2))
print('late2 20x2.9 %.4f vs %.4f; late3 %.4f vs %.4f'%(58*(1+2.9/2),70*2,58*(1+2.9/3),70*(1+2/3)))
