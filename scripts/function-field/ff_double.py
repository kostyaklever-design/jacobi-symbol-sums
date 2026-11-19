# double-zero case of Theorem C_q (q=5, m = T^5+T+1): pole of order 2 gives n*q^(3n/4); usage: python3 ff_double.py q "m coeffs low->high" N
import sys, itertools, pickle
qq=int(sys.argv[1]); mm=[int(c) for c in sys.argv[2].split(',')]; N=int(sys.argv[3])
N=int(sys.argv[3]); D0=len(mm)+1
sys.argv=['ff',str(qq),str(D0),str(N)]
exec(open('ff_jacobi.py').read().split("if __name__ == '__main__':")[0])
import mpmath; mpmath.mp.dps=40
from sympy import mobius, divisors
m=mm; dm=len(m)-1
fm=factor(m); assert all(e==1 for _,e in fm), "m not squarefree"
mdeg=[len(P)-1 for P,_ in fm]
def chi(P):
    r=pmod(m,P); return legendre(r,P) if r else 0
NI={}
def Nirr(d):
    if d not in NI: NI[d]=sum(int(mobius(k))*q**(d//k) for k in divisors(d))//d
    return NI[d]
def even_sum(k): return Nirr(k)-mdeg.count(k)
pi={d:sum(chi(P) for P in irr[d]) for d in range(1,D0+1)}
def s(d): return sum(k*(pi[k] if (d//k)%2 else even_sum(k)) for k in divisors(d))
Ls=[mpmath.mpf(1)]+[0]*D0; lg=[0]+[mpmath.mpf(s(d))/d for d in range(1,D0+1)]
for n in range(1,D0+1): Ls[n]=sum(k*lg[k]*Ls[n-k] for k in range(1,n+1))/n
L=[int(mpmath.nint(c)) for c in Ls]; 
while L[-1]==0: L.pop()
assert all(abs(Ls[i]-L[i])<1e-20 for i in range(len(L))) and all(abs(c)<1e-20 for c in Ls[len(L):])
wz=mpmath.polyroots(L[::-1],maxsteps=200,extraprec=200)
print(f"q={q}, m={m}, factors degrees {mdeg}; L(w) = {L}")
print("  |zeros w|/q^-1/2 =",[mpmath.nstr(abs(w)*mpmath.sqrt(q),12) for w in wz],"  L(q^-1/2) =",mpmath.nstr(mpmath.polyval(L[::-1],q**-0.5),6))
al=[1/w for w in wz]
for d in range(D0+1,N+400):
    sd=int(mpmath.nint(-sum(a**d for a in al).real))
    pi[d]=(sd-sum(k*(pi[k] if (d//k)%2 else even_sum(k)) for k in divisors(d) if k<d))//d
for d in range(1,D0+1): assert s(d)==int(mpmath.nint(-sum(a**d for a in al).real))
# exact series
def smul(a,b):
    r=[0]*(N+1)
    for i,x in enumerate(a):
        if x:
            for j in range(N+1-i): r[i+j]+=x*b[j]
    return r
def spow(a,e):
    r=[1]+[0]*N
    while e:
        if e&1: r=smul(r,a)
        a=smul(a,a); e>>=1
    return r
def G(d,c):
    Qd=q**d; g=[0]*(N+1); g[0]=1
    for e in range(1,N//d+1):
        g[e*d]=(Qd-1)*Qd**(e-1) if c==0 else (-Qd**(e-1) if e%2 else (Qd-1-c)*Qd**(e-1))
    return g
Z=[1]+[0]*N
for d in range(1,N+1):
    nz=mdeg.count(d); npl=(Nirr(d)-nz+pi[d])//2; nmi=(Nirr(d)-nz-pi[d])//2
    Z=smul(Z,spow(G(d,1),npl)); Z=smul(Z,spow(G(d,-1),nmi))
    for _ in range(nz): Z=smul(Z,G(d,0))
def Sm(f):
    fac=factor(f); n=len(f)-1; tot=0
    for x in itertools.product(range(q),repeat=n):
        a=pmod(sub(mul(list(x),list(x)),m),f); j=1
        for P,e in fac:
            j*=(legendre(a,P) if a else 0)**e
            if j==0: break
        tot+=j
    return tot
for n in range(1,4):
    b=sum(Sm(f) for f in monics(n)); print(f"  A({n}): Euler {Z[n]}  brute {b}  {'OK' if b==Z[n] else 'MISMATCH'}")
M=[Z[n]-(q**3*Z[n-2] if n>=2 else 0) for n in range(N+1)]
# analytic side
import sympy as sp
v_,Q_=sp.symbols('v Q'); dif={}
for ch in (1,-1):
    num=1-v_-(1+ch)*Q_*v_**2; den=(1-v_)*(1-Q_*v_**2)*(1-ch*Q_*v_**2)*(1-ch*Q_**2*v_**4)
    dif[ch]=(sp.lambdify((v_,Q_),sp.expand(num-den),'mpmath'),sp.lambdify((v_,Q_),den,'mpmath'))
DMAX=N+399
def logWp(u):
    t=mpmath.mpc(0)
    for d in range(1,DMAX+1):
        nz=mdeg.count(d); npl=(Nirr(d)-nz+pi[d])//2; nmi=(Nirr(d)-nz-pi[d])//2
        V=u**d; Qd=mpmath.mpf(q)**d
        for ch,cnt in ((1,npl),(-1,nmi)):
            if cnt: f,g=dif[ch]; t+=cnt*mpmath.log1p(f(V,Qd)/g(V,Qd))
    return t
Lp=lambda w: mpmath.polyval(L[::-1],w)
def H(u):
    E=1
    for k in mdeg: Qk=mpmath.mpf(q)**k; E*=(1+Qk*u**k)/(1-Qk*u**(2*k))
    return (1-q*u)*(1-q*q*u*u)/Lp(q*q*u**4)*E*mpmath.exp(logWp(u))
Mfun=lambda u: H(u)/Lp(q*u*u)
dLp=lambda w: mpmath.polyval([i*c for i,c in enumerate(L)][1:][::-1],w)
wu=[]
for w in wz:
    if all(abs(w-x)>1e-8 for x in wu): wu.append(w)
poles=[r for w in wu for r in (mpmath.sqrt(w/q),-mpmath.sqrt(w/q))]

ell=lambda u: Lp(q*u*u)
def d2(f,x,h=mpmath.mpf('1e-6')): return (f(x+h)-2*f(x)+f(x-h))/h**2
mpmath.mp.dps=30
print(f"  {len(poles)} distinct poles; multiplicity of w-zeros: {len(wz)//len(wu)}")
a=[2*H(p)/(d2(ell,p)*p*p) for p in poles]
for p,x in zip(poles,a): print(f"    u={mpmath.nstr(p,7):>26}  predicted a_rho={mpmath.nstr(x,8)}")
NP=96; r0=mpmath.mpf('0.02')
circ={}
for j,p in enumerate(poles):
    circ[j]=[(r0*mpmath.expjpi(2*mpmath.mpf(k)/NP), Mfun(p+r0*mpmath.expjpi(2*mpmath.mpf(k)/NP))) for k in range(NP)]
print("\n  residue at each pole (small circle) divided by (n+1) u^-n  ->  should tend to a_rho:")
for n in (10,20,30):
    row=[]
    for j,p in enumerate(poles):
        res=sum(f*(p+z)**(-n-1)*z for z,f in circ[j])/NP     # (1/2pi i)∮ M u^{-n-1} du
        row.append(mpmath.nstr(-res/((n+1)*p**(-n)),7))
    print(f"   n={n}: "+"  ".join(row))
print("\n  M_n/(n q^(3n/4)):")
for n in range(10,N+1,4): print(f"   n={n}: {mpmath.nstr(M[n]/(n*mpmath.mpf(q)**(0.75*n)),6)}")
