# exact pole contribution of 1/L(qu^2) to M = Z_m(1-q^3u^2), m = T^3-T+1, q=3
import sympy as sp, mpmath
from sympy import mobius, divisors
mpmath.mp.dps=50
q=3; dm=3; Qm=q**dm; Lc=[1,3,3]; DMAX=400
v,Q,c=sp.symbols('v Q c')
num=1-v-(1+c)*Q*v**2; den=(1-v)*(1-Q*v**2)*(1-c*Q*v**2)*(1-c*Q**2*v**4)
dif={}
for ch in (1,-1):
    nd=sp.expand((num-den).subs(c,ch)); dd=den.subs(c,ch)
    dif[ch]=(sp.lambdify((v,Q),nd,'mpmath'),sp.lambdify((v,Q),dd,'mpmath'))
def Nirr(d): return sum(int(mobius(k))*q**(d//k) for k in divisors(d))//d
al=[1/w for w in mpmath.polyroots([Lc[2],Lc[1],Lc[0]])]          # inverse zeros of L(w)
def even_sum(k): return Nirr(k)-(1 if k==dm else 0)
pi={}
for d in range(1,DMAX+1):
    sd=int(mpmath.nint(-sum(a**d for a in al).real))
    pi[d]=(sd-sum(k*(pi[k] if (d//k)%2 else even_sum(k)) for k in divisors(d) if k<d))//d
def logWp(u):
    s=mpmath.mpc(0)
    for d in range(1,DMAX+1):
        nz=1 if d==dm else 0; npl=(Nirr(d)-nz+pi[d])//2; nmi=(Nirr(d)-nz-pi[d])//2
        V=u**d; Qd=mpmath.mpf(q)**d
        for ch,cnt in ((1,npl),(-1,nmi)):
            if cnt: f,g=dif[ch]; s+=cnt*mpmath.log1p(f(V,Qd)/g(V,Qd))
    return s
Lw=lambda w: Lc[0]+Lc[1]*w+Lc[2]*w*w
def H(u):  # M(u) = H(u)/L(qu^2)
    return ((1-q*u)*(1-q*q*u*u)/((1-u**3)*(1-Qm*u**6)) /Lw(q*q*u**4) * mpmath.exp(logWp(u))
            * (1-u**3)*(1-Qm**2*u**6)/(1-Qm*u**3))
poles=[r for w in mpmath.polyroots([Lc[2],Lc[1],Lc[0]]) for r in (mpmath.sqrt(w/q),-mpmath.sqrt(w/q))]
res=[]
for p in poles:
    dL=(Lc[1]*q*2*p+Lc[2]*q*q*4*p**3)          # d/du L(qu^2)
    res.append(-H(p)/(dL*p))                       # M_n ~ sum res_j * p_j^-n
    print(f"pole u={mpmath.nstr(p,8)}: H(p)={mpmath.nstr(H(p),8)}  coeff c_j={mpmath.nstr(res[-1],8)}")
# tail stability check
import pickle; pickle.dump([(complex(p),complex(r)) for p,r in zip(poles,res)],open('/tmp/twist_res.pkl','wb'))
