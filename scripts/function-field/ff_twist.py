# S_m over F_q[T], m = T^3 - T + 1 (q=3): exact Euler product, Weil-zero oscillations test
import sys, itertools
sys.argv=['ff','3','8','4']
exec(open('ff_jacobi.py').read().split("if __name__ == '__main__':")[0])
import mpmath; mpmath.mp.dps=60
from sympy import mobius, divisors
m=[1,2,0,1]          # 1 - T + T^3 over F_3 (coeffs low->high, -1 = 2)
D0=8; N=int(sys.argv_n) if hasattr(sys,'argv_n') else 40
def chi(P):
    r=pmod(m,P); return legendre(r,P) if r else 0
def Nirr(d): return sum(int(mobius(k))*q**(d//k) for k in divisors(d))//d
# 1) direct character sums over primes of degree <= D0
pi={d:sum(chi(P) for P in irr[d]) for d in range(1,D0+1)}
dm=len(m)-1
def even_sum(k): return Nirr(k)-(1 if k==dm else 0)
def s(d,pi): return sum(k*(pi[k] if (d//k)%2 else even_sum(k)) for k in divisors(d))
# L(w)=exp(sum s_d w^d/d)
Ls=[mpmath.mpf(1)]+[0]*D0
lg=[0]+[mpmath.mpf(s(d,pi))/d for d in range(1,D0+1)]
for n in range(1,D0+1): Ls[n]=sum(k*lg[k]*Ls[n-k] for k in range(1,n+1))/n
print("L(w) coefficients from primes up to degree 8:", [mpmath.nstr(c,6) for c in Ls])
L=[int(mpmath.nint(c)) for c in Ls[:3]]; assert all(abs(c)<1e-20 for c in Ls[3:])
al=mpmath.polyroots([L[2],L[1],L[0]])           # roots w of L(w)=0
print("  L(w) =",L,"; |zeros| =",[mpmath.nstr(abs(w),10) for w in al]," vs q^-1/2 =",mpmath.nstr(mpmath.mpf(q)**-0.5,10))
# 2) pi_chi(d) for all d from L: s_d = -sum alpha^d, alpha=1/w
for d in range(D0+1,N+1):
    sd=int(mpmath.nint(-sum((1/w)**d for w in al).real))
    pi[d]=(sd-sum(k*(pi[k] if (d//k)%2 else even_sum(k)) for k in divisors(d) if k<d))//d
for d in range(1,D0+1): assert s(d,pi)==int(mpmath.nint(-sum((1/w)**d for w in al).real))
# 3) exact Euler product
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
def G(d,c):  # local series in u for one prime of degree d with chi=c
    Qd=q**d; g=[0]*(N+1); g[0]=1
    for e in range(1,N//d+1):
        g[e*d]=(Qd-1)*Qd**(e-1) if c==0 else (-Qd**(e-1) if e%2 else (Qd-1-c)*Qd**(e-1))
    return g
Z=[1]+[0]*N
for d in range(1,N+1):
    nz=1 if d==dm else 0; npl=(Nirr(d)-nz+pi[d])//2; nmi=(Nirr(d)-nz-pi[d])//2
    Z=smul(Z,spow(G(d,1),npl)); Z=smul(Z,spow(G(d,-1),nmi))
    if nz: Z=smul(Z,G(d,0))
# 4) brute-force check of A(n) for small n
def Sm(f):
    fac=factor(f); n=len(f)-1; tot=0
    for x in itertools.product(range(q),repeat=n):
        a=pmod(sub(mul(list(x),list(x)),m),f); j=1
        for P,e in fac:
            j*=(legendre(a,P) if a else 0)**e
            if j==0: break
        tot+=j
    return tot
for n in range(1,5):
    b=sum(Sm(f) for f in monics(n)); print(f"  A({n}): Euler product {Z[n]}, brute force {b}", "OK" if b==Z[n] else "MISMATCH")
# 5) remove main poles, then remove L(qu^2)
M=[Z[n]-(q**3*Z[n-2] if n>=2 else 0) for n in range(N+1)]
Lq=[1,0,L[1]*q,0,L[2]*q*q]
K=[sum(Lq[i]*M[n-i] for i in range(5) if n-i>=0) for n in range(N+1)]
import pickle; pickle.dump(M,open("/tmp/twist_M.pkl","wb"))
print("\n growth  log_q|c_n|/n :  M = Z(1-q^3u^2)   K = M*L(qu^2)")
for n in range(16,N+1,4): print(f"   n={n}:   {float(mpmath.log(abs(M[n]))/mpmath.log(q))/n:.4f}        {float(mpmath.log(abs(K[n]))/mpmath.log(q))/n if K[n] else float('nan'):.4f}")
# 6) predicted oscillation: poles u^2 = w_j/q
poles=[r for w in al for r in (mpmath.sqrt(w/q),-mpmath.sqrt(w/q))]
print("\n predicted poles: |u| =",mpmath.nstr(abs(poles[0]),8)," = q^-3/4 =",mpmath.nstr(mpmath.mpf(q)**-0.75,8),"; args/pi =",[mpmath.nstr(mpmath.arg(p)/mpmath.pi,5) for p in poles])
# fit M_n ~ sum_j c_j p_j^{-n} on n in [N-12,N] (4 unknowns + check)
A=mpmath.matrix([[p**(-n) for p in poles] for n in range(N-11,N-7)]); bb=mpmath.matrix([M[n] for n in range(N-11,N-7)])
c=mpmath.lu_solve(A,bb)
print(" fit with predicted poles (4 coeffs from n=%d..%d), relative error on later n:"%(N-11,N-8))
for n in range(N-7,N+1):
    pr=sum(c[j]*poles[j]**(-n) for j in range(4)).real
    print(f"   n={n}: M_n/q^(3n/4) = {float(M[n]/mpmath.mpf(q)**(0.75*n)):+.5f}   prediction {float(pr/mpmath.mpf(q)**(0.75*n)):+.5f}")
