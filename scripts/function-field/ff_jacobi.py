"""The Jacobi-sum object over F_q[T] (function-field version of S(n)), q odd prime.
S(f) = sum_{x mod f} ((x^2 - 1)/f)  for monic f, Jacobi symbol of F_q[T].
(1) brute-force check of the local formula S(P^e) = -Q^(e-1) (e odd), Q^(e-1)(Q-2) (e even), Q = q^deg P, and multiplicativity;
(2) exact A(n) = sum_{deg f = n} S(f) up to n = N from the Euler product;
(3) main term from the poles u = ±q^(-3/2) and the size of the remainder against q^(3n/4).
usage: python3 ff_jacobi.py q brute_deg N
"""
import sys, itertools, mpmath

q = int(sys.argv[1]) if len(sys.argv) > 1 else 3
BD = int(sys.argv[2]) if len(sys.argv) > 2 else 5
N = int(sys.argv[3]) if len(sys.argv) > 3 else 30

# ---------- polynomial arithmetic over F_q (lists, low -> high) ----------
def trim(a):
    while a and a[-1] % q == 0: a = a[:-1]
    return a
def sub(a, b):
    n = max(len(a), len(b)); return trim([((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % q for i in range(n)])
def mul(a, b):
    if not a or not b: return []
    r = [0]*(len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): r[i+j] = (r[i+j] + x*y) % q
    return trim(r)
def pmod(a, m):            # m monic
    a = trim(list(a)); dm = len(m) - 1
    while len(a) - 1 >= dm and a:
        c = a[-1]; s = len(a) - 1 - dm
        for i in range(len(m)): a[s+i] = (a[s+i] - c*m[i]) % q
        a = trim(a)
    return a
def powmod(a, e, m):
    r = [1]; a = pmod(a, m)
    while e:
        if e & 1: r = pmod(mul(r, a), m)
        a = pmod(mul(a, a), m); e >>= 1
    return r
def monics(d):
    for c in itertools.product(range(q), repeat=d): yield list(c) + [1]

# irreducibles up to degree BD by sieve
irr = {}
for d in range(1, BD + 1):
    irr[d] = []
    for f in monics(d):
        if all(pmod(f, P) for e in range(1, d//2 + 1) for P in irr[e]): irr[d].append(f)
def factor(f):
    out = []; f = list(f)
    for d in range(1, len(f)):
        for P in irr.get(d, []):
            e = 0
            while len(f) > 1 and not pmod(f, P):
                # exact division
                qt = []; r = list(f); dm = len(P) - 1
                quo = [0]*(len(r) - dm)
                while len(r) - 1 >= dm and r:
                    c = r[-1]; s = len(r) - 1 - dm; quo[s] = c
                    for i in range(len(P)): r[s+i] = (r[s+i] - c*P[i]) % q
                    r = trim(r)
                f = trim(quo); e += 1
            if e: out.append((P, e))
    return out
def legendre(a, P):
    Q = q ** (len(P) - 1); v = powmod(a, (Q - 1)//2, P)
    if not v: return 0
    return 1 if v == [1] else -1
def S_brute(f):
    fac = factor(f); n = len(f) - 1; tot = 0
    for x in itertools.product(range(q), repeat=n):
        a = pmod(sub(mul(list(x), list(x)), [1]), f); j = 1
        for P, e in fac:
            l = legendre(a, P) if a else 0
            j *= l ** e
            if j == 0: break
        tot += j
    return tot
def S_local(f):
    r = 1
    for P, e in factor(f):
        Q = q ** (len(P) - 1); r *= -Q**(e-1) if e % 2 else Q**(e-1)*(Q-2)
    return r

if __name__ == '__main__':
    bad = 0; cnt = 0; Ab = {}
    for d in range(1, BD + 1):
        s = 0
        for f in monics(d):
            v = S_brute(f); cnt += 1; bad += (v != S_local(f)); s += v
        Ab[d] = s
    print(f"q={q}: brute force over all {cnt} monic f of degree 1..{BD}: local formula violations = {bad}")
    # ---------- exact A(n) from the Euler product ----------
    def Nirr(d):
        from sympy import mobius, divisors
        return sum(mobius(k) * q**(d//k) for k in divisors(d)) // d
    def smul(a, b):
        r = [0]*(N+1)
        for i, x in enumerate(a):
            if x:
                for j in range(0, N + 1 - i): r[i+j] += x*b[j]
        return r
    def spow(a, e):
        r = [1] + [0]*N
        while e:
            if e & 1: r = smul(r, a)
            a = smul(a, a); e >>= 1
        return r
    Z = [1] + [0]*N
    for d in range(1, N + 1):
        Q = q**d; G = [0]*(N+1); G[0] = 1
        for e in range(1, N//d + 1): G[e*d] = (-Q**(e-1)) if e % 2 else Q**(e-1)*(Q-2)
        Z = smul(Z, spow(G, Nirr(d)))
    print(f"   A(n) from the Euler product agrees with brute force for n <= {BD}: {all(Z[d] == Ab[d] for d in Ab)}")
    # ---------- main term ----------
    mpmath.mp.dps = 50
    def W(u):
        lw = mpmath.mpf(0)
        for d in range(1, 400):
            v = u**d; Q = mpmath.mpf(q)**d
            Wd = (1 - v - 2*Q*v*v) / ((1 - v)*(1 - Q*v*v)**2)
            lw += Nirr(d) * mpmath.log(Wd)
            if abs(Nirr(d) * (Wd - 1)) < mpmath.mpf(10)**-45: break
        return mpmath.exp(lw)
    r = mpmath.mpf(q)**mpmath.mpf(-1.5)
    F = lambda u: (1 - q*u)*(1 - q*q*u*u)**2 * W(u)
    Fp, Fm = F(r), F(-r)
    print(f"   main term: A(n) ~ 1/2[F(+q^-3/2) q^(3n/2) + F(-q^-3/2) (-q^(3/2))^n],  F(+) = {mpmath.nstr(Fp,10)}, F(-) = {mpmath.nstr(Fm,10)}")
    print("    n        A(n)/q^(3n/2)     main/q^(3n/2)    remainder/q^(3n/4)")
    for n in range(1, N + 1):
        M = (Fp * mpmath.mpf(q)**(1.5*n) + Fm * (-1)**n * mpmath.mpf(q)**(1.5*n)) / 2
        R = Z[n] - M
        if n <= 6 or n % 3 == 0:
            print(f"   {n:>3}  {mpmath.nstr(Z[n]/mpmath.mpf(q)**(1.5*n),10):>16}  {mpmath.nstr(M/mpmath.mpf(q)**(1.5*n),10):>16}  {mpmath.nstr(R/mpmath.mpf(q)**(0.75*n),8):>16}")
