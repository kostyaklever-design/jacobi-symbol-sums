"""Theorem D_m: g_m(p^2) = -(1+chi_m(p)) p (or -p for p | m) by brute force; growth (log z)^(1/2) of the product over
primes p = -1 mod k0 with chi_m(p) = 1; exact sign rule s(q/(k0 e)) = mu(e) c_k0 in rational arithmetic."""
import math
from fractions import Fraction as F
from sympy import primerange, jacobi_symbol as J, divisors, factorint
def Sm(n, m): return sum(J((x * x - m) % n, n) for x in range(n))
chi = lambda m, p: J(m % p, p)
bad = sum((Sm(p * p, m) - p * p) != (-p if m % p == 0 else -(1 + chi(m, p)) * p) for m in (2, 3, 5, -1) for p in primerange(3, 40))
print("violations of g_m(p^2) formula:", bad)
for m in (5, 2, 3, -1):
    k0 = 3 if m == -1 else 4; pr = 1.0
    for p in primerange(5 if k0 == 3 else 3, 10**7):
        if m % p and p % k0 == k0 - 1 and chi(m, p) == 1: pr *= 1 + 2 / p
    print(f"m={m:2d}, k0={k0}: product/(log 1e7)^0.5 = {pr/math.log(1e7)**0.5:.3f}")
def saw(u):
    f = u - math.floor(u); return F(0) if f == 0 else F(1, 2) - f
for m, k0, c in ((5, 4, F(1, 4)), (-1, 3, F(1, 6))):
    P = [p for p in primerange(3, 200) if m % p and p % k0 == k0 - 1 and chi(m, p) == 1][:4]; q = math.prod(P)
    print(f"m={m}, q={q}: sign rule holds for all divisors: {all(saw(F(q, k0*e)) == (-1)**len(factorint(e))*c for e in divisors(q))}")
