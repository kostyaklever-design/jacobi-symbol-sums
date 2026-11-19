"""Theorem D, Lemma D0: exact identity sum_{n<=x odd} S(n) = sum_d g(d)[Y^3/6 + Y^2 s + 2Y B2(s) + 4/3 B3(s)],
Y = sqrt(x/d), s = saw((Y+1)/2); also checks |g(a b^2)| = 2^omega(b) b and g = 0 on cube-divisible d."""
import math, random
from sympy import factorint
from common import sieve
X = 200000; S, g = sieve(X)
bad = 0
for d in range(1, X + 1, 2):
    f = factorint(d)
    if any(e > 2 for e in f.values()): bad += g[d] != 0; continue
    b = math.prod(p for p, e in f.items() if e == 2); om = sum(1 for e in f.values() if e == 2)
    bad += abs(g[d]) != 2**om * b
print("violations of |g(ab^2)| = 2^omega(b) b:", bad)
saw = lambda u: 0.0 if u == int(u) else 0.5 - (u - math.floor(u))
random.seed(1)
for x in [random.uniform(1000, X) for _ in range(4)]:
    lhs = sum(S[n] for n in range(1, int(x) + 1, 2)); rhs = 0.0
    for d in range(1, int(x) + 1, 2):
        if g[d] == 0: continue
        Y = math.sqrt(x / d); s = saw((Y + 1) / 2)
        rhs += g[d] * (Y**3 / 6 + Y * Y * s + 2 * Y * (s * s - 1 / 12) + 4 / 3 * (s**3 - s / 4))
    print(f"x={x:10.3f}: left {int(lhs)}, right {rhs:.4f}")
