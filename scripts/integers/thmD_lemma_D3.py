"""Theorem D, Lemma D3: the average over n of the square terms equals A(q, alpha); sign rule s(q/(4e)) = mu(e)/4;
growth prod_{p<=z, p=3 mod 4} (1+2/p) ~ 0.88 log z."""
import numpy as np, math
from sympy import factorint, divisors, primerange
def saw(u):
    u = np.asarray(u, float); f = u - np.floor(u); return np.where(f == 0, 0.0, 0.5 - f)
def mu_om(n):
    f = factorint(n)
    return (0, 0) if any(e > 1 for e in f.values()) else ((-1)**len(f), len(f))
q = 21; N = 10**6; B = 2000; n = np.arange(N // 2 + 1, N + 1, dtype=float)
far = np.prod([1 - 2 / p**2 for p in primerange(3, 10**6) if q % p])
for alpha in (q / 4, 3 * q / 4):
    v = q * n + alpha; tot = 0.0
    for b in range(1, B + 1, 2):
        m, om = mu_om(b)
        if m: tot += m * 2**om / b * np.mean(saw((v + b) / (2 * b)))
    A = 0.5 * far * sum(mu_om(e)[0] * 2**mu_om(e)[1] / e * float(saw(alpha / e)) for e in divisors(q))
    print(f"alpha={alpha:6.2f}: direct average = {tot:+.6f},  A(q,alpha) = {A:+.6f}")
pr = 1.0; ps = [p for p in primerange(3, 10**7) if p % 4 == 3]; i = 0
for z in (10**2, 10**4, 10**6, 10**7):
    while i < len(ps) and ps[i] <= z: pr *= 1 + 2 / ps[i]; i += 1
    print(f"z=1e{int(math.log10(z))}: product = {pr:.3f}, / log z = {pr/math.log(z):.4f}")
