"""Theorem A (m = 1): logarithmic mean square of E(x)/x tends to 1/pi^2.
E(x) = sum_{n<=x, n odd} S(n) - C x^(3/2) is evaluated exactly on each [n, n+2);
the integral of (E/x)^2 dx/x uses composite Simpson (64 panels per segment below 1e4, 2 above)."""
import numpy as np
from common import sieve, C
N = 10**7
S, _ = sieve(N)
A = np.cumsum(S[1::2]); n = np.arange(1, N + 1, 2, dtype=float)
def simpson(a, Aval, k):
    h = 2.0 / k; tot = 0
    for i in range(k + 1):
        x = a + i * h; w = 1 if i in (0, k) else (4 if i % 2 else 2)
        tot = tot + w * ((Aval - C * x**1.5) / x)**2 / x
    return tot * h / 3
small = n < 1e4
I = np.empty_like(n); I[small] = simpson(n[small], A[small], 64); I[~small] = simpson(n[~small], A[~small], 2)
cum = np.cumsum(I)
for Y in (10**5, 10**6, 10**7):
    j = (Y - 1) // 2; tot = cum[j - 1]; part = tot - cum[(10**4 - 1) // 2 - 1]
    print(f"Y=1e{int(np.log10(Y))}: int_1^Y/logY = {tot/np.log(Y):.4f}   int_1e4^Y/log(Y/1e4) = {part/np.log(Y/1e4):.4f}")
print("1/pi^2 =", round(1 / np.pi**2, 6))
