"""Theorem D, Lemma D1: L(t) = sum_{d<=t} g(d)/d oscillates and decays (heavy cancellation vs sum |g|/d)."""
import numpy as np, math
from common import sieve
N = 10**7; _, g = sieve(N)
d = np.arange(N + 1, dtype=float); d[0] = 1
L = np.cumsum(g / d); La = np.cumsum(np.abs(g) / d)
for k in range(2, 8):
    t = 10**k
    print(f"t=1e{k}: L={L[t]:+.5f}  max|L| on [t/10,t]={np.max(np.abs(L[t//10:t+1])):.5f}  sum|g|/d={La[t]:.2f}  (log t)^3={math.log(t)**3:.1f}")
