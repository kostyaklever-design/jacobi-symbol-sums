"""Shared helpers: S(n), g(n) for odd n by a multiplicative sieve.
S(p^e) = -p^(e-1) (e odd), p^(e-1)(p-2) (e even);  g(p) = -1, g(p^2) = -2p, g(p^k) = 0 (k >= 3)."""
import numpy as np
C = 0.0586383825  # main-term constant C = G(3/2)/6 (m = 1)

def sieve(N):
    """Return arrays S, g on 0..N (only odd indices are meaningful)."""
    spf = np.zeros(N + 1, dtype=np.int32)
    for p in range(3, int(N ** 0.5) + 1, 2):
        if spf[p] == 0:
            blk = spf[p * p::2 * p]; blk[blk == 0] = p
    S = np.zeros(N + 1); g = np.zeros(N + 1); S[1] = g[1] = 1
    for n in range(3, N + 1, 2):
        p = spf[n]
        if p == 0:
            S[n] = g[n] = -1; continue
        m, e = n, 0
        while m % p == 0:
            m //= p; e += 1
        pk = p ** (e - 1)
        S[n] = (-pk if e % 2 else pk * (p - 2)) * S[m]
        g[n] = (-1 if e == 1 else (-2 * p if e == 2 else 0)) * g[m]
    return S, g
