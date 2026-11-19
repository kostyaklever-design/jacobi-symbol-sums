"""Theorem D, Lemma D2: for non-square d and t = hq/(2 sqrt d) >= 1/2, ||t|| >= 1/(6 (hq)^2).
Prints the worst ratio ||t|| * 6 (hq)^2 (must be >= 1)."""
from mpmath import mp, mpf, sqrt, nint, fabs
import math
mp.dps = 50
for q in (3, 15, 21):
    worst = (mpf(10**9), None)
    for h in range(1, 31):
        hq = h * q
        for d in range(2, hq * hq + 1):
            r = math.isqrt(d)
            if r * r == d: continue
            t = mpf(hq) / (2 * sqrt(d))
            if t < 0.5: continue
            k = nint(t); ratio = fabs(t - k) * 6 * hq * hq
            if ratio < worst[0]: worst = (ratio, (h, d, int(k)))
    h, d, k = worst[1]
    print(f"q={q}: min ratio = {float(worst[0]):.4f} at h={h}, d={d}, k={k}; (hq)^2-4k^2d = {(h*q)**2-4*k*k*d}")
