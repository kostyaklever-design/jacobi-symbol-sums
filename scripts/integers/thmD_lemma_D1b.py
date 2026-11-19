"""Theorem D, Lemma D1b: E(x)/x is close to the truncated sawtooth sum sum_{d<=y} g(d)/d * s((sqrt(x/d)+1)/2)."""
import numpy as np, math, random
from common import sieve, C
N = 10**7; S, g = sieve(N); cs = np.cumsum(S)
d = np.arange(1, N + 1, 2, dtype=float); gd = g[1::2]
def trunc(x, y):
    dd = d[d <= y]; gg = gd[:len(dd)]; u = (np.sqrt(x / dd) + 1) / 2
    return np.sum(gg / dd * (0.5 - (u - np.floor(u))))
random.seed(7)
for x in sorted(random.uniform(1e5, N) for _ in range(5)):
    E = (cs[int(x)] - C * x**1.5) / x; y3 = x * math.exp(-math.sqrt(math.log(x)))
    print(f"x={x:11.1f}: E/x={E:+.4f}  y=x/10: {trunc(x,x/10):+.4f}  y=x/100: {trunc(x,x/100):+.4f}  y=x*exp(-sqrt(log x)): {trunc(x,y3):+.4f}")
