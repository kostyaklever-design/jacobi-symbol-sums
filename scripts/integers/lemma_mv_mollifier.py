"""Lemma mv, step (ii): the short sum A_y(s) = sum_{n<=y} mu(n) n^-s approximates 1/zeta on Re s = 1.
Mean of |1/zeta - A_y|^2 over t in [1000, 2000] should decay like 1/y; also prints c_n in zeta*A_10 - 1."""
import numpy as np, mpmath as mp
from sympy import mobius, divisors
U = 1000; ts = np.linspace(U, 2 * U, 3001)
iz = np.array([complex(1 / mp.zeta(mp.mpc(1, t))) for t in ts])
print("mean |1/zeta|^2 =", round(np.mean(abs(iz)**2), 3))
for y in (10, 100, 1000):
    k = np.arange(1, y + 1); mu = np.array([int(mobius(x)) for x in k])
    A = np.array([np.sum(mu * k**(-1.0 - 1j * t)) for t in ts])
    e = np.mean(abs(iz - A)**2)
    print(f"y={y:5d}: mean |1/zeta - A_y|^2 = {e:.5f},  times y = {e*y:.2f}")
print("c_n in zeta*A_10 - 1:", [(m, sum(int(mobius(d)) for d in divisors(m) if d <= 10)) for m in range(1, 31)])
