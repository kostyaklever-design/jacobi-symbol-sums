"""Lemma mv, step (iii): replacing each factor by its short sum changes the mean square of the product little.
Product 1/zeta(s) * 1/zeta(2s-1) on Re s = 1, t in [1000, 2000]; also checks the telescoping identity."""
import numpy as np, mpmath as mp
from sympy import mobius
U = 1000; ts = np.linspace(U, 2 * U, 3001)
f2 = np.array([complex(1 / mp.zeta(mp.mpc(1, t))) for t in ts])
f3 = np.array([complex(1 / mp.zeta(mp.mpc(1, 2 * t))) for t in ts])
P = f2 * f3; MP = np.mean(abs(P)**2); print(f"true mean square = {MP:.4f}")
for y in (10, 100, 1000):
    k = np.arange(1, y + 1); mu = np.array([int(mobius(x)) for x in k], float)
    g2 = np.array([np.sum(mu * k**(-1.0 - 1j * t)) for t in ts]); g3 = np.array([np.sum(mu * k**(-1.0 - 2j * t)) for t in ts])
    Q = g2 * g3; MQ = np.mean(abs(Q)**2)
    tel = np.max(abs((P - Q) - ((f2 - g2) * f3 + g2 * (f3 - g3))))
    print(f"y={y:5d}: approx = {MQ:.4f}, rel. diff = {100*abs(MP-MQ)/MP:.1f}%, identity error = {tel:.1e}")
