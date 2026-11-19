"""Lemma mv, step (iv), diagonal: for Z(t)R(t), Z = sum_{k<=2U} k^-1 k^-2it, R = sum_{l<=10} mu(l)/l l^-it,
the mean square over [U, 2U] is close to the diagonal sum |alpha_n|^2, which tends to the formal-series value."""
import numpy as np
from sympy import mobius
from collections import defaultdict
U = 1000; y = 10
mu = {l: int(mobius(l)) for l in range(1, 200001)}
al = defaultdict(float)
for k in range(1, 2 * U + 1):
    for l in range(1, y + 1):
        if mu[l]: al[k * k * l] += mu[l] / (k * l)
diag = sum(v * v for v in al.values())
N = 200000; a = np.zeros(N + 1)
for k in range(1, int(N**0.5) + 1):
    for l in range(1, N // (k * k) + 1):
        if mu[l]: a[k * k * l] += mu[l] / (k * l)
full = np.sum(a**2)
ts = np.linspace(U, 2 * U, 40001); kk = np.arange(1, 2 * U + 1, dtype=float); Z = np.zeros_like(ts, dtype=complex)
for i in range(0, len(ts), 2000):
    t = ts[i:i + 2000, None]; Z[i:i + 2000] = np.sum(kk**-1.0 * np.exp(-2j * t * np.log(kk)), axis=1)
l = np.arange(1, y + 1); m = np.array([mu[x] for x in l], float)
R = np.array([np.sum(m / l * np.exp(-1j * t * np.log(l))) for t in ts])
mean = np.trapezoid(abs(Z * R)**2, ts) / U
print(f"true mean square = {mean:.4f}\ndiagonal         = {diag:.4f}\noff-diagonal     = {mean-diag:+.4f}\nformal series    = {full:.4f}")
