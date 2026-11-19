"""Lemma mv, step (v): local factors m_p(eps) of M_m(eps); m_p(0) = (p^2+1)/(p^2-1) and (m_p - 1) p^(3/2) stays bounded."""
import numpy as np
def m_p(p, eps, chi, J=600):
    a = p**(-1 - eps); b = (1 + chi) * p**(-1 - 2 * eps); c = p**(-1 + 2 * eps)
    e = np.zeros(J); num = [1, -a, -b]
    for j in range(J): e[j] = (num[j] if j < 3 else 0) + (c * e[j - 2] if j >= 2 else 0)
    return np.sum(e**2)
for p in (3, 5, 7, 101, 1009):
    for chi in (1, -1):
        v = [m_p(p, e, chi) for e in (0, 1 / 16, 1 / 8)]
        print(f"p={p:4d} chi={chi:+d}: m_p(0)={v[0]:.6f} vs {(p*p+1)/(p*p-1):.6f};  (m_p-1)p^1.5 = " + ", ".join(f"{(x-1)*p**1.5:.3f}" for x in v))
