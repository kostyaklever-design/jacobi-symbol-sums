"""Theorem B (under RH): R_2(X) = X^{3/4} sum_rho X^{i gamma/2}(B_rho log X + A_rho) + (integral on Re s = 0.7).
1) Lambda = 2 sum_{gamma>0} |B_rho|^2 from the first NZ zeros (paper: 7.507e-6).
2) The explicit formula is checked exactly: R_2 from a sieve to 1e7 is compared with
   (sum over NZ zeros) + (1/2pi i) int_{Re s=0.7} D(s) kappa(s) X^s ds, the line integral computed numerically up to |t| = T.
Usage: python3 thmB_explicit_formula.py [NZ=100] [T=100] [h=0.01]   (full run ~15 min; quick test: 20 20 0.01).
Keep h <= 0.01: the integrand has a sharp peak near t = 1 (|g| ~ 0.55) that a coarser step misses."""
import sys, numpy as np, mpmath as mp
from common import sieve, C
NZ = int(sys.argv[1]) if len(sys.argv) > 1 else 100
T = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0
h = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01
mp.mp.dps = 15
zo = lambda w: (1 - mp.power(2, -w)) * mp.zeta(w)
kappa = lambda s: 2 / (s * (s + 1) * (s + 2))
P = [p for p in range(3, 3000) if all(p % d for d in range(2, int(p**.5) + 1))]
def H(s):
    t = mp.mpc(0)
    for p in P:
        u = mp.power(p, -s); Wp = (1 - u - 2*p*u**2) / ((1 - u)*(1 - p*u**2)**2)
        t += mp.log(Wp / ((1 - p*u**3)**2 * (1 - p**2*u**4)))
    return mp.exp(t)
W = lambda s: H(s) / (zo(3*s - 1)**2 * zo(4*s - 2))
Q = lambda s: zo(2*s - 2) * W(s) * kappa(s) / zo(s)
zs = []
for k in range(1, NZ + 1):
    rho = mp.zetazero(k); s0 = (1 + rho)/2
    c1 = 2*mp.diff(zo, rho); c2 = 4*mp.diff(zo, rho, 2); q = Q(s0); qp = mp.diff(Q, s0)
    zs.append((float(rho.imag), complex(q/c1**2), complex((qp - (c2/c1)*q)/c1**2)))
print(f"Lambda from {NZ} zeros = {2*sum(abs(B)**2 for _, B, _ in zs):.4e}   (paper: 7.507e-6)")
ts = np.arange(0, T + h/2, h)
g = np.array([complex(zo(2*s - 2)*W(s)*kappa(s)/(zo(s)*zo(2*s - 1)**2)) for s in (mp.mpc(0.7, t) for t in ts)])
N = 10**7; S, _ = sieve(N); n = np.arange(N + 1, dtype=float)
c0 = S.cumsum(); c1s = (n*S).cumsum(); c2s = (n*n*S).cumsum()
R2 = lambda X: c0[int(X)] - 2*c1s[int(X)]/X + c2s[int(X)]/X**2 - 4*6*C/105*X**1.5
zsum = lambda X: X**0.75*2*sum((np.exp(1j*g_/2*np.log(X))*(B*np.log(X) + A)).real for g_, B, A in zs)
def integral(X):
    w = np.full(len(ts), h); w[0] = w[-1] = h/2
    return X**0.7/np.pi*np.sum(w*g*np.exp(1j*ts*np.log(X))).real
print("      X      R2/norm   zeros/norm  integral/norm  zeros+int")
for X in np.geomspace(1e5, 1e7, 5):
    nrm = X**0.75*np.log(X)
    print(f"{X:11.0f}  {R2(X)/nrm:+.5f}   {zsum(X)/nrm:+.5f}     {integral(X)/nrm:+.5f}     {(zsum(X)+integral(X))/nrm:+.5f}")
Xs = np.geomspace(1e5, 1e7, 300); a = np.array([R2(X) for X in Xs]); b = np.array([zsum(X) + integral(X) for X in Xs])
print("relative rms of R2 - (zeros + integral):", round(float(np.sqrt(np.mean((a-b)**2))/np.sqrt(np.mean(a**2))), 4))
