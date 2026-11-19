"""Theorem C: the coefficient B_rho = Q(s_rho)/c_1^2 of the leading (double-pole) wave for the first zero of zeta,
s_rho = (1+rho)/2, Q(s) = zeta_o(2s-2) W(s) kappa(s) / zeta_o(s), c_1 = 2 zeta_o'(rho), kappa(s) = 2/(s(s+1)(s+2)).
On Re s = 3/4 W is taken as H / [zeta_o(3s-1)^2 zeta_o(4s-2)], H an absolutely convergent Euler product.
Also prints |omega_0|, the loss exponent 1 + log(4/3)/log|omega_0|, the gap Delta_* to the second zero,
and the real zeros of the local factors 1 - u - 2p u^2 (only p = 3, 5 lie in Re s >= 3/4)."""
import mpmath as mp, numpy as np
mp.mp.dps = 30
zo = lambda w: (1 - mp.power(2, -w)) * mp.zeta(w)
kappa = lambda s: 2 / (s * (s + 1) * (s + 2))
def odd_primes(N):
    sv = np.ones(N + 1, bool); sv[:2] = False
    for i in range(2, int(N**.5) + 1):
        if sv[i]: sv[i*i::i] = False
    return [int(p) for p in np.nonzero(sv)[0] if p > 2]
def H(s, N):
    tot = mp.mpc(0)
    for p in odd_primes(N):
        u = mp.power(p, -s)
        Wp = (1 - u - 2*p*u**2) / ((1 - u) * (1 - p*u**2)**2)
        tot += mp.log(Wp / ((1 - p*u**3)**2 * (1 - p**2*u**4)))
    return mp.exp(tot)
rho = mp.zetazero(1); s0 = (1 + rho) / 2; c1 = 2 * mp.diff(zo, rho)
for N in (10**4, 10**5):
    W = H(s0, N) / (zo(3*s0 - 1)**2 * zo(4*s0 - 2))
    B = zo(2*s0 - 2) * W * kappa(s0) / zo(s0) / c1**2
    print(f"primes <= 1e{int(np.log10(N))}: |B_rho0| = {mp.nstr(abs(B), 6)}")
w0 = abs(s0); w1 = abs((1 + mp.zetazero(2)) / 2)
print("|omega_0| =", mp.nstr(w0, 8), "  loss exponent =", mp.nstr(1 + mp.log(4/3) / mp.log(w0), 6), "  Delta_* =", mp.nstr(mp.log(w1 / w0), 6))
for p in (3, 5, 7, 11):
    u = (-1 + mp.sqrt(1 + 8*p)) / (4*p); print(f"p={p}: zero of 1-u-2pu^2 at Re s = {mp.nstr(-mp.log(u)/mp.log(p), 6)}")
