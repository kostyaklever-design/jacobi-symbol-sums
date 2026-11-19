"""Lemma mv, step (iv), off-diagonal: solutions of l k^2 - l' k'^2 = Delta (k, k' <= 2000) are few:
count <= tau(Delta) * log U, and they form chains under the unit of Q(sqrt(l l'))."""
import math
from sympy import divisor_count as tau
K = 2000; logU = math.log(1000)
def sols(l, lp, D):
    out = []
    for kp in range(1, K + 1):
        v = D + lp * kp * kp
        if v <= 0 or v % l: continue
        k = math.isqrt(v // l)
        if k <= K and l * k * k == v: out.append((k, kp))
    return out
print("l=1, l'=2, Delta=1:", sols(1, 2, 1), " unit^2 =", round((1 + 2**.5)**2, 3))
print("l=1, l'=2, Delta=7:", sols(1, 2, 7))
for l, lp in ((1, 2), (1, 1)):
    cnt = {D: len(sols(l, lp, D)) for D in range(1, 1001)}
    ratio = max(c / (tau(D) * logU) for D, c in cnt.items() if c)
    print(f"l={l}, l'={lp}: max count/(tau*logU) = {ratio:.2f}, sum count/Delta = {sum(c/D for D,c in cnt.items()):.2f}, (log U)^2 = {logU**2:.1f}")
