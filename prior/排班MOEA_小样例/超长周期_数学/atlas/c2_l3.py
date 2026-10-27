import math, json
from fractions import Fraction as F
from net import solve
def cval(g, K):
    S = sum(g); L = len(g); lb = math.ceil(F(K * S, L)); q = L // math.gcd(K, L); best = None
    for P in range(q, 2 * S + 1, q):
        v, x = solve(g, K, P, time_limit=30)
        if v is None: continue
        N = int(round(v)); best = N if best is None else min(best, N)
        if best == lb: break
    return best - lb, lb * L - K * S
for g in [[3,3,2],[3,6,5],[6,3,8],[3,3,5],[6,6,2],[3,9,2],[9,3,5],[3,3,11]]:
    print(g, [(K,) + cval(g, K) for K in [2, 5]], flush=True)
