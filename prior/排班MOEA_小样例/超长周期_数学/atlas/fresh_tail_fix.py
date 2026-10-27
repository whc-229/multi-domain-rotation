import math, json, sys, time
from fresh_tail import fresh_tail
cases = [([3,10],3,20),([2,9],5,28),([3,6,3,6,9],5,30),([3,10],9,60),([2,9],10,56),([3,6,3,6,9],10,60)]
for g, K, target in cases:
    S = sum(g); L = len(g); q = L // math.gcd(K, L); best = None; t0 = time.time()
    for T in [S, 2*S, 3*S, 4*S]:
        for P in range(q, 3*S + 1, q):
            if P < 2: continue
            v = fresh_tail(g, K, T, P, time_limit=60)
            if v is not None and (best is None or v < best[0]): best = (v, T, P)
            if best and best[0] <= target: break
        if best and best[0] <= target: break
    print(json.dumps(dict(g=g, K=K, lower=target, upper=best, sec=round(time.time()-t0))), flush=True)
