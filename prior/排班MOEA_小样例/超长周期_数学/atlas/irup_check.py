"""检验猜想：全新启动最少台数 N_fresh(K) = max(⌈α·K⌉, N_steady(K))，α 为全新启动流体 LP 值（H=3S 的精确对偶下界）。
下界：⌈αK⌉（对偶证书）与 N_steady（任何全新无穷排班都含一个稳态周期段）；上界：fresh_tail 构造的可无限延拓排班。"""
import math, json, sys, random
from fractions import Fraction as F
from cert import dual_cert
from fresh_tail import fresh_tail
from net import solve
def steady(g, K, reset=False):
    S = sum(g); L = len(g); lb = math.ceil(K * F(S, L)); q = L // math.gcd(K, L); best = None
    for P in range(q, max(2 * S, 24) + 1, q):
        v, x = solve(g, K, P, reset=reset, time_limit=30)
        if v is None: continue
        N = int(round(v)); best = N if best is None else min(best, N)
        if best <= lb: break
    return best
random.seed(7)
seqs = [[1,4],[1,6],[2,4],[2,3],[1,1,5],[2,2,7],[3,4],[2,7],[1,2,6],[3,3,4,8]]
while len(seqs) < 30:
    L = random.choice([2, 3, 3, 4]); g = [random.randint(1, 7) for _ in range(L - 1)] + [random.randint(4, 14)]
    if g not in seqs: seqs.append(g)
reset = len(sys.argv) > 1 and sys.argv[1] == "reset"
Kmax = 6; tot = 0; ok = 0
for g in seqs:
    S = sum(g); L = len(g)
    alpha, lpval, _, _ = dual_cert(g, 3 * S, reset)
    rig = alpha is not None
    if not rig: alpha = F(lpval).limit_denominator(1000) - F(1, 10**6)   # 对偶有理化失败：近似值，非严格
    row = []
    for K in range(1, Kmax + 1):
        st = steady(g, K, reset); lower = max(math.ceil(alpha * K), st if st else 0)
        q = L // math.gcd(K, L); up = None
        for T in [S, 2 * S, 3 * S]:
            for P in range(max(q, 2), 2 * S + 1, q):
                v = fresh_tail(g, K, T, P, reset, time_limit=30)
                if v is not None and (up is None or v < up): up = v
                if up is not None and up <= lower: break
            if up is not None and up <= lower: break
        tot += 1; ok += (up == lower); row.append((K, lower, up))
    print(json.dumps(dict(g=g, reset=reset, rigorous=rig, alpha=str(alpha), abar=str(F(S, L)), rows=row)), flush=True)
print("相等", ok, "/", tot)
