"""针对猜想 C1（c≤1）与 C2（正余量 c=1 ⇒ W=1）的定向搜索：
随机抽 L=3,4,5、间隔 1..12 的序列，K=1、2，记录 W 与 c（周期 P ≤ 2S 的网络 ILP）。"""
import random, math, json
from fractions import Fraction as F
from net import solve
random.seed(11)
out = []; seen = set()
while len(out) < 160:
    L = random.choice([3, 4, 5]); g = [random.randint(1, 12) for _ in range(L)]
    K = random.choice([1, 2]); S = sum(g)
    key = (tuple(g), K)
    if key in seen: continue
    seen.add(key)
    lb = math.ceil(F(K * S, L)); W = lb * L - K * S
    if W == 0: continue                     # 只看正余量
    q = L // math.gcd(K, L); best = None
    for P in range(q, 2 * S + 1, q):
        v, x = solve(g, K, P, time_limit=30)
        if v is None: continue
        N = int(round(v)); best = N if best is None else min(best, N)
        if best == lb: break
    c = best - lb if best is not None else None
    out.append(dict(g=g, K=K, W=W, c=c)); print(json.dumps(out[-1]), flush=True)
from collections import Counter
print("按 W 统计 c：", {w: dict(Counter(r["c"] for r in out if r["W"] == w)) for w in sorted(set(r["W"] for r in out))})
