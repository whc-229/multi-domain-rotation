"""对剩余均衡拉伸无法给出 c≤1 的余数模式，抽具体序列直接算 c。"""
import itertools, math, json, random
from fractions import Fraction as F
from net import solve
from balanced_stretch import emin
def cval(g, K, Pmax=None):
    S = sum(g); L = len(g); lb = math.ceil(F(K * S, L)); q = L // math.gcd(K, L); best = None
    for P in range(q, (Pmax or 2 * S) + 1, q):
        v, x = solve(g, K, P, time_limit=30)
        if v is None: continue
        N = int(round(v)); best = N if best is None else min(best, N)
        if best == lb: break
    return best - lb
random.seed(5)
for L, j in [(4, 2), (5, 1), (5, 2), (5, 3)]:
    pats = []
    seen = set()
    for res in itertools.product(range(L), repeat=L):
        key = min(tuple(res[i:] + res[:i]) for i in range(L))
        if key in seen: continue
        seen.add(key)
        S = sum(res); W = (-j * S) % L
        if (emin(res, j, L) - W) // L >= 2: pats.append(res)
    for res in pats:
        for _ in range(3):
            g = [r + L * random.randint(0 if r else 1, 2) for r in res]
            print(json.dumps(dict(L=L, K=j, residues=res, g=g, c=cval(g, j))), flush=True)
