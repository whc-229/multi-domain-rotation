"""一致性检查：构造性上界（引理 6′ 与 6″ 取小）必须 ≥ 网络 ILP 得到的 c。"""
import json, math
from fractions import Fraction as F
from balanced_stretch import emin
from mix_stretch import mix_bound
rows = []
for it in json.load(open("atlas_random_80.json")):
    for r in it["rows"]: rows.append((it["g"], r["K"], r["N_wait"]))
for key, rs in json.load(open("net_wait_K12.json")).items():
    for r in rs: rows.append((json.loads(key), r["K"], r["N"]))
viol = 0; tight = 0; n = 0; cache = {}
for g, K, N in rows:
    L = len(g); S = sum(g); lb = math.ceil(F(K * S, L)); c = N - lb
    k = K % L
    if k == 0: b = 0
    else:
        res = tuple(x % L for x in g); key = (min(tuple(res[i:] + res[:i]) for i in range(L)), k)
        if key not in cache:
            W = (-k * sum(key[0])) % L
            b1 = (emin(key[0], k, L) - W) // L; b2 = mix_bound(key[0], k, L) if k >= 2 else None
            cache[key] = min(x for x in [b1, b2] if x is not None)
        b = cache[key]
    n += 1; viol += (c > b); tight += (c == b)
    if c > b: print("违反！", g, K, N, c, b)
print("检查", n, "组；违反", viol, "；上界取等", tight)
