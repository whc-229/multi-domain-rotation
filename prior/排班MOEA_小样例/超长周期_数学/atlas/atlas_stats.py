import json, math
from fractions import Fraction as F
from collections import Counter
rows = []
for it in json.load(open("atlas_random_80.json")):
    g = it["g"]; S = sum(g); L = len(g)
    for r in it["rows"]:
        K = r["K"]; lb = math.ceil(F(K * S, L)); W = lb * L - K * S
        rows.append(dict(g=g, K=K, gap=r["N_wait"] - lb, fixgap=(r["N_fix"] - lb) if r.get("N_fix") is not None else None,
                         W=W, zero=(W == 0), tiling=r.get("tiling")))
d = json.load(open("net_wait_K12.json"))
for key, rs in d.items():
    g = json.loads(key); S = sum(g); L = len(g)
    for r in rs:
        K = r["K"]; lb = math.ceil(F(K * S, L)); W = lb * L - K * S
        rows.append(dict(g=g, K=K, gap=r["N"] - lb, fixgap=None, W=W, zero=(W == 0), tiling=None))
print("(g,K) 组合数", len(rows), "序列数", len(set(tuple(r["g"]) for r in rows)))
print("稳态等待 gap 分布", Counter(r["gap"] for r in rows))
print("零余量：gap 分布", Counter(r["gap"] for r in rows if r["zero"]))
print("零余量 gap 与铺砌判定一致：", sum(1 for r in rows if r["zero"] and r["tiling"] is not None and ((r["gap"] == 0) == r["tiling"])), "/", sum(1 for r in rows if r["zero"] and r["tiling"] is not None))
pos = [r for r in rows if not r["zero"]]
print("正余量：按 W 分 gap", {w: dict(Counter(r["gap"] for r in pos if r["W"] == w)) for w in sorted(set(r["W"] for r in pos))})
fx = [r for r in rows if r["fixgap"] is not None]
print("固定模式 − 等待 的台数差分布", Counter(r["fixgap"] - r["gap"] for r in fx))
