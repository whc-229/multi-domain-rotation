"""T6：分组损失。把 K 个点分成 m 组各自独立排班（组间不共享机器人）：loss = Σ N(K_i) − N(ΣK_i)。
初等界：loss ≤ Σ(⌈K_i ā⌉ + c(K_i)) − ⌈Kā⌉ ≤ (m−1) + Σ c(K_i)（因 Σ⌈x_i⌉ − ⌈Σx_i⌉ ≤ m−1）。"""
import json, math
from fractions import Fraction as F
d = json.load(open("net_wait_K12.json"))
def partitions(n, m, mn=1):
    if m == 1:
        if n >= mn: yield (n,)
        return
    for k in range(mn, n // m + 1):
        for rest in partitions(n - k, m - 1, k): yield (k,) + rest
for key, rows in d.items():
    g = json.loads(key); S = sum(g); L = len(g); a = F(S, L)
    N = {r["K"]: r["N"] for r in rows}; c = {k: N[k] - math.ceil(k * a) for k in N}
    K = max(N); line = []
    for m in [2, 3, 4, K]:
        worst = None
        for p in partitions(K, m):
            loss = sum(N[k] for k in p) - N[K]; bound = (m - 1) + sum(c[k] for k in p)
            assert loss <= bound
            if worst is None or loss > worst[0]: worst = (loss, p, bound)
        line.append(f"m={m}: 最大损失 {worst[0]}（{'+'.join(map(str, worst[1]))}，界 {worst[2]}）")
    print(g, "K=12;", "；".join(line), "；专点 K·N(1)−N(K) =", K * N[1] - N[K])
