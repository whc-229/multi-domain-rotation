"""稳态可行性的严格判定（任意周期）：状态 = N 台机器人的多重集 {(下一相位 p, 距可用还剩 r 格)}。
每格从 r=0 的机器人中选 K 台上岗（按相位的多重集选择），上岗者变为 (p+1, g_p − 1)，其余 r ← max(r−1, 0)。
有限图中存在无穷路径 ⇔ 存在环 ⇔ 反复删除"无存活后继"的状态后仍非空。这给出不依赖周期上限的精确判定。"""
import itertools, sys, json
from collections import Counter
def feasible(g, K, N):
    L = len(g)
    types = [(p, r) for p in range(L) for r in range(g[(p - 1) % L])]   # 刚做完相位 p−1，剩余 r < g_{p−1}
    idx = {t: i for i, t in enumerate(types)}
    states = list(itertools.combinations_with_replacement(range(len(types)), N))
    sid = {s: i for i, s in enumerate(states)}
    succ = []
    for s in states:
        rob = [types[i] for i in s]
        avail = Counter(p for p, r in rob if r == 0)
        others = [(p, max(r - 1, 0)) for p, r in rob if r > 0]
        outs = set()
        phases = sorted(avail)
        # 枚举从可用机器人中选 K 台（按相位计数）
        def rec(k, i, chosen):
            if k == 0:
                new = list(others)
                for p in phases:
                    stay = avail[p] - chosen.get(p, 0)
                    new += [(p, 0)] * stay
                    new += [((p + 1) % L, g[p] - 1)] * chosen.get(p, 0)
                outs.add(sid[tuple(sorted(idx[t] for t in new))])
                return
            if i == len(phases): return
            p = phases[i]
            for c in range(min(avail[p], k), -1, -1):
                if c: chosen[p] = c
                elif p in chosen: del chosen[p]
                rec(k - c, i + 1, chosen)
            chosen.pop(p, None)
        if sum(avail.values()) >= K: rec(K, 0, {})
        succ.append(list(outs))
    alive = [True] * len(states); changed = True
    while changed:
        changed = False
        for i in range(len(states)):
            if alive[i] and not any(alive[j] for j in succ[i]):
                alive[i] = False; changed = True
    return any(alive), len(states)
if __name__ == "__main__":
    tests = [([3, 10], 1, 6), ([3, 10], 1, 7), ([4, 6, 5, 4], 1, 5), ([4, 6, 5, 4], 1, 6), ([5, 2, 2, 11, 3], 1, 5), ([5, 2, 2, 11, 3], 1, 6),
             ([3, 3, 15], 1, 7), ([3, 3, 15], 1, 8), ([2, 4], 1, 3), ([2, 4], 1, 4), ([3, 7, 5, 4], 1, 5), ([3, 7, 5, 4], 1, 6)]
    for g, K, N in tests:
        ok, ns = feasible(g, K, N)
        print(json.dumps(dict(g=g, K=K, N=N, feasible=ok, states=ns)), flush=True)
