"""结构分解启发式：K 个值班点 → 把机器人分成 K 组，每组专守一个点（单点问题用精确 DFS 找环）。
外层：分组（按速率 LPT 贪心起步，失败则在组间交换机器人的局部搜索）；内层：精确、带证书。
分组限制了跨点共享，所以只给上界；但每个内层问题都是小规模、可精确求解的。"""
import time, numpy as np
from slotmodel import rate
import kpost

def solve_group(Rg, cache, node_limit):
    key = tuple(sorted(id(r) for r in Rg))
    if key not in cache:
        if sum(rate(r) for r in Rg) < 1 - 1e-12: cache[key] = (False, 0)
        else: cache[key] = kpost.exact(Rg, 1, node_limit)
    return cache[key]

def decomp(pool_N, K, rng, tries=60, node_limit=200_000):
    t = time.time(); cache = {}; nodes = 0
    groups = [[] for _ in range(K)]
    for r in sorted(pool_N, key=lambda r: -rate(r)):        # LPT：给当前速率和最小的组
        min(groups, key=lambda g: sum(rate(x) for x in g)).append(r)
    for it in range(tries):
        res = [solve_group(g, cache, node_limit) for g in groups]
        if all(ok for ok, _ in res): return True, it, round(time.time() - t, 1)
        bad = [i for i, (ok, _) in enumerate(res) if not ok]
        a = bad[rng.integers(len(bad))]; b = rng.integers(K)
        if a == b or not groups[a] or not groups[b]: continue
        i, j = rng.integers(len(groups[a])), rng.integers(len(groups[b]))
        groups[a][i], groups[b][j] = groups[b][j], groups[a][i]
    return False, tries, round(time.time() - t, 1)

if __name__ == "__main__":
    from slotmodel import hetero_pool
    from exp_fn import lbK
    for K in [1, 2, 3, 4]:
        for ps in range(3):
            pool = sorted(hetero_pool(8 * K, 100 + 10 * K + ps), key=lambda r: -rate(r))
            lb = lbK(pool, K); best = None
            for N in range(lb, min(len(pool), lb + 6) + 1):
                ok, it, secs = decomp(pool[:N], K, np.random.default_rng(ps))
                if ok: best = (N, it, secs); break
            print(f"K={K} 池{ps} LB={lb} 分解启发式最少={best}", flush=True)
