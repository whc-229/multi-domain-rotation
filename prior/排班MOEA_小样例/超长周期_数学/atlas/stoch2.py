"""T5（改进）：随机大修时长下的最少机队 N_ε(K)。两种派工策略：
FIFO（就绪最久优先）；PB（相位平衡：每格按各相位累计派工缺额分配 K 个名额，同相位内 FIFO）。
确定性情形下 PB 应接近 ⌈Kā⌉；随机情形看超额 N_ε − Kā 是否 ∝ √K。"""
import numpy as np, math, sys, json, heapq, collections

def simulate(g, G0, spread, N, K, T, seed, policy, warm=3000):
    rng = np.random.default_rng(seed); L = len(g) + 1
    q = [collections.deque() for _ in range(L)]      # 各相位就绪队列（存就绪时刻，FIFO）
    for i in range(N): q[0].append((0, i))
    transit = []; D = np.zeros(L); short = 0; cnt = 0
    for t in range(T):
        while transit and transit[0][0] <= t:
            rt, i, ph = heapq.heappop(transit); q[ph].append((rt, i))
        avail = [len(x) for x in q]; tot = sum(avail)
        if tot < K:
            if t >= warm: short += 1
            take = avail[:]
        elif policy == "FIFO":
            take = [0] * L; heads = []
            # 合并各相位队头按就绪时刻取 K 个
            idx = [0] * L
            for _ in range(K):
                best = None
                for ph in range(L):
                    if idx[ph] < avail[ph]:
                        r = q[ph][idx[ph]][0]
                        if best is None or r < best[0]: best = (r, ph)
                idx[best[1]] += 1
            take = idx
        else:  # PB
            target = (t + 1) * K / L - D   # 各相位累计派工缺额
            take = [0] * L; need = K
            order = sorted(range(L), key=lambda ph: -target[ph])
            # 先按缺额向下取整分配，再按缺额大小补足
            for ph in order:
                x = min(avail[ph], max(0, int(math.floor(target[ph]))), need); take[ph] += x; need -= x
            while need > 0:
                for ph in order:
                    if need and take[ph] < avail[ph]: take[ph] += 1; need -= 1
        for ph in range(L):
            for _ in range(take[ph]):
                _, i = q[ph].popleft()
                gap = g[ph] if ph < L - 1 else max(1, int(rng.integers(G0 - spread, G0 + spread + 1)))
                heapq.heappush(transit, (t + gap, i, (ph + 1) % L))
            D[ph] += take[ph]
        if t >= warm: cnt += 1
    return short / max(1, cnt)

def min_fleet(g, G0, spread, K, eps, T, policy, seed=0):
    abar = (sum(g) + G0) / (len(g) + 1); lo = math.ceil(K * abar - 1e-9); hi = lo
    while simulate(g, G0, spread, hi, K, T, seed, policy) > eps: hi += max(1, hi // 8)
    while lo < hi:
        mid = (lo + hi) // 2
        if simulate(g, G0, spread, mid, K, T, seed, policy) <= eps: hi = mid
        else: lo = mid + 1
    return lo, abar

if __name__ == "__main__":
    configs = [([3], 10, [0, 2, 4, 6]), ([4, 4], 5, [0, 1, 2]), ([2, 2], 5, [0, 1, 2])]
    out = []
    for gs, G0, spreads in configs:
        for spread in spreads:
            for policy in ["PB"]:
                for K in [1, 2, 4, 8, 16, 32, 64, 128]:
                    T = 40000 if K <= 16 else (20000 if K <= 64 else 12000)
                    N, a = min_fleet(gs, G0, spread, K, 1e-3, T, policy)
                    line = dict(g=gs + [G0], spread=spread, policy=policy, K=K, N=N, excess=round(N - K * a, 2), ratio=round((N - K * a) / math.sqrt(K), 2))
                    out.append(line); print(json.dumps(line), flush=True)
    json.dump(out, open("stoch2.json", "w"), indent=1)
