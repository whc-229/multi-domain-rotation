"""T5：随机维修时长。间隔序列中小间隔确定，大修间隔 G 随机（均值 G0，取整分布）。
策略：每格从就绪机器人中选 K 台，按"就绪最久优先"（FIFO）。
服务水平：长时域内出现缺岗（就绪 < K）的格子比例 ≤ ε。二分最少台数 N_ε(K)，看 N_ε − K·ā 随 K 的增长。"""
import numpy as np, math, sys, json, heapq

def simulate(gsmall, G0, spread, N, K, T, rng, warm=2000):
    L = len(gsmall) + 1
    # 状态：每台的就绪时刻与相位；用堆维护在途机器人
    ready = []   # (ready_since, id)
    transit = []  # (return_time, id, phase)
    phase = [0] * N
    for i in range(N): ready.append((0, i))
    short = 0; ready.sort()
    import collections
    q = collections.deque(sorted(ready))
    for t in range(T):
        while transit and transit[0][0] <= t:
            rt, i, ph = heapq.heappop(transit); q.append((rt, i))
        if len(q) < K:
            if t >= warm: short += 1
            take = len(q)
        else: take = K
        for _ in range(take):
            _, i = q.popleft(); ph = phase[i]
            if ph < L - 1: gap = gsmall[ph]
            else: gap = max(1, int(round(rng.uniform(G0 - spread, G0 + spread))))
            phase[i] = (ph + 1) % L
            heapq.heappush(transit, (t + gap, i, phase[i]))
    return short / max(1, T - warm)

def min_fleet(gsmall, G0, spread, K, eps, T, seed=0):
    abar = (sum(gsmall) + G0) / (len(gsmall) + 1)
    lo = math.ceil(K * abar); hi = lo
    rng = np.random.default_rng(seed)
    while simulate(gsmall, G0, spread, hi, K, T, np.random.default_rng(seed)) > eps: hi = hi + max(1, hi // 10)
    while lo < hi:
        mid = (lo + hi) // 2
        if simulate(gsmall, G0, spread, mid, K, T, np.random.default_rng(seed)) <= eps: hi = mid
        else: lo = mid + 1
    return lo, abar

if __name__ == "__main__":
    out = []
    for (gs, G0, spread) in [([3], 10, 0), ([3], 10, 3), ([3], 10, 6), ([4, 4, 5, 4], 9, 0), ([4, 4, 5, 4], 9, 4)]:
        for K in [1, 2, 4, 8, 16, 32, 64]:
            N, a = min_fleet(gs, G0, spread, K, eps=1e-3, T=60000)
            line = dict(g=gs + [G0], spread=spread, K=K, N=N, Kabar=round(K * a, 2), excess=round(N - K * a, 2), sqrtK=round(math.sqrt(K), 2))
            out.append(line); print(json.dumps(line, ensure_ascii=False), flush=True)
    json.dump(out, open("stoch.json", "w"), ensure_ascii=False, indent=1)
