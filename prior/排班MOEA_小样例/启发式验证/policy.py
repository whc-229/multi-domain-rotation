"""策略编码：不编码排班表，而编码"派谁"的打分规则（5 个权重）。确定性仿真 → 有限状态 → 必然进入环或断档。
特征（对每台可用机器人）：
 f0 下一班后要等的槽数 g（越小越"便宜"）
 f1 选它之后未来 W 槽内的最小供需余量（前瞻）
 f2 距离下一次长维修还剩几班（越大越不急）
 f3 已空闲槽数
 f4 该机器人长期速率（越高越值得用）"""
import numpy as np
from kpost import start, stepK
from slotmodel import canon, rate

def features(R, st, j, K, W=30):
    rb = R[j]; ph = st[j][1]; g = rb.g[ph]
    per = rb.g; big = max(rb.g[rb.mu:])
    k = 0; p = ph
    while rb.g[p] < big and k < rb.L: k += 1; p = rb.nph(p)
    return g, k, rate(rb)

def run_policy(R, K, wts, maxstep=6000, W=30):
    st = start(R); seen = {}; idle = [0] * len(R)
    for t in range(maxstep):
        c = canon(R, st)
        if c in seen: return True, t, t - seen[c]
        seen[c] = t
        av = [j for j, (w, _) in enumerate(st) if w == 0]
        if len(av) < K: return False, t, 0
        # 前瞻余量：按"每槽需 K 台"估计，选中后各机器人返回时刻的累计供给 − 需求
        ret = np.zeros(W + 2)
        for (w, ph) in st:
            if 0 < w <= W: ret[w] += 1
        base_av = len(av)
        sc = []
        for j in av:
            g, k, r = features(R, st, j, K)
            sc.append((wts[0] * (-g) + wts[1] * k + wts[2] * idle[j] + wts[3] * r * 10
                       + wts[4] * (1.0 if g > W else 0.0), j))
        sc.sort(reverse=True); ch = [j for _, j in sc[:K]]
        for j in range(len(R)): idle[j] = 0 if j in ch else idle[j] + (st[j][0] == 0)
        st = stepK(R, st, ch)
    return False, maxstep, 0   # 未在上限内进入环：不算证书

def policy_ga(R, K, budget, rng, pop=20):
    n = 5
    def ev(x):
        ok, t, per = run_policy(R, K, x)
        return (2 if ok else 0) * 10**6 + t, per
    X = rng.normal(0, 1, (pop, n)); F = [ev(x) for x in X]; used = pop
    while used < budget:
        b = max(range(pop), key=lambda i: F[i][0])
        if F[b][0] >= 2 * 10**6: return True, used, F[b][1]
        kids = []
        for _ in range(pop):
            a, c = [max(rng.integers(pop, size=2), key=lambda i: F[i][0]) for _ in range(2)]
            ch = np.where(rng.random(n) < 0.5, X[a], X[c]) + rng.normal(0, 0.3, n) * (rng.random(n) < 0.4)
            kids.append(ch)
        FK = [ev(x) for x in kids]; used += pop
        allp = list(zip(list(X) + kids, F + FK)); allp.sort(key=lambda z: -z[1][0])
        X = [z[0] for z in allp[:pop]]; F = [z[1] for z in allp[:pop]]
    return False, used, 0
