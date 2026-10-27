"""多型号 × 多地长时域派工仿真（只算机器人）。与 mtms.Inst 同一网格模型：
机器人状态 (机型 k, 相位 i, 虚拟就绪时刻 v)；站点 s 在 t 时刻可接收满足 v + lead_s ≤ t 的允许机型机器人；
值岗首尾相接：站点当前值岗结束时刻 free[s] 到达时必须立刻有人接班，否则该格缺岗（下一格再试）。
派工策略决定"在可用机器人中选谁"；策略本身与时域长度无关。"""
import math, heapq
import numpy as np
from mtms import fluid_lb

def simulate(I, Nk, policy, H, warm, fresh=True, params=None, return_trace=False, cycle=False):
    nk = len(I.T); ns = len(I.sites)
    v0 = [math.ceil(t["off_h"] / 2.0 - 1e-9) for t in I.T]
    robots = []   # [k, i, v]
    for k in range(nk):
        for _ in range(int(Nk[k])): robots.append([k, 0, v0[k] if fresh else 0])
    c0 = max(v0) if fresh else 0
    free = [c0 + I.lead(s) for s in range(ns)]
    lead = [I.lead(s) for s in range(ns)]
    # 预计算 step
    st = {(k, i, s): I.step(k, i, s) for k in range(nk) for i in range(I.L[k]) for s in I.elig[k]}
    D = [[0] * I.L[k] for k in range(nk)]    # 各型各相位累计派工数（相位平衡策略用）
    tot = [0] * nk
    gaps = 0; gap_sites = [0] * ns
    nelig_site = [sum(1 for k in range(nk) if s in I.elig[k]) for s in range(ns)]
    flex = [len(I.elig[k]) for k in range(nk)]
    t = min(free)
    seen = {}; cyc = None; gaps_at = {}
    while t < H:
        if cycle:
            key = (tuple((r[0], r[1], r[2] - t) for r in robots), tuple(f - t for f in free),
                   tuple(tuple(tot[k] - I.L[k] * D[k][i] for i in range(I.L[k])) for k in range(nk)))
            if key in seen:
                cyc = (seen[key], t); break
            seen[key] = t; gaps_at[t] = gaps
        due = [s for s in range(ns) if free[s] <= t]
        # 限制多的站点优先挑
        due.sort(key=lambda s: nelig_site[s])
        for s in due:
            cand = [r for r in robots if s in I.elig[r[0]] and r[2] + lead[s] <= t]
            if not cand:
                if t >= warm: gaps += 1; gap_sites[s] += 1
                free[s] = t + 1; continue
            if policy == "fifo":
                r = min(cand, key=lambda r: (r[2], r[0]))
            elif policy == "short":     # 下一次占用最短者优先（吞吐量贪心）
                r = min(cand, key=lambda r: (st[(r[0], r[1], s)][1] / I.d[r[0]], r[2]))
            elif policy == "pb":        # 相位平衡：同型内补派工缺额最大的相位；型间选灵活度低、等待久的
                def key(r):
                    k, i = r[0], r[1]
                    deficit = tot[k] / I.L[k] - D[k][i]
                    return (flex[k], -deficit, r[2])
                r = min(cand, key=key)
            elif policy == "score":     # 参数化打分（供进化）：w·特征，越小越优
                w = params
                def key(r):
                    k, i = r[0], r[1]
                    nxt = st[(k, i, s)][1]
                    deficit = tot[k] / I.L[k] - D[k][i]
                    return (w[0] * (r[2] - t) + w[1] * nxt / I.d[k] + w[2] * flex[k] - w[3] * deficit
                            + w[4] * I.d[k] + w[5] * i / I.L[k])
                r = min(cand, key=key)
            else:
                raise ValueError(policy)
            k, i = r[0], r[1]
            ld, nx = st[(k, i, s)]
            start = t
            r[2] = start - ld + nx    # 新的虚拟就绪时刻：值岗开始 D = v' + ld → 下一节点 = D − ld + nx
            r[1] = (i + 1) % I.L[k]
            D[k][i] += 1; tot[k] += 1
            free[s] = start + I.d[k]
        t = min(free)
    if cycle:
        if cyc is None: return None, None, None            # 时域内未出现状态重复
        t1, t2 = cyc
        return (gaps_at[t1] == gaps), (t1, t2 - t1), gaps  # 环内无新增缺岗 ⇒ 无穷时域可行（确定性系统按周期 t2−t1 重复）
    return gaps, gap_sites

def size_fleet(I, policy, H, warm, fresh=True, params=None, start=None, max_add=40):
    """从整数下界的机型配比出发，贪心逐台加机器人直到无缺岗。返回 (总台数, 各型台数)。"""
    if start is None:
        _, x = fluid_lb(I, integer=True); start = [int(round(v)) for v in x]
    Nk = list(start)
    g, _ = simulate(I, Nk, policy, H, warm, fresh, params)
    added = 0
    while g > 0 and added < max_add:
        best = None
        for k in range(len(I.T)):
            Nk[k] += 1; gg, _ = simulate(I, Nk, policy, H, warm, fresh, params); Nk[k] -= 1
            if best is None or gg < best[0]: best = (gg, k)
        Nk[best[1]] += 1; g = best[0]; added += 1
    # 回删：尝试逐台去掉仍无缺岗
    improved = True
    while improved and g == 0:
        improved = False
        for k in range(len(I.T)):
            if Nk[k] == 0: continue
            Nk[k] -= 1; gg, _ = simulate(I, Nk, policy, H, warm, fresh, params)
            if gg == 0: improved = True; break
            Nk[k] += 1
    return (sum(Nk) if g == 0 else None), Nk
