import os, sys; _h=os.path.dirname(os.path.abspath(__file__)); sys.path[:0]=[os.path.join(_h,'..','启发式验证'), os.path.join(_h,'..')]
"""D4 修正：先算单台机器人的最优长期间隔 ā*（允许提前维修时 = 单机状态图上的最小平均环，Karp 算法），
下界 ⌈K·ā*⌉；再从 min(强制下界, 提前下界) 开始精确搜索。"""
import math, json, sys
from d4 import params, options, feasible

def single_graph(P, early):
    full = tuple(m[0] for m in P["maint"])
    S = {full: []}; fr = [full]
    while fr:
        s = fr.pop()
        for gap, nl in options(P, s, early):
            S[s].append((gap, nl))
            if nl not in S: S[nl] = []; fr.append(nl)
    return S

def min_mean_cycle(S):
    nodes = list(S); idx = {v: i for i, v in enumerate(nodes)}; n = len(nodes)
    INF = float("inf"); D = [[INF] * n for _ in range(n + 1)]
    for v in range(n): D[0][v] = 0
    for k in range(1, n + 1):
        for u in nodes:
            du = D[k - 1][idx[u]]
            if du == INF: continue
            for w, v in S[u]:
                if du + w < D[k][idx[v]]: D[k][idx[v]] = du + w
    best = INF
    for v in range(n):
        if D[n][v] == INF: continue
        best = min(best, max((D[n][v] - D[k][v]) / (n - k) for k in range(n) if D[k][v] < INF))
    return best

def nstar(P, K, early, lo, span=6):
    for N in range(lo, lo + span):
        r = feasible(P, N, K, early, lim=600_000)
        if r is None: return ("?", N)
        if r: return N
    return ("超出范围", lo + span)

out = []
for name in ["L1-720", "L2-720", "L3-720"]:
    for r in [0, 1, 2]:
        P = params(name, r)
        a_f = min_mean_cycle(single_graph(P, False)); a_e = min_mean_cycle(single_graph(P, True))
        for K in [1, 2]:
            lb_f = math.ceil(K * a_f - 1e-9); lb_e = math.ceil(K * a_e - 1e-9)
            nf = nstar(P, K, False, lb_f); ne = nstar(P, K, True, lb_e)
            line = dict(inst=name, type=r, K=K, 强制_平均间隔=round(a_f, 3), 提前_平均间隔=round(a_e, 3),
                        强制_下界=lb_f, 提前_下界=lb_e, 强制_最少=nf, 提前_最少=ne)
            out.append(line); print(json.dumps(line, ensure_ascii=False), flush=True)
json.dump(out, open("d4b_results.json", "w"), ensure_ascii=False, indent=1)
