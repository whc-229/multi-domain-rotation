"""启动代价试点：全新机队同时上线（全部相位 0）与稳态（相位已错开）的最少台数差，以及"提前大修"能否消掉这个差。
模型同 mtms.Inst（只算机器人）。三类弧：
  z  原地等 1 格；
  w  在站点 s 值岗，按本相位的占用 occ[i] 回来，相位 +1；
  r  （reset=True 时）在相位 i ≠ L−1 值岗后直接做最高级维修，占用 occ[L−1]，回到相位 0。
稳态：周期网络（周期 P）。全新启动：
  fresh_lb(H)        只要求 [c0, H) 覆盖 —— 无限期全新启动最少台数的下界；
  fresh_ub(T, P)     时间展开到 T+P，强制 [T−G, T) 与 [T+P−G, T+P) 的弧流量相同，于是 T 之后按周期 P 永远重复 —— 可无限延拓，是上界。"""
import math, sys, json, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
from mtms import Inst, U

def arc_list(I, k, i, reset):
    """(kind, s, dv, cov_start, cov_len, next_phase)"""
    L = I.L[k]; out = [("z", -1, 1, 0, 0, i)]
    for s in I.elig[k]:
        ld, nx = I.step(k, i, s)
        out.append(("w", s, nx, ld, I.d[k], (i + 1) % L))
        if reset and i != L - 1:
            nxr = ld + math.ceil((I.T[k]["occ_h"][L - 1] + I.sites[s]["delta"]) / U - 1e-9)
            out.append(("r", s, nxr, ld, I.d[k], 0))
    return out

def steady(I, P, reset=False, time_limit=120):
    nk = len(I.T); arcs = []
    for k in range(nk):
        for i in range(I.L[k]):
            for a in arc_list(I, k, i, reset):
                for v in range(P): arcs.append((k, i, v) + a)
    node = {}
    for k in range(nk):
        for i in range(I.L[k]):
            for v in range(P): node[(k, i, v)] = len(node)
    n = len(arcs); A = lil_matrix((len(node), n)); C = lil_matrix((len(I.sites) * P, n))
    for j, (k, i, v, kind, s, dv, cs, cl, ni) in enumerate(arcs):
        A[node[(k, i, v)], j] -= 1; A[node[(k, ni, (v + dv) % P)], j] += 1
        for q in range(cl): C[s * P + (v + cs + q) % P, j] += 1
    Kv = np.repeat([s["K"] for s in I.sites], P).astype(float)
    cost = np.array([a[5] / P for a in arcs])
    r = milp(cost, constraints=[LinearConstraint(A.tocsr(), 0, 0), LinearConstraint(C.tocsr(), Kv, Kv)],
             integrality=np.ones(n), bounds=Bounds(0, np.inf), options={"time_limit": time_limit, "mip_rel_gap": 0})
    return None if r.x is None else round(r.fun, 6)

def fresh(I, T, P=0, reset=False, time_limit=300, tail=True, return_sol=False):
    """tail=False：下界模型，只要求 [c0, T) 覆盖。tail=True：上界模型，覆盖到 T+P 并强制末段周期 P。"""
    nk = len(I.T); H = T + P if tail else T
    v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]; c0 = max(v0)
    G = max(a[2] for k in range(nk) for i in range(I.L[k]) for a in arc_list(I, k, i, reset))
    if tail: assert T - G >= max(v0), (T, G)
    arcs = []; idx = {}
    for k in range(nk):
        for i in range(I.L[k]):
            for ai, a in enumerate(arc_list(I, k, i, reset)):
                for v in range(v0[k], H):
                    idx[(k, i, ai, v)] = len(arcs); arcs.append((k, i, v) + a)
    node = {}
    for k in range(nk):
        for i in range(I.L[k]):
            for v in range(v0[k], H): node[(k, i, v)] = len(node)
    na = len(arcs); n = na + nk
    A = lil_matrix((len(node), n))
    rows = {}
    for s in range(len(I.sites)):
        for t in range(c0 + I.lead(s), H): rows[(s, t)] = len(rows)
    C = lil_matrix((len(rows), n))
    for j, (k, i, v, kind, s, dv, cs, cl, ni) in enumerate(arcs):
        A[node[(k, i, v)], j] -= 1
        if v + dv < H: A[node[(k, ni, v + dv)], j] += 1
        for q in range(cl):
            if (s, v + cs + q) in rows: C[rows[(s, v + cs + q)], j] += 1
    for k in range(nk): A[node[(k, 0, v0[k])], na + k] = 1
    Kv = np.array([I.sites[s]["K"] for (s, t) in rows], float)
    cons = [LinearConstraint(A.tocsr(), 0, 0), LinearConstraint(C.tocsr(), Kv, Kv)]
    if tail:
        R = lil_matrix((len(idx), n)); nr = 0
        for (k, i, ai, v), j in idx.items():
            if T - G <= v < T:
                R[nr, j] = 1; R[nr, idx[(k, i, ai, v + P)]] = -1; nr += 1
        cons.append(LinearConstraint(R[:nr].tocsr(), 0, 0))
    cost = np.zeros(n); cost[na:] = 1
    r = milp(cost, constraints=cons, integrality=np.ones(n), bounds=Bounds(0, np.inf),
             options={"time_limit": time_limit, "mip_rel_gap": 0})
    bound = getattr(r, "mip_dual_bound", None)
    if r.x is None: return dict(obj=None, bound=bound, status=r.status)
    out = dict(obj=round(r.fun), bound=None if bound is None else round(bound, 3), status=r.status,
               N=[int(round(x)) for x in r.x[na:]])
    if return_sol:
        out["arcs"] = [(a, int(round(x))) for a, x in zip(arcs, r.x[:na]) if x > 0.5 and a[3] != "z"]
    return out

def G_of(I, reset):
    return max(a[2] for k in range(len(I.T)) for i in range(I.L[k]) for a in arc_list(I, k, i, reset))
