"""全新启动流体下界的精确对偶证书。
原问题（K=1 归一化，时域 H）：min N，s.t. 流守恒（(0,0) 注入 N），每格工作弧流量和 = 1，x>=0。
对偶：max Σ_t y_t  s.t. π_00 <= 1；对每条弧 u->v（v 越界则 π_v=0）：-π_u + π_v + [工作弧] y_t <= 0。
任一可行对偶解给出 N_fresh(K) >= K·Σy（对所有 K、无穷时域都成立，因为截断是松弛）。
这里用 HiGHS 解对偶，把解有理化，再用 Fraction 精确核对可行性与目标值。"""
import numpy as np, json, sys
from fractions import Fraction as F
from scipy.optimize import linprog
from scipy.sparse import lil_matrix

def arcs_of(g, H, reset=False):
    L = len(g); arcs = []
    for t in range(H):
        for f in range(L):
            arcs.append((t, f, t + g[f], (f + 1) % L, True))
            arcs.append((t, f, t + 1, f, False))
            if reset and f != L - 1: arcs.append((t, f, t + g[L - 1], 0, True))
    return arcs

def dual_cert(g, H, reset=False, maxden=None):
    L = len(g); arcs = arcs_of(g, H, reset)
    nn = H * L; nv = nn + H          # π 变量在前，y 在后
    idx = lambda t, f: t * L + f
    A = lil_matrix((len(arcs) + 1, nv)); b = np.zeros(len(arcs) + 1)
    for j, (t, f, ht, hf, work) in enumerate(arcs):
        A[j, idx(t, f)] = -1
        if ht < H: A[j, idx(ht, hf)] = 1
        if work: A[j, nn + t] = 1
    A[len(arcs), idx(0, 0)] = 1; b[-1] = 1
    c = np.zeros(nv); c[nn:] = -1
    res = linprog(c, A_ub=A.tocsr(), b_ub=b, bounds=[(None, None)] * nv, method="highs")
    z = res.x; val = -res.fun
    for md in ([maxden] if maxden else [12, 60, 420, 2520, 27720, 10**6]):
        q = [F(v).limit_denominator(md) for v in z]
        ok = all(-q[idx(t, f)] + (q[idx(ht, hf)] if ht < H else 0) + (q[nn + t] if work else 0) <= 0 for (t, f, ht, hf, work) in arcs) and q[idx(0, 0)] <= 1
        if ok:
            return sum(q[nn:]), val, md, q
    return None, val, None, None

if __name__ == "__main__":
    seqs = [[3,10],[2,9],[4,4,5],[3,5,2,6],[2,2,5],[3,3,15],[4,4,5,4,9],[3,6,3,6,9],[4,6,5,4],[2,6,3],[5,5,5,21],[1,4],[2,4],[1,1,5],[2,2,7]]
    out = {}
    for g in seqs:
        S = sum(g)
        for reset in [False, True]:
            bound, val, md, q = dual_cert(g, 2 * S, reset)
            abar = F(S, len(g))
            out[f"{g}{'R' if reset else ''}"] = dict(bound=str(bound), lp=val, abar=str(abar))
            print(g, "提前大修" if reset else "仅等待", "ā=", abar, "精确对偶下界=", bound, "(LP", round(val, 6), ")", "超出ā:" , (bound - abar) if bound is not None else None, flush=True)
    json.dump(out, open("cert.json", "w"), ensure_ascii=False, indent=1)
