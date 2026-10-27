"""全新启动、有限时域 H 的时间展开 ILP：第 0 格 N 台机器人全部在相位 0 就绪；每格恰好 K 台工作；
工作弧 (t,φ)→(t+g_φ, φ+1)，等待弧 (t,φ)→(t+1,φ)，可选提前大修弧 (t,φ)→(t+g_{L-1}, 0)；越过 H 的弧流出到汇点。
min N。N_H 随 H 单调不减；H 足够大时的值记为 N_fresh。"""
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix

def fresh_min(g, K, H, reset=False, time_limit=120):
    L = len(g); arcs = []
    for t in range(H):
        for f in range(L):
            arcs.append(("w", t, f, t + g[f], (f + 1) % L))
            arcs.append(("z", t, f, t + 1, f))
            if reset and f != L - 1: arcs.append(("r", t, f, t + g[L - 1], 0))
    n = len(arcs) + 1   # 最后一个变量是 N
    node = {(t, f): t * L + f for t in range(H) for f in range(L)}
    A = lil_matrix((H * L, n))
    for j, (k, t, f, ht, hf) in enumerate(arcs):
        A[node[(t, f)], j] -= 1
        if ht < H: A[node[(ht, hf)], j] += 1
    A[node[(0, 0)], n - 1] = 1            # 初始供给 N 台（相位 0）
    C = lil_matrix((H, n))
    for j, (k, t, f, ht, hf) in enumerate(arcs):
        if k in ("w", "r"): C[t, j] = 1
    c = np.zeros(n); c[-1] = 1
    res = milp(c=c, constraints=[LinearConstraint(A.tocsr(), 0, 0), LinearConstraint(C.tocsr(), K, K)],
               integrality=np.ones(n), bounds=Bounds(0, np.inf), options={"time_limit": time_limit})
    if res.x is None: return None
    return int(round(res.x[-1]))
