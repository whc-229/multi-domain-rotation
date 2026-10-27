import numpy as np, sys
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from fresh import fresh_min
import fresh as F
def fresh_lp(g, K, H):
    # 复用 fresh_min 的构造，但求 LP 松弛
    import scipy.optimize as so
    orig = so.milp
    L = len(g); arcs = []
    for t in range(H):
        for f in range(L):
            arcs.append(("w", t, f, t + g[f], (f + 1) % L)); arcs.append(("z", t, f, t + 1, f))
    from scipy.sparse import lil_matrix
    n = len(arcs) + 1; node = {(t, f): t * L + f for t in range(H) for f in range(L)}
    A = lil_matrix((H * L, n))
    for j, (k, t, f, ht, hf) in enumerate(arcs):
        A[node[(t, f)], j] -= 1
        if ht < H: A[node[(ht, hf)], j] += 1
    A[node[(0, 0)], n - 1] = 1
    C = lil_matrix((H, n))
    for j, (k, t, f, ht, hf) in enumerate(arcs):
        if k == "w": C[t, j] = 1
    c = np.zeros(n); c[-1] = 1
    res = milp(c=c, constraints=[LinearConstraint(A.tocsr(), 0, 0), LinearConstraint(C.tocsr(), K, K)],
               integrality=np.zeros(n), bounds=Bounds(0, np.inf))
    return res.fun
for g in [[4,4,5],[3,10],[2,2,5],[3,5,2,6],[4,4,5,4,9]]:
    S = sum(g); a = S / len(g)
    print(g, "ā=", round(a,3), "全新启动 LP(K=1) 随 H:", [round(fresh_lp(g, 1, H), 3) for H in (S, 2*S, 4*S, 8*S)])
