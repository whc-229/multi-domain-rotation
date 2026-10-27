"""全新启动的无穷时域上界：时间展开 [0, T+P+gmax)，并强制出发时刻 t∈[T, T+gmax) 的每条弧流量等于 t+P 处的流量。
这样的解可按周期 P 无限延拓（证明见笔记），所以给出 N_fresh(K) 的严格上界；与有限时域下界相等时即为精确值。"""
import numpy as np, math, json, sys
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix

def fresh_tail(g, K, T, P, reset=False, time_limit=120):
    L = len(g); gmax = max(g + [1]); H = T + P + gmax; arcs = []
    for t in range(H):
        for f in range(L):
            arcs.append(("w", t, f, t + g[f], (f + 1) % L))
            arcs.append(("z", t, f, t + 1, f))
            if reset and f != L - 1: arcs.append(("r", t, f, t + g[L - 1], 0))
    pos = {(k, t, f): j for j, (k, t, f, ht, hf) in enumerate(arcs)}
    n = len(arcs) + 1; node = lambda t, f: t * L + f
    A = lil_matrix((H * L, n))
    for j, (k, t, f, ht, hf) in enumerate(arcs):
        A[node(t, f), j] -= 1
        if ht < H: A[node(ht, hf), j] += 1
    A[node(0, 0), n - 1] = 1
    C = lil_matrix((H, n))
    for j, (k, t, f, ht, hf) in enumerate(arcs):
        if k in ("w", "r"): C[t, j] = 1
    rows = []
    for (k, t, f), j in pos.items():
        if T <= t < T + gmax and (k, t + P, f) in pos: rows.append((j, pos[(k, t + P, f)]))
    E = lil_matrix((len(rows), n))
    for r, (j1, j2) in enumerate(rows): E[r, j1] = 1; E[r, j2] = -1
    c = np.zeros(n); c[-1] = 1
    res = milp(c=c, constraints=[LinearConstraint(A.tocsr(), 0, 0), LinearConstraint(C.tocsr(), K, K), LinearConstraint(E.tocsr(), 0, 0)],
               integrality=np.ones(n), bounds=Bounds(0, np.inf), options={"time_limit": time_limit})
    if res.x is None: return None
    return int(round(res.x[-1]))

if __name__ == "__main__":
    reset = len(sys.argv) > 1 and sys.argv[1] == "reset"
    Kmax = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    seqs = [[3,10],[2,9],[4,4,5],[3,5,2,6],[2,2,5],[3,3,15],[4,6,5,4],[2,6,3],[4,4,5,4,9],[3,6,3,6,9]]
    for g in seqs:
        S = sum(g); L = len(g); res = []
        for K in range(1, Kmax + 1):
            q = L // math.gcd(K, L); best = None
            for T in [S, 2 * S]:
                for P in sorted(set([q * math.ceil(S / q), q * math.ceil(max(g) / q), 2 * q * math.ceil(S / q)])):
                    v = fresh_tail(g, K, T, P, reset)
                    if v is not None and (best is None or v < best[0]): best = (v, T, P)
            res.append(best)
        print(json.dumps(dict(g=g, reset=reset, upper=[b[0] if b else None for b in res], TP=[b[1:] if b else None for b in res])), flush=True)
