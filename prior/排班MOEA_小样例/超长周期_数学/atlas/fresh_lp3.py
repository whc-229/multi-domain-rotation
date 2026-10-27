"""全新启动的流体值 α（时间展开 LP，K=1，H=6S）在四种灵活性下：仅等待 / 提前大修 / 空跑（覆盖≥K）/ 两者。"""
import numpy as np, json
from fractions import Fraction
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
def fresh_lp(g, H, reset=False, dummy=False):
    L = len(g); arcs = []
    for t in range(H):
        for f in range(L):
            arcs.append(("w", t, f, t + g[f], (f + 1) % L)); arcs.append(("z", t, f, t + 1, f))
            if reset and f != L - 1: arcs.append(("r", t, f, t + g[L - 1], 0))
    n = len(arcs) + 1; node = {(t, f): t * L + f for t in range(H) for f in range(L)}
    A = lil_matrix((H * L, n))
    for j, (k, t, f, ht, hf) in enumerate(arcs):
        A[node[(t, f)], j] -= 1
        if ht < H: A[node[(ht, hf)], j] += 1
    A[node[(0, 0)], n - 1] = 1
    C = lil_matrix((H, n))
    for j, (k, t, f, ht, hf) in enumerate(arcs):
        if k in ("w", "r"): C[t, j] = 1
    c = np.zeros(n); c[-1] = 1
    res = milp(c=c, constraints=[LinearConstraint(A.tocsr(), 0, 0), LinearConstraint(C.tocsr(), 1, np.inf if dummy else 1)],
               integrality=np.zeros(n), bounds=Bounds(0, np.inf))
    return res.fun
if __name__ == "__main__":
    out = {}
    for g in [[3,10],[2,9],[4,4,5],[3,5,2,6],[2,2,5],[3,3,15],[4,4,5,4,9],[3,6,3,6,9],[4,6,5,4],[2,6,3],[5,5,5,21]]:
        S = sum(g); a = S / len(g); H = 6 * S
        vals = {k: fresh_lp(g, H, r, d) for k, (r, d) in dict(仅等待=(False, False), 提前大修=(True, False), 空跑=(False, True), 两者=(True, True)).items()}
        out[str(g)] = {k: round(v, 4) for k, v in vals.items()}
        print(g, "ā=", round(a, 4), {k: f"{v:.4f}({Fraction(v).limit_denominator(300)})" for k, v in vals.items()}, flush=True)
    json.dump(out, open("fresh_lp_variants.json", "w"), ensure_ascii=False, indent=1)
    