import numpy as np, sys, json
from fractions import Fraction
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
def fresh_lp(g, K, H, reset=False):
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
    res = milp(c=c, constraints=[LinearConstraint(A.tocsr(), 0, 0), LinearConstraint(C.tocsr(), K, K)],
               integrality=np.zeros(n), bounds=Bounds(0, np.inf))
    return res.fun
if __name__ == "__main__":
    out = []
    for a in range(1, 7):
        for G in range(a, 21):
            S = a + G; v = fresh_lp([a, G], 1, 6 * S)
            fr = Fraction(v).limit_denominator(200)
            out.append((a, G, round(v, 5), str(fr), round(v - S / 2, 4)))
    for o in out: print(o)
    json.dump(out, open("fresh_lp_L2.json", "w"))
