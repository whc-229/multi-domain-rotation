"""假设 H（均匀拉伸）：K=1 时，稳态最少台数达到 ⌈ā⌉ 当且仅当存在拉伸 g' ≥ g（总拉伸 = nL − S），
使 A(g') 的 n 个平移恰好铺满 Z_{nL}。用图谱结果检验；一般 K 时检验 M = N·L/K 为整数的情形。"""
import json, math, itertools, sys
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

def tiles(gp, M, K, N):
    A = [0]
    for x in gp[:-1]: A.append(A[-1] + x)
    if len(set(a % M for a in A)) < len(A): return False
    Mx = np.zeros((M, M))
    for s in range(M):
        for a in A: Mx[s, (s - a) % M] += 1
    res = milp(c=np.ones(M), constraints=[LinearConstraint(Mx, K, K), LinearConstraint(np.ones((1, M)), N, N)],
               integrality=np.ones(M), bounds=Bounds(0, np.inf), options={"time_limit": 20})
    return res.x is not None

def uniform_exists(g, K):
    L, S = len(g), sum(g); N = math.ceil(K * S / L - 1e-9)
    if (N * L) % K: return None          # 均匀模式不可能（周期不是整数）
    M = N * L // K; B = M - S
    for cut in itertools.combinations_with_replacement(range(L), B):
        e = [0] * L
        for c in cut: e[c] += 1
        if tiles([gi + ei for gi, ei in zip(g, e)], M, K, N): return True
    return False

if __name__ == "__main__":
    data = json.load(open(sys.argv[1]))
    agree = dis = na = 0; examples = []
    for r in data:
        g = r["g"]
        for row in r["rows"]:
            K = row["K"]; gap = row["N_wait"] - row["LB"] if row["N_wait"] is not None else None
            if gap is None: continue
            u = uniform_exists(g, K)
            if u is None: na += 1; continue
            if (u and gap == 0) or ((not u) and gap >= 1): agree += 1
            else: dis += 1; examples.append((g, K, gap, u))
    print("一致", agree, "不一致", dis, "不适用(M 非整数)", na)
    print("不一致例子", examples[:30])
