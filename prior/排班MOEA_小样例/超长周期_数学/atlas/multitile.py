"""零余量：N_w(K) = Kā ⇔ A 在 Z_S 上 K 重平移铺满。用小 ILP 直接判定 K 重铺满（B_t 为非负整数，Σ_a B_{x-a} = K）。
枚举 L=3,4,5,6 的间隔序列（循环同构去重），统计 T_A = {K ≤ L : K 重可铺} 的形状。"""
import itertools, math, json, sys
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from cm_check import cm_tiles_ZS
def ktile(A, S, K):
    M = np.zeros((S, S))
    for x in range(S):
        for a in A: M[x, (x - a) % S] += 1
    r = milp(c=np.zeros(S), constraints=[LinearConstraint(M, K, K)], integrality=np.ones(S), bounds=Bounds(0, np.inf), options={"time_limit": 20})
    return r.status == 0
def canon(g):
    L = len(g); return min(tuple(g[i:] + g[:i]) for i in range(L))
res = {}; shapes = {}
for L, Smax in [(3, 24), (4, 24), (5, 20), (6, 24)]:
    seen = set()
    for S in range(L, Smax + 1, L):
        for comp in itertools.product(range(1, S), repeat=L - 1):
            last = S - sum(comp)
            if last < 1: continue
            g = list(comp) + [last]; c = canon(g)
            if c in seen: continue
            seen.add(c)
            A = [sum(g[:i]) for i in range(L)]
            T = tuple(K for K in range(1, L + 1) if ktile(A, S, K))
            cm, SA = cm_tiles_ZS(A, S)
            key = (L, T); shapes[key] = shapes.get(key, 0) + 1
            if (1 in T) != cm: print("CM 不一致", g, T, SA)
            # 形状是否为 m·N
            m = T[0] if T else None
            if m and T != tuple(range(m, L + 1, m)): print("非 mN 形状", g, T)
    print(L, {k[1]: v for k, v in shapes.items() if k[0] == L}, flush=True)
