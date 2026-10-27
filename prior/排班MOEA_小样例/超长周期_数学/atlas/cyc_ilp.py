"""固定模式稳态：循环覆盖整数规划。x_t = 在偏移 t 处摆放的机器人数（按间隔序列 g 严格循环，不等待）。
min Σx_t  s.t. 每格 s：Σ_{a∈A} x_{(s−a) mod S} ≥ K。用 scipy.milp（HiGHS）。"""
import numpy as np, math, json, sys
from scipy.optimize import milp, LinearConstraint, Bounds

def n_fixed(g, K):
    S = sum(g); A = [0]
    for x in g[:-1]: A.append(A[-1] + x)
    M = np.zeros((S, S))
    for s in range(S):
        for a in A: M[s, (s - a) % S] += 1
    res = milp(c=np.ones(S), constraints=LinearConstraint(M, lb=K * np.ones(S), ub=np.inf),
               integrality=np.ones(S), bounds=Bounds(0, np.inf))
    return int(round(res.fun))

if __name__ == "__main__":
    seqs = [[3,10],[2,9],[4,4,5],[3,5,2,6],[2,2,5],[3,3,15],[4,4,5,4,9],[3,6,3,6,9],[5,5,5,21],[3,3,15,4]]
    out = {}
    for g in seqs:
        S, L = sum(g), len(g); a = S / L
        Ns = [n_fixed(g, K) for K in range(1, 61)]
        gaps = [N - math.ceil(K * a - 1e-9) for K, N in zip(range(1, 61), Ns)]
        # 最终周期性：N(K+L) = N(K)+S 从哪个 K0 起对所有 K ≤ 60-L 成立
        K0 = None
        for k0 in range(1, 61 - L):
            if all(Ns[k + L - 1] == Ns[k - 1] + S for k in range(k0, 61 - L)): K0 = k0; break
        out[str(g)] = dict(S=S, L=L, abar=round(a, 3), N_1_12=Ns[:12], gap_1_24=gaps[:24],
                           limsup_gap=max(gaps[30:]), liminf_gap=min(gaps[30:]), 周期性起点K0=K0)
        print(json.dumps({"g": g, **out[str(g)]}, ensure_ascii=False), flush=True)
    json.dump(out, open("atlas_fixed.json", "w"), ensure_ascii=False, indent=1)
