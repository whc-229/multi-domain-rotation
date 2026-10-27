"""相位提升的周期时间展开网络（稳态，允许等待）。
节点 (t, φ)：第 t 格（mod P）开始时、处于相位 φ 且已就绪的机器人。
弧：工作 x[t,φ]：(t,φ) → (t+g_φ, φ+1)，长度 g_φ；等待 z[t,φ]：(t,φ) → (t+1, φ)，长度 1。
可选：提前大修 r[t,φ]（工作后直接进大修回到相位 0）：(t,φ) → (t+g_{L-1}, 0)，长度 g_{L-1}，同样计入覆盖。
覆盖：Σ_φ (x+r)[t,φ] = K（exact=True）或 ≥ K（允许空跑）。
机器人数 N = Σ 流量·长度 / P（循环流的每个圈长度都是 P 的倍数）。"""
import numpy as np, math
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix

def build(g, K, P, exact=True, reset=False):
    L = len(g)
    arcs = []   # (kind, t, φ, head_t, head_φ, length)
    for t in range(P):
        for f in range(L):
            arcs.append(("w", t, f, (t + g[f]) % P, (f + 1) % L, g[f]))
            arcs.append(("z", t, f, (t + 1) % P, f, 1))
            if reset and f != L - 1:
                arcs.append(("r", t, f, (t + g[L - 1]) % P, 0, g[L - 1]))
    n = len(arcs); idx = {(t, f): t * L + f for t in range(P) for f in range(L)}
    A_eq = lil_matrix((P * L, n)); 
    for j, (k, t, f, ht, hf, ln) in enumerate(arcs):
        A_eq[idx[(t, f)], j] -= 1; A_eq[idx[(ht, hf)], j] += 1
    C = lil_matrix((P, n))
    for j, (k, t, f, ht, hf, ln) in enumerate(arcs):
        if k in ("w", "r"): C[t, j] = 1
    cost = np.array([a[5] for a in arcs], float) / P
    cons = [LinearConstraint(A_eq.tocsr(), 0, 0),
            LinearConstraint(C.tocsr(), K, K if exact else np.inf)]
    return arcs, cost, cons

def solve(g, K, P, exact=True, reset=False, integer=True, time_limit=60):
    arcs, cost, cons = build(g, K, P, exact, reset)
    res = milp(c=cost, constraints=cons, integrality=np.ones(len(cost)) if integer else np.zeros(len(cost)),
               bounds=Bounds(0, np.inf), options={"time_limit": time_limit})
    if res.status not in (0,) or res.x is None: return None, None
    return res.fun, res.x

def decompose(g, P, x_arcs, arcs):
    """把整数循环流拆成机器人轨迹并逐格核验：返回 (机器人数, 是否通过)"""
    L = len(g)
    flow = {}
    for v, a in zip(np.round(x_arcs).astype(int), arcs):
        if v: flow[a] = flow.get(a, 0) + v
    out = {}
    for a, v in flow.items(): out.setdefault((a[1], a[2]), []).extend([a] * v)
    total_len = 0; robots = 0
    while any(out.values()):
        start = next(k for k, v in out.items() if v); node = start; length = 0
        while True:
            a = out[node].pop(); length += a[5]; node = (a[3], a[4])
            if node == start: break
        if length % P: return None, False
        robots += length // P; total_len += length
    return robots, True
