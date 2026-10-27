"""多型号 × 多地 × 站点限制机型（MT-MS）排班模型，只算机器人（工具以后再加）。

时间网格 u = 2h（值岗 10h = 5 格，6h = 3 格）。时间变量取"虚拟时刻" v：机器人可以在一个零额外路程站点开始值岗的最早时刻。
- 机型 τ：值岗格数 d_τ，占用序列 occ_τ[i]（小时，含准备、装工具、去程、值岗、收尾、返程、充电/维修，按维修超周期循环），相位 i。
- 站点 s：允许机型集合 E_s，额外单程路程 δ_s（小时，去回各一次），每格需要 K_s 台在岗（恰好覆盖，值岗首尾相接）。
- 机器人在节点 (v, τ, i) 选站点 s 值岗：值岗开始 D = v + ⌈δ_s/u⌉，覆盖 [D, D + d_τ)，下一节点 v′ = D + ⌈(occ_τ[i] + δ_s)/u⌉，相位 i+1。
  所有取整都向上，所以网格模型是连续时间模型的保守近似（每班最多损失 1 格）。
- 也可以原地等待：v → v+1。
"""
import math, json
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix
from roster_sim import INSTANCES
from lb import robot_occ_seq

U = 2.0

def real_types():
    """从 model-crow 的三个实例取出 9 种机型（同一实例内 3 种）。"""
    out = {}
    for nm in ["L1-720", "L2-720", "L3-720"]:
        inst = INSTANCES[nm]()
        for r, p in enumerate(inst.robots):
            seq, mu = robot_occ_seq(inst, r)
            out[f"{nm[:2]}-{r}"] = dict(duty_h=p["duty"], off_h=p["prep"] + p["max_install"] + p["go"], occ_h=list(seq[mu:]),
                                        max_num=p["max_num"])
    return out

class Inst:
    def __init__(self, types, sites, cost=None, name=""):
        """types: dict 名 -> 机型参数；sites: list of dict(E=[机型名], delta=小时, K=台数)"""
        self.tn = list(types); self.T = [types[t] for t in self.tn]; self.sites = sites; self.name = name
        self.cost = cost or {t: 1.0 for t in self.tn}
        self.d = [int(round(t["duty_h"] / U)) for t in self.T]
        assert all(abs(t["duty_h"] / U - dd) < 1e-9 for t, dd in zip(self.T, self.d))
        self.L = [len(t["occ_h"]) for t in self.T]
        self.elig = [[j for j, s in enumerate(sites) if self.tn[k] in s["E"]] for k in range(len(self.T))]
    def lead(self, s):
        return math.ceil(self.sites[s]["delta"] / U - 1e-9)
    def step(self, k, i, s):
        """(值岗开始相对 v 的偏移, 下一节点相对 v 的偏移)"""
        ld = self.lead(s)
        return ld, ld + math.ceil((self.T[k]["occ_h"][i] + self.sites[s]["delta"]) / U - 1e-9)
    def to_json(self):
        return dict(name=self.name, types={n: t for n, t in zip(self.tn, self.T)}, sites=self.sites, cost=self.cost)

# ---------------- 下界 ----------------
def fluid_lb(I, integer=False):
    """流体下界：y[k,s] = 单位时间内 k 型在 s 的班数。覆盖 Σ_k y[k,s]·d_k = K_s；台数 N_k ≥ Σ_s y[k,s]·(平均占用 + δ_s 项)/u（不取整，所以对连续时间模型也成立）。
    integer=True 时 N_k 取整（仍是下界）。"""
    nk, ns = len(I.T), len(I.sites)
    var = [(k, s) for k in range(nk) for s in I.elig[k]]
    n = len(var) + nk
    A = lil_matrix((ns + nk, n)); lo = np.zeros(ns + nk); hi = np.zeros(ns + nk)
    for j, (k, s) in enumerate(var):
        A[s, j] = I.d[k]
        avg = np.mean(I.T[k]["occ_h"]) + 2 * I.sites[s]["delta"]
        A[ns + k, j] = avg / U
    for s in range(ns): lo[s] = hi[s] = I.sites[s]["K"]
    for k in range(nk):
        A[ns + k, len(var) + k] = -1; lo[ns + k] = -np.inf; hi[ns + k] = 0
    c = np.zeros(n); c[len(var):] = [I.cost[t] for t in I.tn]
    integ = np.zeros(n);
    if integer: integ[len(var):] = 1
    r = milp(c, constraints=[LinearConstraint(A.tocsr(), lo, hi)], integrality=integ, bounds=Bounds(0, np.inf))
    return (r.fun, r.x[len(var):]) if r.x is not None else (None, None)

# ---------------- 周期网络 ILP（稳态，无限期） ----------------
def periodic_ilp(I, P, time_limit=120, integer=True, max_num=None):
    """节点 (k, i, v mod P)。工作弧覆盖站点 s 的 d_k 格。目标 Σ_k c_k N_k，N_k = Σ 弧流量×弧长 / P。"""
    nk = len(I.T)
    arcs = []   # (k, i, v, kind, s, dv_next, cov_start, cov_len)
    for k in range(nk):
        for i in range(I.L[k]):
            for v in range(P):
                arcs.append((k, i, v, "z", -1, 1, 0, 0))
                for s in I.elig[k]:
                    ld, nx = I.step(k, i, s)
                    arcs.append((k, i, v, "w", s, nx, ld, I.d[k]))
    node = {}
    for k in range(nk):
        for i in range(I.L[k]):
            for v in range(P): node[(k, i, v)] = len(node)
    na = len(arcs); n = na
    A = lil_matrix((len(node), n))
    for j, (k, i, v, kind, s, dv, cs, cl) in enumerate(arcs):
        A[node[(k, i, v)], j] -= 1
        ni = i if kind == "z" else (i + 1) % I.L[k]
        A[node[(k, ni, (v + dv) % P)], j] += 1
    ns = len(I.sites)
    C = lil_matrix((ns * P, n))
    for j, (k, i, v, kind, s, dv, cs, cl) in enumerate(arcs):
        if kind == "w":
            for q in range(cl): C[s * P + (v + cs + q) % P, j] += 1
    Kvec = np.repeat([s["K"] for s in I.sites], P).astype(float)
    cost = np.array([I.cost[I.tn[a[0]]] * a[5] / P for a in arcs])
    cons = [LinearConstraint(A.tocsr(), 0, 0), LinearConstraint(C.tocsr(), Kvec, Kvec)]
    if max_num:
        M = lil_matrix((nk, n))
        for j, a in enumerate(arcs): M[a[0], j] = a[5] / P
        cons.append(LinearConstraint(M.tocsr(), 0, max_num))
    r = milp(cost, constraints=cons, integrality=np.ones(n) if integer else np.zeros(n), bounds=Bounds(0, np.inf),
             options={"time_limit": time_limit, "mip_rel_gap": 0})
    if r.x is None: return dict(status=r.status, obj=None, N=None, bound=getattr(r, "mip_dual_bound", None))
    Nk = np.zeros(nk)
    for j, a in enumerate(arcs): Nk[a[0]] += r.x[j] * a[5] / P
    return dict(status=r.status, obj=r.fun, N=[round(x, 3) for x in Nk], bound=getattr(r, "mip_dual_bound", None), nvar=n)

# ---------------- 全新启动、有限时域 H 的时间展开 ILP（原问题的设定） ----------------
def fresh_ilp(I, H, time_limit=300, cover_from=None):
    """所有机器人 t=0 时相位 0、满额度；k 型的第一个虚拟时刻 v0_k = ⌈off_k/u⌉。站点 s 从 c0_s 起每格恰好 K_s 台（c0 默认 = ⌈max off/u⌉ + lead_s）。
    越过 H 的弧流向汇点。变量 N_k 为各型投入台数。"""
    nk = len(I.T)
    v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]
    c0 = cover_from if cover_from is not None else max(v0)
    arcs = []
    for k in range(nk):
        for i in range(I.L[k]):
            for v in range(v0[k], H):
                arcs.append((k, i, v, "z", -1, 1, 0, 0))
                for s in I.elig[k]:
                    ld, nx = I.step(k, i, s)
                    if v + ld + I.d[k] <= H or v + ld < H:
                        arcs.append((k, i, v, "w", s, nx, ld, I.d[k]))
    node = {}
    for k in range(nk):
        for i in range(I.L[k]):
            for v in range(v0[k], H): node[(k, i, v)] = len(node)
    na = len(arcs); n = na + nk
    A = lil_matrix((len(node), n))
    for j, (k, i, v, kind, s, dv, cs, cl) in enumerate(arcs):
        A[node[(k, i, v)], j] -= 1
        ni = i if kind == "z" else (i + 1) % I.L[k]
        if v + dv < H: A[node[(k, ni, v + dv)], j] += 1
    for k in range(nk): A[node[(k, 0, v0[k])], na + k] = 1
    ns = len(I.sites)
    rows = {}
    for s in range(ns):
        for t in range(c0 + I.lead(s), H): rows[(s, t)] = len(rows)
    C = lil_matrix((len(rows), n))
    for j, (k, i, v, kind, s, dv, cs, cl) in enumerate(arcs):
        if kind == "w":
            for q in range(cl):
                t = v + cs + q
                if (s, t) in rows: C[rows[(s, t)], j] += 1
    Kvec = np.array([I.sites[s]["K"] for (s, t) in rows], float)
    cost = np.zeros(n); cost[na:] = [I.cost[t] for t in I.tn]
    # 覆盖开始前的工作弧不允许"部分越界"地覆盖未定义的格子：只要求 rows 内恰好 K
    r = milp(cost, constraints=[LinearConstraint(A.tocsr(), 0, 0), LinearConstraint(C.tocsr(), Kvec, Kvec)],
             integrality=np.ones(n), bounds=Bounds(0, np.inf), options={"time_limit": time_limit, "mip_rel_gap": 0})
    if r.x is None: return dict(status=r.status, obj=None, N=None, bound=getattr(r, "mip_dual_bound", None), nvar=n)
    return dict(status=r.status, obj=r.fun, N=[int(round(x)) for x in r.x[na:]], bound=getattr(r, "mip_dual_bound", None), nvar=n)
