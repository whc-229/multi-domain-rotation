"""全新启动（全部机器人相位 0、同时上线）的时间展开下界，机器人 + 共享工具，可选"提前大修"弧。
与 C++ 仿真 rost3/rost4 的规则一致：
- 机器人 k 在虚拟时刻 v 选站点 s：值岗开始 D = v + lead_s，覆盖 [D, D+d_k)，下一虚拟时刻 v + nx（nx 来自步长表），相位 +1；
  提前大修弧：相位 i ≠ L−1 时按相位 L−1 的步长回来，相位回到 0。
- 每格每站恰好 K_s 台在岗，从 c0 + lead_s 开始（c0 = max_k v0_k）。
- 每班需要 q_k 把适配工具，工具 j 被占用 [max(0, D − pre_{k,s}), D + d_k + post)，其中 post 用 q_k·(最小卸载时长) 代替实际卸载总和（更短，保证是下界）。
只要求 [0, H) 内满足，越过 H 的弧流向汇点：任何无限期可行排班截到 H 都满足这些约束，所以 LP（或整数）最优值是下界。"""
import math, sys, json, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
from mtms import U
from tools_model import ToolInst
from export_inst import build

def fresh_lb(I, H, tools=True, reset=False, integer=False, time_limit=600, method="highs-ipm", want_plan=False, full_int=False):
    TI = ToolInst(I) if tools else None
    nk, ns = len(I.T), len(I.sites); nt = TI.nt if tools else 0
    v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]; c0 = max(v0)
    # 变量：机器人弧、工具弧、N_k、M_j
    cols = []   # (kind, data)
    rob_node = {}
    for k in range(nk):
        for i in range(I.L[k]):
            for v in range(v0[k], H): rob_node[(k, i, v)] = len(rob_node)
    tool_node = {}
    for j in range(nt):
        for v in range(H): tool_node[(j, v)] = len(tool_node)
    cov_row = {}
    for s in range(ns):
        for t in range(c0 + I.lead(s), H): cov_row[(s, t)] = len(cov_row)
    duty_key = {}   # (k, s, D) -> 耦合行
    Ar, Ac, Av = [], [], []   # 守恒行：先机器人节点，再工具节点
    Cr, Cc, Cv = [], [], []   # 覆盖行
    Qr, Qc, Qv = [], [], []   # 耦合行：Σ_j 工具弧 − q_k Σ 机器人值岗弧 = 0
    nrn = len(rob_node)
    col = 0; rcols = []; dcols = []   # 提前大修弧：(列, k, i, 值岗开始 D)；所有值岗弧：(列, k, s, D)
    for k in range(nk):
        L = I.L[k]
        for i in range(L):
            for v in range(v0[k], H):
                # 等待
                Ar.append(rob_node[(k, i, v)]); Ac.append(col); Av.append(-1)
                if v + 1 < H: Ar.append(rob_node[(k, i, v + 1)]); Ac.append(col); Av.append(1)
                col += 1
                for s in I.elig[k]:
                    ld, nx = I.step(k, i, s)
                    opts = [(nx, (i + 1) % L)]
                    if reset and i != L - 1: opts.append((I.step(k, L - 1, s)[1], 0))
                    D = v + ld
                    for oi, (dv, ni) in enumerate(opts):
                        if oi == 1: rcols.append((col, k, i, D))
                        dcols.append((col, k, s, D))
                        Ar.append(rob_node[(k, i, v)]); Ac.append(col); Av.append(-1)
                        if v + dv < H: Ar.append(rob_node[(k, ni, v + dv)]); Ac.append(col); Av.append(1)
                        for q in range(I.d[k]):
                            if (s, D + q) in cov_row: Cr.append(cov_row[(s, D + q)]); Cc.append(col); Cv.append(1)
                        if tools and D < H:
                            key = (k, s, D)
                            if key not in duty_key: duty_key[key] = len(duty_key)
                            Qr.append(duty_key[key]); Qc.append(col); Qv.append(-TI.R[k]["q"])
                        col += 1
    n_rob_arcs = col
    if tools:
        for j in range(nt):
            for v in range(H):
                Ar.append(nrn + tool_node[(j, v)]); Ac.append(col); Av.append(-1)
                if v + 1 < H: Ar.append(nrn + tool_node[(j, v + 1)]); Ac.append(col); Av.append(1)
                col += 1
        for (k, s, D), row in duty_key.items():
            R = TI.R[k]; pre = TI.pre(k, s)
            unin_lb = R["q"] * min(TI.TT[jj]["uninstall"] for jj in R["C"])
            for j in R["C"]:
                post = TI.post(k, s, j, unin_lb)
                a, b = max(0, D - pre), D + I.d[k] + post
                Ar.append(nrn + tool_node[(j, a)]); Ac.append(col); Av.append(-1)
                if b < H: Ar.append(nrn + tool_node[(j, b)]); Ac.append(col); Av.append(1)
                Qr.append(row); Qc.append(col); Qv.append(1)
                col += 1
    # 供给变量
    iN = col
    for k in range(nk):
        Ar.append(rob_node[(k, 0, v0[k])]); Ac.append(col); Av.append(1); col += 1
    iM = col
    for j in range(nt):
        Ar.append(nrn + tool_node[(j, 0)]); Ac.append(col); Av.append(1); col += 1
    n = col
    nA = nrn + len(tool_node)
    A = coo_matrix((Av, (Ar, Ac)), shape=(nA, n)).tocsr()
    C = coo_matrix((Cv, (Cr, Cc)), shape=(len(cov_row), n)).tocsr()
    Kv = np.array([I.sites[s]["K"] for (s, t) in cov_row], float)
    cons = [LinearConstraint(A, 0, 0), LinearConstraint(C, Kv, Kv)]
    if tools and duty_key:
        Q = coo_matrix((Qv, (Qr, Qc)), shape=(len(duty_key), n)).tocsr(); cons.append(LinearConstraint(Q, 0, 0))
    c = np.zeros(n); c[iN:] = 1
    integ = np.zeros(n)
    if integer: integ[iN:] = 1
    if full_int: integ[:] = 1; integer = True
    t0 = time.time()
    if integer:
        r = milp(c, constraints=cons, integrality=integ, bounds=Bounds(0, np.inf), options={"time_limit": time_limit, "disp": False})
        if r.x is None: return dict(obj=None, status=r.status, bound=getattr(r, "mip_dual_bound", None), sec=round(time.time() - t0, 1), n=n)
    else:
        from scipy.optimize import linprog
        from scipy.sparse import vstack
        Aeq = vstack([cc.A for cc in cons]).tocsr(); beq = np.concatenate([np.atleast_1d(cc.lb) * np.ones(cc.A.shape[0]) for cc in cons])
        r = linprog(c, A_eq=Aeq, b_eq=beq, bounds=(0, None), method=method, options={"time_limit": time_limit})
    if r.x is None: return dict(obj=None, status=r.status, sec=round(time.time() - t0, 1), n=n)
    out = dict(obj=round(r.fun, 4), N=[round(x, 3) for x in r.x[iN:iM]], M=[round(x, 3) for x in r.x[iM:]],
               status=r.status, sec=round(time.time() - t0, 1), n=n, bound=getattr(r, "mip_dual_bound", None))
    if want_plan:
        from collections import defaultdict
        fl = defaultdict(float)
        for (cc, k, i, D) in rcols:
            if r.x[cc] > 1e-6: fl[(k, i, D)] += r.x[cc]
        out["reset_flow"] = sorted((k, i, D, round(x, 4)) for (k, i, D), x in fl.items())
        dz = defaultdict(float)
        for (cc, k, s, D) in dcols:
            if D < H // 2 and r.x[cc] > 1e-9: dz[(k, s)] += r.x[cc]
        out["duty_flow"] = {f"{k},{s}": round(x, 4) for (k, s), x in dz.items()}
    return out

def round_plan(flows, H_cut=None):
    """按 (k, i) 沿时间累计流量，四舍五入成整数次数：第 t 格的次数 = round(F(t)) − round(F(t 之前))。"""
    from collections import defaultdict
    by = defaultdict(list)
    for k, i, D, x in flows:
        if H_cut is None or D < H_cut: by[(k, i)].append((D, x))
    plan = []
    for (k, i), L in by.items():
        L.sort(); F = 0.0; prev = 0
        for D, x in L:
            F += x; c = int(round(F)) - prev
            if c > 0: plan.append((k, i, D, c)); prev += c
    return sorted(plan, key=lambda e: (e[0], e[1], e[2]))

def S_max(I):
    return max(sum(-(-int(round(o * 10)) // 20) for o in t["occ_h"]) for t in I.T)

if __name__ == "__main__":
    for spec in sys.argv[1].split(","):
        I = build(spec); S = S_max(I)
        mult = float(sys.argv[2]) if len(sys.argv) > 2 else 4
        H = int(mult * S)
        out = dict(spec=spec, inst=I.name, S=S, H=H)
        for reset in [False, True]:
            r = fresh_lb(I, H, tools=True, reset=reset, time_limit=3600)
            out["lp" + ("R" if reset else "")] = r["obj"]; out["sec" + ("R" if reset else "")] = r["sec"]; out["n"] = r["n"]
        print(json.dumps(out, ensure_ascii=False), flush=True)
