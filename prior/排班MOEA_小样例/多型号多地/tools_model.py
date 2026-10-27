"""第二阶段：加入基地共享工具池（按原 roster_sim 规则）。
- 工具型 j：安装、卸载、保养时长（小时）。机器人型 k：适配工具型集合 C_k、每班需 q_k 把、go/back/post、max_install。
- 一班（k 型、站点 s、值岗开始 D）占用所选工具的区间：
    需在 D − ⌈(go + δ_s + max_install)/u⌉ 之前空闲；占用到 D + d_k + ⌈(post + back + δ_s + Σ卸载 + 保养_j)/u⌉。
  向上取整，保守。工具没有维修超周期，但被所有站点共享，所以把各站耦合在一起。
- 资源数 = 机器人数 + 工具数（与原问题 f1 一致）。"""
import math
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
from roster_sim import ROBOTS, TOOLS, Instance
from mtms import real_types, Inst, U

def tool_params():
    return [dict(install=t["install"], uninstall=t["uninstall"], maintain=t["maintain"]) for t in TOOLS]

def robot_tool_params():
    """机型名 → 适配工具与时长参数（同一序号的机型在 L1/L2/L3 中工具参数相同，只有维修参数不同）。"""
    out = {}
    inst = Instance("x", 720, 3)
    for lv in ["L1", "L2", "L3"]:
        for r, p in enumerate(inst.robots):
            out[f"{lv}-{r}"] = dict(C=list(p["tools"]), q=p["tool_num"], go=p["go"], back=p["back"], post=p["post_consume"],
                                    max_install=p["max_install"])
    return out

class ToolInst:
    def __init__(self, I):
        self.I = I; self.TT = tool_params(); RP = robot_tool_params()
        self.R = [RP[n] for n in I.tn]
        self.nt = len(self.TT)
    def pre(self, k, s):
        r = self.R[k]; return math.ceil((r["go"] + self.I.sites[s]["delta"] + r["max_install"]) / U - 1e-9)
    def post(self, k, s, j, unin_sum):
        r = self.R[k]; return math.ceil((r["post"] + r["back"] + self.I.sites[s]["delta"] + unin_sum + self.TT[j]["maintain"]) / U - 1e-9)

def tool_lb(TI, integer=True):
    """机器人 + 工具的流体下界（LP 或整数台数）。"""
    I = TI.I; nk, ns, nt = len(I.T), len(I.sites), TI.nt
    yv = [(k, s) for k in range(nk) for s in I.elig[k]]
    zv = [(k, s, j) for (k, s) in yv for j in TI.R[k]["C"]]
    n = len(yv) + len(zv) + nk + nt
    iy = {v: a for a, v in enumerate(yv)}; iz = {v: len(yv) + a for a, v in enumerate(zv)}
    iN = len(yv) + len(zv); iM = iN + nk
    rows = []; A = lil_matrix((ns + len(yv) + nk + nt, n)); lo = []; hi = []
    r = 0
    for s in range(ns):
        for (k, ss) in yv:
            if ss == s: A[r, iy[(k, ss)]] = I.d[k]
        lo.append(I.sites[s]["K"]); hi.append(I.sites[s]["K"]); r += 1
    for (k, s) in yv:   # Σ_j z = q_k y
        for j in TI.R[k]["C"]: A[r, iz[(k, s, j)]] = 1
        A[r, iy[(k, s)]] = -TI.R[k]["q"]; lo.append(0); hi.append(0); r += 1
    for k in range(nk):
        for (kk, s) in yv:
            if kk == k: A[r, iy[(kk, s)]] = (np.mean(I.T[k]["occ_h"]) + 2 * I.sites[s]["delta"]) / U
        A[r, iN + k] = -1; lo.append(-np.inf); hi.append(0); r += 1
    for j in range(nt):
        for (k, s, jj) in zv:
            if jj == j:
                rk = TI.R[k]
                busy = (rk["max_install"] * 0 + rk["go"] + 2 * I.sites[s]["delta"] + I.T[k]["duty_h"] + rk["post"] + rk["back"]
                        + TI.TT[j]["install"] + TI.TT[j]["uninstall"] + TI.TT[j]["maintain"]) / U   # 不取整、只计本工具自己的装卸，保证是下界
                A[r, iz[(k, s, j)]] = busy
        A[r, iM + j] = -1; lo.append(-np.inf); hi.append(0); r += 1
    c = np.zeros(n); c[iN:] = 1
    integ = np.zeros(n)
    if integer: integ[iN:] = 1
    res = milp(c, constraints=[LinearConstraint(A.tocsr(), lo, hi)], integrality=integ, bounds=Bounds(0, np.inf))
    return res.fun, res.x[iN:iM], res.x[iM:]

def simulate_tools(TI, Nk, Mj, policy, H, warm, params=None, fresh=True):
    """机器人 + 工具联合派工仿真。返回 (缺岗格数, 各站缺岗)。工具选择：适配且按时空闲的工具中，按策略选 q 把。"""
    I = TI.I; nk, ns = len(I.T), len(I.sites)
    v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]
    robots = [[k, 0, v0[k] if fresh else 0] for k in range(nk) for _ in range(int(Nk[k]))]
    tools = [[j, -10**9] for j in range(TI.nt) for _ in range(int(Mj[j]))]
    c0 = max(v0) if fresh else 0
    lead = [I.lead(s) for s in range(ns)]
    free = [c0 + lead[s] for s in range(ns)]
    st = {(k, i, s): I.step(k, i, s) for k in range(nk) for i in range(I.L[k]) for s in I.elig[k]}
    pre = {(k, s): TI.pre(k, s) for k in range(nk) for s in I.elig[k]}
    D = [[0] * I.L[k] for k in range(nk)]; tot = [0] * nk
    flex = [len(I.elig[k]) for k in range(nk)]
    nelig = [sum(1 for k in range(nk) if s in I.elig[k]) for s in range(ns)]
    gaps = 0; gap_sites = [0] * ns
    w = params
    why = {'robot': [0] * nk, 'tool': [0] * TI.nt}
    t = min(free)
    while t < H:
        due = sorted([s for s in range(ns) if free[s] <= t], key=lambda s: nelig[s])
        for s in due:
            cand = [r for r in robots if s in I.elig[r[0]] and r[2] + lead[s] <= t]
            if policy == "fifo":
                cand.sort(key=lambda r: (r[2], r[0]))
            elif policy == "pb":
                cand.sort(key=lambda r: (flex[r[0]], -(tot[r[0]] / I.L[r[0]] - D[r[0]][r[1]]), r[2]))
            else:   # score
                def key(r):
                    k, i = r[0], r[1]; nxt = st[(k, i, s)][1]
                    return (w[0] * (r[2] - t) + w[1] * nxt / I.d[k] + w[2] * flex[k] - w[3] * (tot[k] / I.L[k] - D[k][i])
                            + w[4] * I.d[k] + w[5] * i / I.L[k])
                cand.sort(key=key)
            done = False
            for r in cand:
                k, i = r[0], r[1]
                need_free = t - pre[(k, s)]
                tc = [tl for tl in tools if tl[0] in TI.R[k]["C"] and tl[1] <= need_free]
                if len(tc) < TI.R[k]["q"]: continue
                if policy == "score" and len(w) > 6:
                    tc.sort(key=lambda tl: (w[6] * TI.TT[tl[0]]["maintain"] + w[7] * (need_free - tl[1])))
                else:
                    tc.sort(key=lambda tl: (TI.TT[tl[0]]["maintain"], tl[1]))   # 保养短、空闲早者优先
                ch = tc[:TI.R[k]["q"]]
                unin = sum(TI.TT[tl[0]]["uninstall"] for tl in ch)
                for tl in ch: tl[1] = t + I.d[k] + TI.post(k, s, tl[0], unin)
                ld, nx = st[(k, i, s)]
                r[2] = t - ld + nx; r[1] = (i + 1) % I.L[k]
                D[k][i] += 1; tot[k] += 1
                free[s] = t + I.d[k]; done = True; break
            if not done:
                if t >= warm:
                    gaps += 1; gap_sites[s] += 1
                    if cand:
                        for r in cand[:1]:
                            for j in TI.R[r[0]]["C"]: why['tool'][j] += 1
                    else:
                        for k in range(nk):
                            if s in I.elig[k]: why['robot'][k] += 1
                free[s] = t + 1
        t = min(free)
    simulate_tools.why = why
    return gaps, gap_sites

def size_resources(TI, policy, H, warm, params=None, fresh=True, max_add=60):
    _, xN, xM = tool_lb(TI, integer=True)
    Nk = [int(round(v)) for v in xN]; Mj = [int(round(v)) for v in xM]
    g, _ = simulate_tools(TI, Nk, Mj, policy, H, warm, params, fresh)
    added = 0
    while g > 0 and added < max_add:
        simulate_tools(TI, Nk, Mj, policy, H, warm, params, fresh); why0 = dict(simulate_tools.why)
        best = None
        for which in range(len(Nk) + len(Mj)):
            if which < len(Nk): Nk[which] += 1
            else: Mj[which - len(Nk)] += 1
            gg, _ = simulate_tools(TI, Nk, Mj, policy, H, warm, params, fresh)
            if which < len(Nk): Nk[which] -= 1
            else: Mj[which - len(Nk)] -= 1
            # 缺岗数相同则按瓶颈信号打破平局：该资源在缺岗原因中出现得越多越优先
            sig = why0['robot'][which] if which < len(Nk) else why0['tool'][which - len(Nk)]
            if best is None or (gg, -sig) < (best[0], -best[2]): best = (gg, which, sig)
        if best[1] < len(Nk): Nk[best[1]] += 1
        else: Mj[best[1] - len(Nk)] += 1
        g = best[0]; added += 1
    improved = True
    while improved and g == 0:
        improved = False
        for which in range(len(Nk) + len(Mj)):
            arr, idx = (Nk, which) if which < len(Nk) else (Mj, which - len(Nk))
            if arr[idx] == 0: continue
            arr[idx] -= 1
            gg, _ = simulate_tools(TI, Nk, Mj, policy, H, warm, params, fresh)
            if gg == 0: improved = True; break
            arr[idx] += 1
    return (sum(Nk) + sum(Mj) if g == 0 else None), Nk, Mj
