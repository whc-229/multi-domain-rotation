"""可证明无限期可行的派工仿真（机器人 + 可选共享工具池），以及"证书适应度"。

状态有限化：机器人的等待特征与工具的空闲特征都截断在 −C 格（C 大于任何提前量），截断后仍可用的机器人/工具彼此不可区分，
派工只依赖截断值，所以截断状态决定未来。于是确定性仿真一定最终进入循环：
  在 t₁ 与 t₂ 两个决策时刻状态完全相同 ⇒ 从 t₁ 起以 t₂ − t₁ 为周期永远重复；
  若 [warm, t₂) 内没有缺岗 ⇒ 严格证明从 warm 起永远不缺岗（无限期可行证书）。
本试点 warm = 0：从第一班起就必须覆盖（与原问题一致，含开局）。"""
import math, json
import numpy as np
from mtms import U, fluid_lb
from tools_model import tool_lb

C = 40   # 截断（格）。所有提前量 pre ≤ 8 格、lead ≤ 2 格，远小于 C。

def run(I, Nk, Mj, w, Tmax, TI=None, warm=0, certify=True):
    nk, ns = len(I.T), len(I.sites)
    v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]
    rk, rph, rv = [], [], []
    for k in range(nk):
        for _ in range(int(Nk[k])): rk.append(k); rph.append(0); rv.append(v0[k])
    nr = len(rk)
    use_tools = TI is not None
    tf = [[-10**6] * int(Mj[j]) for j in range(TI.nt)] if use_tools else []   # 各工具型的空闲时刻列表
    c0 = max(v0)
    lead = [I.lead(s) for s in range(ns)]
    free = [c0 + lead[s] for s in range(ns)]
    st = {(k, i, s): I.step(k, i, s) for k in range(nk) for i in range(I.L[k]) for s in I.elig[k]}
    pre = {(k, s): TI.pre(k, s) for k in range(nk) for s in I.elig[k]} if use_tools else {}
    D = [[0] * I.L[k] for k in range(nk)]; tot = [0] * nk
    flex = [len(I.elig[k]) for k in range(nk)]
    nelig = [sum(1 for k in range(nk) if s in I.elig[k]) for s in range(ns)]
    elig_set = [set(I.elig[k]) for k in range(nk)]
    gaps = 0; why_r = [0] * nk; why_t = [0] * (TI.nt if use_tools else 0)
    seen = {}; gaps_at = {}
    t = min(free)
    while t < Tmax:
        if certify and t >= warm:
            key = (tuple(zip(rk, rph, [max(v - t, -C) for v in rv])), tuple(f - t for f in free),
                   tuple(tuple(tot[k] - I.L[k] * D[k][i] for i in range(I.L[k])) for k in range(nk)),
                   tuple(tuple(sorted(max(f - t, -C) for f in lst)) for lst in tf))
            if key in seen:
                t1 = seen[key]
                return dict(gaps=gaps, certified=(gaps == 0), enter=t1, period=t - t1, t_end=t, why_r=why_r, why_t=why_t)
            seen[key] = t
        due = sorted([s for s in range(ns) if free[s] <= t], key=lambda s: nelig[s])
        for s in due:
            cand = [r for r in range(nr) if s in elig_set[rk[r]] and rv[r] + lead[s] <= t]
            def score(r):
                k, i = rk[r], rph[r]
                return (w[0] * max(rv[r] - t, -C) + w[1] * st[(k, i, s)][1] / I.d[k] + w[2] * flex[k]
                        - w[3] * (tot[k] / I.L[k] - D[k][i]) + w[4] * I.d[k] + w[5] * i / I.L[k], r)
            cand.sort(key=score)
            done = False
            for r in cand:
                k, i = rk[r], rph[r]
                if use_tools:
                    need = t - pre[(k, s)]
                    opts = []
                    for j in TI.R[k]["C"]:
                        for idx, f in enumerate(tf[j]):
                            if f <= need: opts.append((w[6] * TI.TT[j]["maintain"] + w[7] * max(f - t, -C), max(f - t, -C), j, idx))
                    q = TI.R[k]["q"]
                    if len(opts) < q: continue
                    opts.sort()
                    ch = opts[:q]
                    unin = sum(TI.TT[j]["uninstall"] for (_, _, j, _) in ch)
                    for (_, _, j, idx) in ch: tf[j][idx] = t + I.d[k] + TI.post(k, s, j, unin)
                ld, nx = st[(k, i, s)]
                rv[r] = t - ld + nx; rph[r] = (i + 1) % I.L[k]
                D[k][i] += 1; tot[k] += 1
                free[s] = t + I.d[k]; done = True; break
            if not done:
                if t >= warm:
                    gaps += 1
                    if cand:
                        for j in TI.R[rk[cand[0]]]["C"]: why_t[j] += 1
                    else:
                        for k in range(nk):
                            if s in elig_set[k]: why_r[k] += 1
                free[s] = t + 1
        t = min(free)
    return dict(gaps=gaps, certified=None, enter=None, period=None, t_end=t, why_r=why_r, why_t=why_t)

def start_point(I, TI):
    if TI is None:
        _, x = fluid_lb(I, integer=True); return [int(round(v)) for v in x], []
    _, xN, xM = tool_lb(TI, integer=True)
    return [int(round(v)) for v in xN], [int(round(v)) for v in xM]

def size(I, w, H, TI=None, mode="finite", Tcert=20000, max_add=60, start=None):
    """mode='finite'：可行 ⇔ [0,H) 内无缺岗（常规做法）。
    mode='cert'：先按 finite 定台数，再要求无限期证书；证书失败（环内或之前有缺岗）则按瓶颈补资源再证；
                 到 Tcert 仍未出现状态重复且无缺岗，记为"未证"（返回 status='uncert'）。"""
    Nk, Mj = (list(start[0]), list(start[1])) if start else start_point(I, TI)
    nres = len(Nk) + len(Mj)
    def tot(): return sum(Nk) + sum(Mj)
    def bump(idx, d):
        if idx < len(Nk): Nk[idx] += d
        else: Mj[idx - len(Nk)] += d
    def feas_short():
        return run(I, Nk, Mj, w, H, TI, certify=False)
    r = feas_short(); added = 0
    while r["gaps"] > 0 and added < max_add:
        sig = r["why_r"] + r["why_t"]; best = None
        for idx in range(nres):
            bump(idx, 1); g = run(I, Nk, Mj, w, H, TI, certify=False)["gaps"]; bump(idx, -1)
            if best is None or (g, -sig[idx]) < (best[0], -best[2]): best = (g, idx, sig[idx])
        bump(best[1], 1); added += 1; r = feas_short()
    if r["gaps"] > 0: return dict(R=None, Nk=Nk, Mj=Mj, status="fail")
    improved = True
    while improved:
        improved = False
        for idx in range(nres):
            arr = Nk if idx < len(Nk) else Mj; ii = idx if idx < len(Nk) else idx - len(Nk)
            if arr[ii] == 0: continue
            arr[ii] -= 1
            if feas_short()["gaps"] == 0: improved = True; break
            arr[ii] += 1
    if mode == "finite":
        return dict(R=tot(), Nk=list(Nk), Mj=list(Mj), status="finite")
    for _ in range(max_add):
        c = run(I, Nk, Mj, w, Tcert, TI, certify=True)
        if c["certified"]: return dict(R=tot(), Nk=list(Nk), Mj=list(Mj), status="cert", period=c["period"], enter=c["enter"])
        if c["certified"] is None and c["gaps"] == 0:
            return dict(R=tot(), Nk=list(Nk), Mj=list(Mj), status="uncert")
        sig = c["why_r"] + c["why_t"]; bump(int(np.argmax(sig)), 1)
    return dict(R=None, Nk=Nk, Mj=Mj, status="fail")

def es(I, fit_mode, TI=None, H=3000, gens=10, mu=4, lam=12, seed=0, Tcert=20000):
    """(μ+λ)-ES。适应度 = 资源数；cert 模式下"未证"罚 0.5，失败罚 999。"""
    rng = np.random.default_rng(seed); dim = 8 if TI is not None else 6
    base = [np.array([1, 0, 0, 0, 0, 0, 1, 0.][:dim]), np.array([0.1, 0, 10, 1, 0, 0, 1, 0.1][:dim]), np.array([1, 0.5, 1, 0.5, 0, 0, 1, -0.1][:dim])]
    st0 = start_point(I, TI)
    def fit(wv):
        r = size(I, wv, H, TI, fit_mode, Tcert, start=st0)
        f = 999 if r["R"] is None else r["R"] + (0.5 if r["status"] == "uncert" else 0)
        return f, r
    pop = base + [rng.normal(0, 1, dim) for _ in range(mu + lam - len(base))]
    sc = sorted([(fit(wv), wv) for wv in pop], key=lambda z: z[0][0]); par = sc[:mu]; hist = [par[0][0][0]]
    for g in range(gens):
        kids = []
        for _ in range(lam):
            a = par[rng.integers(len(par))][1]
            kids.append(a + 0.5 * rng.normal(0, 1, dim) * (np.abs(a) + 0.2))
        par = sorted(par + [(fit(wv), wv) for wv in kids], key=lambda z: z[0][0])[:mu]; hist.append(par[0][0][0])
    (f, r), wv = par[0]
    return f, r, wv.tolist(), hist
