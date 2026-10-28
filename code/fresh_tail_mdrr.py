"""T4 上界：把 超长周期_数学/atlas/fresh_tail.py 推广到多站、多机型、带共享工具。

模型（整数，highspy MIP）：在 [0, Hx)，Hx = T + P + G 上做全新启动时间展开：
- 机器人：节点 (k, i, v)，v ∈ [v0_k, Hx)；等待弧、值岗弧（站 s，下一节点 v + nx；只在 v ≥ c0 时有，与仿真的开班时刻一致），可选提前大修弧（相位 i ≠ L−1 时按相位 L−1 的步长回到相位 0）。
  弧头 ≥ Hx 的弧流向汇点。N_k 在 (k, 0, v0_k) 注入。
- 覆盖：(s, t)，t ∈ [c0 + lead_s, Hx)，恰好 K_s。
- 工具（与仿真规则一致，不是下界里的近似）：每班从 C_k 里选 q_k 把（可重复的组合 c），实际卸载总时长 unin_c 决定每把工具的占用
  [D − pre_{k,s}, D + d_k + post(k, s, j, unin_c))。变量 z_{k,s,D,c} = 在 (k, s, D) 用组合 c 的班数：
  Σ_c z_{k,s,D,c} = 机器人在 (k, s, D) 的值岗弧流量；z 驱动工具节点 (j, v) 上的流（每用一把 j 型工具，从起点流到终点）。M_j 在 (j, 0) 注入。
- 周期：每个变量有时间下标（机器人弧 = 出发 v，z = 值岗开始 D，工具等待弧 = v）。下标 τ ∈ [T, T + G) 的变量等于下标 τ + P 的同类变量。
  G = (各约束里变量下标相对约束时刻向后的最大跨度) + (向前的最大跨度) + 1。于是窗口外的每条约束都是窗口内某条完整约束的平移，
  解按周期 P 无限延拓后处处可行（证明见 报告/T4_小实例精确.md）。所以最优值是无限期全新启动问题的严格上界。
输出：目标值、状态、各型台数、非零的值岗弧流和 z，用于独立重放。
用法：python3 code/fresh_tail_mdrr.py <实例名> <T> <P> <reset 0/1> <限时秒> <输出 json>"""
import os, sys, json, math, time, itertools
import numpy as np
from scipy.sparse import coo_matrix
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import U, ToolInst, load_train  # noqa: E402


def combos(TI, k):
    """k 型每班可选的工具组合（多重集）及每个组合的卸载总时长。"""
    R = TI.R[k]
    out = []
    for c in itertools.combinations_with_replacement(sorted(R["C"]), R["q"]):
        cnt = {}
        for j in c: cnt[j] = cnt.get(j, 0) + 1
        out.append((tuple(sorted(cnt.items())), sum(TI.TT[j]["uninstall"] for j in c)))
    return out


def spans(I, TI, reset):
    """返回 (向后跨度, 向前跨度)。"""
    nk = len(I.T); back = 1; fwd = 0
    for k in range(nk):
        for s in I.elig[k]:
            ld = I.lead(s)
            for i in range(I.L[k]): back = max(back, I.step(k, i, s)[1])
            back = max(back, ld + I.d[k])
            fwd = max(fwd, TI.pre(k, s))
            for (cnt, unin) in combos(TI, k):
                for j, _ in cnt: back = max(back, I.d[k] + TI.post(k, s, j, unin))
    return back, fwd


def build(I, T, P, reset):
    """建模；返回 (Highs 模型数据, 变量说明)。"""
    TI = ToolInst(I)
    nk, ns, nt = len(I.T), len(I.sites), TI.nt
    back, fwd = spans(I, TI, reset)
    G = back + fwd + 1
    Hx = T + P + G
    v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]; c0 = max(v0)
    lead = [I.lead(s) for s in range(ns)]
    rows = {}
    def row(key):
        r = rows.get(key)
        if r is None: r = rows[key] = len(rows)
        return r
    R, C, V = [], [], []; cost = []; integ = []; meta = []
    tindex = []     # 变量的 (类别键, 时间下标)，用于周期约束
    def add_col(entries, c, m, cls, tau):
        j = len(cost)
        for (key, v) in entries: R.append(row(key)); C.append(j); V.append(v)
        cost.append(c); meta.append(m); tindex.append((cls, tau)); return j
    # 先登记等式行（覆盖行 rhs = K）
    for s in range(ns):
        for t in range(c0 + lead[s], Hx): row(("cov", s, t))
    # 机器人弧
    for k in range(nk):
        L = I.L[k]
        for i in range(L):
            for v in range(v0[k], Hx):
                ent = [(("rn", k, i, v), -1.0)]
                if v + 1 < Hx: ent.append((("rn", k, i, v + 1), 1.0))
                add_col(ent, 0.0, ("w", k, i, v), ("w", k, i), v)
                if v < c0: continue      # 与仿真一致：岗位从 c0 + lead_s 起才开班，值岗不能提前（freshlb 的下界允许，只会更松）
                for s in I.elig[k]:
                    ld, nx = I.step(k, i, s); D = v + ld
                    opts = [(nx, (i + 1) % L, 0)]
                    if reset and i != L - 1: opts.append((I.step(k, L - 1, s)[1], 0, 1))
                    for (dv, ni, r) in opts:
                        ent = [(("rn", k, i, v), -1.0)]
                        if v + dv < Hx: ent.append((("rn", k, ni, v + dv), 1.0))
                        for q in range(I.d[k]):
                            if c0 + lead[s] <= D + q < Hx: ent.append((("cov", s, D + q), 1.0))
                        if D < Hx: ent.append((("cp", k, s, D), -1.0))
                        add_col(ent, 0.0, ("d", k, i, v, s, r), ("d", k, i, s, r), v)
    # 工具组合变量 z
    for k in range(nk):
        cb = combos(TI, k)
        for s in I.elig[k]:
            pre = TI.pre(k, s)
            for D in range(c0 + lead[s], Hx):
                for ci, (cnt, unin) in enumerate(cb):
                    ent = [(("cp", k, s, D), 1.0)]
                    a = max(0, D - pre)
                    for j, n in cnt:
                        b = D + I.d[k] + TI.post(k, s, j, unin)
                        ent.append((("tn", j, a), -float(n)))
                        if b < Hx: ent.append((("tn", j, b), float(n)))
                    add_col(ent, 0.0, ("z", k, s, D, ci), ("z", k, s, ci), D)
    # 工具等待弧
    for j in range(nt):
        for v in range(Hx):
            ent = [(("tn", j, v), -1.0)]
            if v + 1 < Hx: ent.append((("tn", j, v + 1), 1.0))
            add_col(ent, 0.0, ("tw", j, v), ("tw", j), v)
    iN = len(cost)
    for k in range(nk): add_col([(("rn", k, 0, v0[k]), 1.0)], I.cost[I.tn[k]], ("N", k), ("N", k), -1)
    iM = len(cost)
    for j in range(nt): add_col([(("tn", j, 0), 1.0)], 1.0, ("M", j), ("M", j), -1)
    # 周期约束：下标 τ ∈ [T, T+G) 的变量 = 下标 τ+P 的同类变量
    pos = {ti: j for j, ti in enumerate(tindex) if ti[1] >= 0}
    nper = 0
    for (cls, tau), j in list(pos.items()):
        if T <= tau < T + G:
            j2 = pos.get((cls, tau + P))
            if j2 is None: continue
            rr = row(("per", nper)); nper += 1
            R.extend([rr, rr]); C.extend([j, j2]); V.extend([1.0, -1.0])
    ncol = len(cost); nrow = len(rows)
    lo = np.zeros(nrow)
    for key, r in rows.items():
        if key[0] == "cov": lo[r] = I.sites[key[1]]["K"]
    A = coo_matrix((V, (R, C)), shape=(nrow, ncol)).tocsc()
    return dict(A=A, rhs=lo, cost=np.array(cost), meta=meta, iN=iN, iM=iM, Hx=Hx, G=G, back=back, fwd=fwd, nper=nper, v0=v0, c0=c0)


def solve(I, T, P, reset, time_limit, integer=True, opts=None):
    """求解；返回结果 dict（含非零流）。"""
    import highspy      # 放在函数内：code/t4_cpsat.py 只用 build()，而 OR-Tools 与 highspy 不能在同一进程里同时加载
    t0 = time.time()
    M = build(I, T, P, reset)
    A = M["A"]; n, m = A.shape[1], A.shape[0]
    h = highspy.Highs(); h.setOptionValue("output_flag", False); h.setOptionValue("threads", 1)
    h.setOptionValue("time_limit", float(time_limit)); h.setOptionValue("mip_rel_gap", 0.0)
    for key, val in (opts or {}).items(): h.setOptionValue(key, val)
    lp = highspy.HighsLp()
    lp.num_col_ = n; lp.num_row_ = m
    lp.col_cost_ = M["cost"]; lp.col_lower_ = np.zeros(n); lp.col_upper_ = np.full(n, highspy.kHighsInf)
    lp.row_lower_ = M["rhs"]; lp.row_upper_ = M["rhs"].copy()
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A.indptr; lp.a_matrix_.index_ = A.indices; lp.a_matrix_.value_ = A.data
    lp.a_matrix_.num_col_ = n; lp.a_matrix_.num_row_ = m
    if integer: lp.integrality_ = [highspy.HighsVarType.kInteger] * n
    h.passModel(lp)
    tb = time.time() - t0
    h.run()
    st = h.modelStatusToString(h.getModelStatus()); info = h.getInfo()
    out = dict(T=T, P=P, reset=int(reset), Hx=M["Hx"], G=M["G"], back=M["back"], fwd=M["fwd"], ncol=n, nrow=m, status=st,
               build_sec=round(tb, 1), sec=round(time.time() - t0, 1))
    if info.primal_solution_status >= 2 or st == "Optimal":
        x = np.array(h.getSolution().col_value)
        out["obj"] = float(round(info.objective_function_value, 6))
        out["bound"] = float(info.mip_dual_bound) if integer else None
        out["N"] = [int(round(v)) for v in x[M["iN"]:M["iM"]]]
        out["M"] = [int(round(v)) for v in x[M["iM"]:]]
        nz = [(M["meta"][j], int(round(x[j]))) for j in range(n) if x[j] > 0.5 and M["meta"][j][0] in ("d", "z")]
        out["duties"] = [list(mm[1:]) + [c] for mm, c in nz if mm[0] == "d"]
        out["z"] = [list(mm[1:]) + [c] for mm, c in nz if mm[0] == "z"]
        out["frac"] = float(np.max(np.abs(x - np.round(x)))) if integer else None
    else:
        out["obj"] = None; out["bound"] = float(info.mip_dual_bound) if integer else None
    return out


def main():
    name, T, P, rs, tl, outp = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), sys.argv[6]
    if os.path.exists(outp): print("已完成，跳过", outp); return
    I = load_train(name)
    r = solve(I, T, P, bool(rs), tl)
    r["inst"] = name
    tmp = outp + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(r) + "\n")
    os.replace(tmp, outp)
    print(json.dumps({k: v for k, v in r.items() if k not in ("duties", "z")}))


if __name__ == "__main__":
    main()
