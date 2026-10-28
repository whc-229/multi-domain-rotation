"""T6 单个作业：下界口径修正。结果写进 results/T6/jobs/<名>.json（原子写入，已存在就跳过）。

种类：
- rate：网格口径的速率 LP（稳态，无限期），与 freshlb 的步长表、工具占用口径一致：
  * 变量 y[k,i,s]：k 型在相位 i 去 s 站值岗的班率（每格班数）；允许稳态提前大修时另有 z[k,i,s]（值完这班回到相位 0，i ≠ L−1）；
  * 覆盖：Σ_k d_k Σ_i (y + z)[k,i,s] = K_s；
  * 相位流守恒：每型每相位 流入 = 流出（正常值岗 i → i+1，提前大修 i → 0）；
  * 机器人：N_k ≥ Σ y·nx(k,i,s) + Σ z·nx(k,L−1,s)，nx 取 I.step 的步长（虚拟时刻到下一虚拟时刻，含 lead）；
  * 工具：Σ_j x[k,s,j] = q_k Σ_i (y + z)[k,i,s]；M_j ≥ Σ x[k,s,j]·(pre + d_k + post)，pre / post 与 freshlb 相同（卸载按 q_k·最短卸载）。
  稳态排班里每台机器人的时间被它的各班步长（加等待）铺满，所以 N_k ≥ Σ 班率 × 步长；这是周期排班的有效下界。
  输出 R_rate_g（不提前大修）、R_rate_g_oldtools（网格机器人 + tool_lb 的工具口径，只用于拆分"网格"与"工具口径"两部分差异）、R_rate_g^R（允许稳态提前大修，只作参照）及各自的整数台数版本，
  每型每站的平均周期 ā（每班格数，不提前大修）与最小平均周期 ā*（允许提前大修），以及 LP 解的班率 y[k,s]（fresh_rlp_w 外推用）。
- frRw：全新启动 LP_fr^Rw(H)：freshlb 的模型（允许提前大修），但提前大修弧只在值岗开始 D < W = 2·S_max 时保留（删列）。
  求解、索引核对、对偶证书、保存对偶与值岗流都与 T2 的 frR 相同；证书的定价也只在窗口内允许提前大修（code/cpp/price.cpp 的 longest_path_w）。
用法：python3 code/t6_job.py <实例名> rate
      python3 code/t6_job.py <实例名> frRw <H 倍数>"""
import os, sys, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, U, ToolInst, load_train  # noqa: E402
import highspy  # noqa: E402
from freshlb import S_max  # noqa: E402

OUT = os.environ.get("T6_OUT", os.path.join(ROOT, "results", "T6", "jobs"))       # 环境变量只用于开发期试算
DUAL = os.environ.get("T6_DUAL", os.path.join(ROOT, "results", "T6_duals"))


def solve_lp(c, rows, lo, hi, n, integ=None):
    """min c·x，lo ≤ A x ≤ hi，x ≥ 0；rows 为 [(行号, 列号, 值)]。integ 为整数列下标。返回 (状态, 目标值, x)。"""
    from scipy.sparse import coo_matrix
    R, C, V = zip(*rows)
    A = coo_matrix((V, (R, C)), shape=(len(lo), n)).tocsc()
    h = highspy.Highs(); h.setOptionValue("output_flag", False); h.setOptionValue("threads", 1)
    h.setOptionValue("mip_rel_gap", 0.0)
    lp = highspy.HighsLp(); lp.num_col_ = n; lp.num_row_ = len(lo)
    lp.col_cost_ = np.asarray(c, float); lp.col_lower_ = np.zeros(n); lp.col_upper_ = np.full(n, highspy.kHighsInf)
    lp.row_lower_ = np.asarray(lo, float); lp.row_upper_ = np.asarray(hi, float)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A.indptr; lp.a_matrix_.index_ = A.indices; lp.a_matrix_.value_ = A.data
    lp.a_matrix_.num_col_ = n; lp.a_matrix_.num_row_ = len(lo)
    if integ is not None:
        t = [highspy.HighsVarType.kContinuous] * n
        for a in integ: t[a] = highspy.HighsVarType.kInteger
        lp.integrality_ = t
    h.passModel(lp); h.run()
    st = h.modelStatusToString(h.getModelStatus())
    if st != "Optimal": return st, None, None
    return st, h.getInfo().objective_function_value, np.array(h.getSolution().col_value)


def rate_g(I, reset, integer=False, tools="fresh"):
    """网格口径速率 LP（见文件头）。返回 dict(obj, N, M, y_ks, z_flow)。
    tools = "fresh"：工具占用与 freshlb 相同（主口径）；"old"：工具占用沿用 tool_lb（小时、不取整、只计本工具装卸），只用来拆分口径差异。"""
    TI = ToolInst(I); nk, ns, nt = len(I.T), len(I.sites), TI.nt
    col = {}; c = []
    def var(key, cost=0.0):
        col[key] = len(c); c.append(cost); return col[key]
    for k in range(nk):
        L = I.L[k]
        for i in range(L):
            for s in I.elig[k]:
                var(("y", k, i, s))
                if reset and i != L - 1: var(("z", k, i, s))
        for s in I.elig[k]:
            for j in TI.R[k]["C"]: var(("x", k, s, j))
    for k in range(nk): var(("N", k), 1.0)
    for j in range(nt): var(("M", j), 1.0)
    rows, lo, hi = [], [], []
    def row(terms, a, b):
        r = len(lo); lo.append(a); hi.append(b)
        for key, v in terms: rows.append((r, col[key], v))
    shift = lambda k, i, s: [("y", k, i, s)] + ([("z", k, i, s)] if reset and i != I.L[k] - 1 else [])
    for s in range(ns):      # 覆盖
        row([(key, I.d[k]) for k in range(nk) if s in I.elig[k] for i in range(I.L[k]) for key in shift(k, i, s)], I.sites[s]["K"], I.sites[s]["K"])
    for k in range(nk):      # 相位流守恒：流出 − 流入 = 0
        L = I.L[k]
        for i in range(L):
            t = [(key, 1.0) for s in I.elig[k] for key in shift(k, i, s)]
            if i == 0:
                t += [(("y", k, L - 1, s), -1.0) for s in I.elig[k]]
                if reset: t += [(("z", k, i2, s), -1.0) for i2 in range(L - 1) for s in I.elig[k]]
            else:
                t += [(("y", k, i - 1, s), -1.0) for s in I.elig[k]]
            if L == 1: t = []      # 单相位：流出 = 流入恒成立
            if t: row(t, 0, 0)
        # 机器人台数
        t = [(("N", k), -1.0)]
        for i in range(L):
            for s in I.elig[k]:
                t.append((("y", k, i, s), I.step(k, i, s)[1]))
                if reset and i != L - 1: t.append((("z", k, i, s), I.step(k, L - 1, s)[1]))
        row(t, -np.inf, 0)
        # 工具耦合
        for s in I.elig[k]:
            row([(("x", k, s, j), 1.0) for j in TI.R[k]["C"]] + [(key, -TI.R[k]["q"]) for i in range(L) for key in shift(k, i, s)], 0, 0)
    for j in range(nt):      # 工具台数
        t = [(("M", j), -1.0)]
        for k in range(nk):
            if j not in TI.R[k]["C"]: continue
            Rk = TI.R[k]; unin_lb = Rk["q"] * min(TI.TT[jj]["uninstall"] for jj in Rk["C"])
            for s in I.elig[k]:
                if tools == "fresh": occ = TI.pre(k, s) + I.d[k] + TI.post(k, s, j, unin_lb)
                else: occ = (Rk["go"] + 2 * I.sites[s]["delta"] + I.T[k]["duty_h"] + Rk["post"] + Rk["back"]
                             + TI.TT[j]["install"] + TI.TT[j]["uninstall"] + TI.TT[j]["maintain"]) / U
                t.append((("x", k, s, j), occ))
        row(t, -np.inf, 0)
    integ = [col[("N", k)] for k in range(nk)] + [col[("M", j)] for j in range(nt)] if integer else None
    st, obj, x = solve_lp(c, rows, lo, hi, len(c), integ)
    if obj is None: return dict(status=st, obj=None)
    yks = {}; zf = {}
    for key, a in col.items():
        if key[0] in ("y", "z") and x[a] > 1e-12:
            k, i, s = key[1:]; yks[f"{k},{s}"] = yks.get(f"{k},{s}", 0.0) + float(x[a])
            if key[0] == "z": zf[f"{k},{i},{s}"] = float(x[a])
    return dict(status=st, obj=float(obj), N=[float(x[col[("N", k)]]) for k in range(nk)], M=[float(x[col[("M", j)]]) for j in range(nt)],
                y=yks, z=zf)


def cycle_table(I):
    """每型每站：不提前大修的平均周期 ā = mean_i nx(k,i,s)，与允许提前大修的最小平均周期 ā*（单型单站的小 LP：
    min Σ (y_i nx_i + z_i nx_{L−1})，Σ (y + z) = 1，相位流守恒）。单位：格 / 班。"""
    out = []
    for k in range(len(I.T)):
        L = I.L[k]
        for s in I.elig[k]:
            nx = [I.step(k, i, s)[1] for i in range(L)]; nxR = I.step(k, L - 1, s)[1]
            col = {}; c = []
            for i in range(L): col[("y", i)] = len(c); c.append(nx[i])
            for i in range(L - 1): col[("z", i)] = len(c); c.append(nxR)
            rows, lo, hi = [], [], []
            r = 0; lo.append(1); hi.append(1)
            for a in range(len(c)): rows.append((r, a, 1.0))
            for i in range(L):
                if L == 1: break
                r = len(lo); lo.append(0); hi.append(0)
                rows.append((r, col[("y", i)], 1.0))
                if i < L - 1: rows.append((r, col[("z", i)], 1.0))
                if i == 0:
                    rows.append((r, col[("y", L - 1)], -1.0))
                    for i2 in range(L - 1): rows.append((r, col[("z", i2)], -1.0))
                else:
                    rows.append((r, col[("y", i - 1)], -1.0))
            _, obj, _ = solve_lp(c, rows, lo, hi, len(c))
            out.append(dict(k=k, type=I.tn[k], s=s, delta=I.sites[s]["delta"], abar=float(np.mean(nx)), astar=float(obj)))
    return out


def rate_job(I):
    """rate 种类：四个速率下界（LP / 整数，不提前大修 / 允许稳态提前大修）+ 周期表。"""
    t0 = time.time(); out = {}
    for reset, tag in ((False, "g"), (True, "gR")):
        a = rate_g(I, reset); b = rate_g(I, reset, integer=True)
        out[f"R_rate_{tag}"] = a["obj"]; out[f"R_rate_{tag}_int"] = b["obj"]
        out[f"N_{tag}"] = a.get("N"); out[f"M_{tag}"] = a.get("M"); out[f"y_{tag}"] = a.get("y")
        out[f"N_{tag}_int"] = b.get("N"); out[f"M_{tag}_int"] = b.get("M")
        if reset: out["z_gR"] = a.get("z")
    a = rate_g(I, False, tools="old"); out["R_rate_g_oldtools"] = a["obj"]     # 拆分用：网格机器人 + 旧工具口径
    out["cycles"] = cycle_table(I)
    out["has_L3_2"] = int("L3-2" in I.tn)
    out["sec"] = round(time.time() - t0, 2)
    return out


def frRw_job(I, H, W, dual_path):
    """LP_fr^Rw：删去 D ≥ W 的提前大修弧后求解，其余与 t2_job.fresh_job 相同（索引核对、对偶证书、保存流量）。"""
    import t2_job
    from fresh_full import fresh_full
    idx = t2_job.fresh_index(I, H, True)
    drop = [col for (col, k, i, v, s, r) in idx["arcs"] if r == 1 and v + I.lead(s) >= W]
    keep = {}
    t0 = time.time()
    r = fresh_full(I, H, reset=True, solver="ipm", time_limit=6 * 3600, keep=keep, drop=drop)
    out = dict(obj=r["obj"], N=r.get("N"), M=r.get("M"), lp_status=keep.get("status"), n=r.get("n"), n_drop=len(drop), sec=round(time.time() - t0, 1))
    if r["obj"] is None or keep.get("row_dual") is None: return out
    x = keep["x"]; obj = float(keep["c"] @ x); out["obj_full"] = obj
    out["drop_max_x"] = float(np.max(np.abs(x[drop]))) if drop else 0.0
    bad, rows_ok = t2_job.verify_index(I, H, idx, keep["A"])
    out["index_mismatch"] = bad; out["rows_ok"] = rows_ok
    # 对偶证书：定价只在窗口内允许提前大修
    import cg_fresh
    old = cg_fresh.FreshCG
    cg_fresh.FreshCG = lambda I_, H_, reset_: old(I_, H_, reset_, W=W)
    try:
        t1 = time.time(); cert = t2_job.certificate(I, H, True, idx, keep["row_dual"], obj)
    finally:
        cg_fresh.FreshCG = old
    out.update({k: cert[k] for k in ("path_viol", "path_viol_k", "tool_viol", "wait_viol", "M_viol", "dual_obj", "dual_gap")})
    out["cert_sec"] = round(time.time() - t1, 1)
    fl = {}
    for (col, k, i, v, s, rr) in idx["arcs"]:
        if s < 0 or x[col] <= 1e-9: continue
        key = (k, i, s, v + I.lead(s), rr); fl[key] = fl.get(key, 0.0) + float(x[col])
    keys = sorted(fl)
    duty_flow = np.array([list(kk) + [fl[kk]] for kk in keys]) if keys else np.zeros((0, 6))
    cg = cert["cg"]
    yv = np.array([[k, s, D, float(cert["mu"][(k, s)][D])] for (k, s, D) in cg.duty_key]) if cg.duty_key else np.zeros((0, 4))
    np.savez_compressed(dual_path, pi=cert["pi"], cov_first=np.array(cg.cov_first), mu=yv, duty_flow=duty_flow,
                        N=np.array(r["N"]), M=np.array(r["M"]), H=H, W=W, reset=1, obj=obj, v0=np.array(idx["v0"]), c0=idx["c0"])
    out["n_flow_arcs"] = len(keys)
    out["resets_flow"] = float(sum(fl[kk] for kk in keys if kk[4] == 1))
    out["resets_flow_late"] = float(sum(fl[kk] for kk in keys if kk[4] == 1 and kk[3] >= W))   # 应为 0
    return out


def main():
    inst, kind = sys.argv[1], sys.argv[2]
    mult = float(sys.argv[3]) if len(sys.argv) > 3 else 0
    os.makedirs(OUT, exist_ok=True); os.makedirs(DUAL, exist_ok=True)
    name = f"{inst}__{kind}" + (f"__m{mult:g}" if kind != "rate" else "")
    path = os.path.join(OUT, name + ".json")
    if os.path.exists(path): print("已完成，跳过", name); return
    I = load_train(inst); S = S_max(I)
    rec = dict(inst=inst, m=I.meta["m"], family=I.meta["family"], kind=kind, S=S)
    t0 = time.time()
    if kind == "rate":
        rec.update(rate_job(I))
    elif kind == "frRw":
        H = int(mult * S); W = 2 * S
        rec.update(mult=mult, H=H, W=W)
        rec.update(frRw_job(I, H, W, os.path.join(DUAL, name + ".npz")))
    else:
        sys.exit(f"未知种类 {kind}")
    rec["wall"] = round(time.time() - t0, 1)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: v for k, v in rec.items() if k not in ("y_g", "y_gR", "z_gR", "cycles", "path_viol_k")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
