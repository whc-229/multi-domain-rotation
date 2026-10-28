"""T3c 诊断：LP 提前大修计划为什么没省下资源。只记录观察，不改规则。
对每个实例、每个变体（fresh_rlp_u 为主，fresh_noreset_u 作对照）输出四样东西：
1. 计划的提前大修（型、相位、时刻、次数）与实际执行的逐条对照（先到期先执行，与仿真器消耗计划的顺序相同）；
   没执行的给出原因线索：到期后该型该相位有没有被派过工。
2. 各型"维修中"台数随时间的曲线，与 LP 对应量（Σ 流量）对比；同时给出值岗中、空闲的台数。
   维修中 = 值岗结束到下一次就绪基准 rv 之间（含赶路）；LP 侧对每条值岗流用同一个 step 算同一区间。
3. 定台数循环逐步记录：Python 按 C++ size_res 的同样步骤重放（补足 → 回删 → 证书补台，起点倍数 100/115/130/150），
   每次补台记录触发它的缺岗（站、时刻、归因到哪型或哪种工具）和各候选资源加一台后的缺岗数。最后核对台数与 C++ 结果相同。
4. 在 LP 取整台数（定台数起点）下第一处缺岗前后 20 格内的派工与缺岗记录。
仿真用 code/rost5_ref.py（与 rost5 逐组一致，trace 只记录、不改决策）。
输出：results/T3c/<实例>__<变体>.json（全部记录）、results/T3c/<实例>__<变体>__occ.csv（曲线）、报告/T3c_figs/<实例>.png。
用法：python3 code/t3c_diag.py <实例名> [变体,...] [plot（只重画图）]"""
import os, sys, json, csv
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train  # noqa: E402
import rost5_ref as ref  # noqa: E402
import t3_make_rlp  # noqa: E402

RDIR = t3_make_rlp.OUT
OUT = os.path.join(ROOT, "results", "T3c")
FIG = os.path.join(ROOT, "报告", "T3c_figs")
HS, TCERT, MAX_ADD = 365, 60000, 60          # 与 t3_job 的 rost5 rlp 参数相同


def load_rule(I, name, tag):
    """读规则文件和（fresh_rlp* 时的）planf。"""
    R = ref.parse_rlp(I, os.path.join(RDIR, f"{name}__{tag}.rlp"))
    pf = os.path.join(RDIR, f"{name}__{tag}.planf")
    return R, (ref.parse_planf(I, pf) if tag.startswith("fresh_rlp") else None)


def why_vec(I, tr):
    """把 trace 里的缺岗归因汇总成 C++ 的 sig 向量（前 nk 个是机器人型，后 nt 个是工具）。"""
    nk, nt = len(I.T), len(tr_tools(I))
    sig = [0] * (nk + nt)
    for (_, _, _, why) in tr["gaps"]:
        for x in why[1]: sig[x if why[0] == "r" else nk + x] += 1
    return sig


def tr_tools(I):
    """工具种类列表（只为取种类数）。"""
    return ref.ToolInst(I).TT


def res_names(I):
    """资源下标 → 名字（机器人型名、工具 j）。"""
    return list(I.tn) + [f"工具{j}" for j in range(len(tr_tools(I)))]


def gap_summary(I, tr, n=5):
    """前 n 个缺岗：(t, 站, 候选数, 归因名字)。"""
    nm = res_names(I); nk = len(I.T)
    return [dict(t=t, site=s, cand=c, why=[nm[x] if why[0] == "r" else nm[nk + x] for x in why[1]]) for (t, s, c, why) in tr["gaps"][:n]]


def size_py(I, R, planf, log):
    """按 C++ size_res 的步骤重放定台数循环，log 收集每一步。返回 (Nk, Mj, 状态, 起点倍数)。"""
    nk = len(I.T); nt = len(tr_tools(I)); nm = res_names(I)

    def sim(Nk, Mj, T, cert, light=False):
        tr = {"light": light}
        r = ref.run(I, Nk, Mj, T, cert, rlp=R, planf=planf, trace=tr)
        return r, tr

    for fac in (100, 115, 130, 150):
        Nk = [(x * fac + 99) // 100 for x in R["N0"]]; Mj = [(x * fac + 99) // 100 for x in R["M0"]]
        log.append(dict(phase="起点", factor=fac, Nk=list(Nk), Mj=list(Mj)))
        V = lambda: Nk + Mj

        def bump(idx, d):
            if idx < nk: Nk[idx] += d
            else: Mj[idx - nk] += d
        r, tr = sim(Nk, Mj, HS, False); added = 0
        while r["gaps"] > 0 and added < MAX_ADD:
            sig = why_vec(I, tr); bg = bi = bs = None; cand = []
            for idx in range(nk + nt):
                bump(idx, 1); g = sim(Nk, Mj, HS, False)[0]["gaps"]; bump(idx, -1); cand.append(g)
                if bi is None or g < bg or (g == bg and sig[idx] > bs): bg, bi, bs = g, idx, sig[idx]
            log.append(dict(phase="补足", gaps=r["gaps"], first_gaps=gap_summary(I, tr), sig=dict(zip(nm, sig)),
                            cand_gaps=dict(zip(nm, cand)), chosen=nm[bi], gaps_after=bg))
            bump(bi, 1); added += 1; r, tr = sim(Nk, Mj, HS, False)
        if r["gaps"] > 0:
            log.append(dict(phase="补足失败", gaps=r["gaps"])); continue
        improved = True
        while improved:
            improved = False
            for idx in range(nk + nt):
                v = V()[idx]
                if v == 0: continue
                bump(idx, -1)
                if sim(Nk, Mj, HS, False)[0]["gaps"] == 0:
                    improved = True; log.append(dict(phase="回删", removed=nm[idx])); break
                bump(idx, 1)
        for _ in range(MAX_ADD):
            c, tr = sim(Nk, Mj, TCERT, True, light=True)
            if c["cert"] == -1 and c["gaps"] == 0: c, tr = sim(Nk, Mj, 4 * TCERT, True, light=True)
            if c["cert"] == 1: log.append(dict(phase="证书通过", period=c["period"])); return Nk, Mj, "cert", fac
            if c["cert"] == -1 and c["gaps"] == 0: log.append(dict(phase="判不出循环、无缺岗")); return Nk, Mj, "uncert", fac
            sig = why_vec(I, tr); bi = int(np.argmax(sig))
            log.append(dict(phase="证书补台", cert=c["cert"], gaps=c["gaps"], first_gaps=gap_summary(I, tr), sig=dict(zip(nm, sig)), chosen=nm[bi]))
            bump(bi, 1)
        log.append(dict(phase="证书补台失败"))
    return None, None, "fail", None


def plan_vs_exec(I, tr, H):
    """计划提前大修与实际执行的逐条对照。计划按 (型, 相位) 先到期先执行，与仿真器一致。"""
    out = []
    ex = {}
    for (t, s, q, r, k, i, planned, rv, tend) in tr["disp"]:
        if planned: ex.setdefault((k, i), []).append(t)
    for (k, i), lst in sorted(tr["plan"].items()):
        units = [t0 for (t0, c) in lst for _ in range(c)]
        done = ex.get((k, i), [])
        for u, t0 in enumerate(units):
            if u < len(done):
                out.append(dict(type=I.tn[k], phase=i, t_plan=t0, t_exec=done[u], delay=done[u] - t0))
            else:
                # 没执行：到期后该型该相位有没有被派过工（派过就会执行，所以这里一定是没派过）；该型其他相位派了多少次
                same = sum(1 for d in tr["disp"] if d[4] == k and d[5] == i and d[0] >= t0)
                other = sum(1 for d in tr["disp"] if d[4] == k and d[5] != i and d[0] >= t0)
                out.append(dict(type=I.tn[k], phase=i, t_plan=t0, t_exec=None, delay=None,
                                why=f"到期后该型相位 {i} 被派工 {same} 次，其他相位 {other} 次"))
    return out


def occupancy(I, tr, Nk, H):
    """仿真侧每格每型：值岗中、维修中（值岗结束到 rv）、空闲（其余；含初始准备）。"""
    nk = len(I.T)
    duty = np.zeros((nk, H)); mnt = np.zeros((nk, H))
    for (t, s, q, r, k, i, planned, rv, tend) in tr["disp"]:
        duty[k, max(t, 0):min(tend, H)] += 1
        mnt[k, max(tend, 0):min(rv, H)] += 1
    idle = np.array(Nk, float)[:, None] - duty - mnt
    return duty, mnt, idle


def occupancy_lp(I, name, kind, H):
    """LP 侧：对每条值岗流 (k, i, s, D, r, x)，值岗 [D, D + d_k)、维修 [D + d_k, D − ld + nx)，nx 取相位（提前大修时取 L−1）的 step。"""
    d = np.load(t3_make_rlp.dual_path(name, kind))
    nk = len(I.T); duty = np.zeros((nk, H)); mnt = np.zeros((nk, H))
    for (k, i, s, D, r, x) in d["duty_flow"]:
        k, i, s, D, r = int(k), int(i), int(s), int(D), int(r)
        ld, nx = I.step(k, I.L[k] - 1 if r else i, s)
        duty[k, max(D, 0):min(D + I.d[k], H)] += x
        mnt[k, max(D + I.d[k], 0):min(D - ld + nx, H)] += x
    return duty, mnt, d["N"]


def first_gap_window(I, R, planf, w=20, T=3000):
    """LP 取整台数下第一处缺岗前后 w 格的派工与缺岗记录。"""
    Nk, Mj = list(R["N0"]), list(R["M0"]); tr = {}; gl = []
    ref.run(I, Nk, Mj, T, False, rlp=R, planf=planf, trace=tr, gap_log=gl)
    if not tr["gaps"]: return dict(Nk=Nk, Mj=Mj, first_gap=None)
    tg = tr["gaps"][0][0]; nm = res_names(I); nk = len(I.T)
    disp = [dict(t=t, site=s, robot=r, type=I.tn[k], phase=i, reset=bool(p), next_rv=rv) for (t, s, q, r, k, i, p, rv, te) in tr["disp"] if tg - w <= t <= tg + w]
    gaps = [dict(t=t, site=s, cand=c, why=[nm[x] if why[0] == "r" else nm[nk + x] for x in why[1]]) for (t, s, c, why) in tr["gaps"] if tg - w <= t <= tg + w]
    lags = [g for g in gl if tg - w <= g["t"] <= tg + w]
    return dict(Nk=Nk, Mj=Mj, first_gap=dict(t=tg, site=tr["gaps"][0][1]), n_gaps_T=len(tr["gaps"]), dispatch=disp, gaps=gaps, gap_lags=lags)


def diag(name, tag):
    """一个 (实例, 变体) 的全部诊断。"""
    I = load_train(name); R, planf = load_rule(I, name, tag); H = R["H"]
    kind = "frR" if tag.startswith("fresh_rlp") else "fr"
    cjob = json.load(open(os.path.join(ROOT, "results", "T3b" if tag.endswith("_u") else "T3", "jobs", f"{name}__{tag}.json")))
    log = []; Nk, Mj, st, fac = size_py(I, R, planf, log)
    same = (Nk == cjob["Nk"] and Mj == cjob["Mj"])
    res = dict(inst=name, tag=tag, H=H, sizing_log=log, py_Nk=Nk, py_Mj=Mj, py_status=st, py_factor=fac,
               c_Nk=cjob["Nk"], c_Mj=cjob["Mj"], c_status=cjob["status"], c_factor=cjob["start_factor"], same_as_cpp=same)
    views = {}
    for lab, (nk_, mj_) in (("起点", (list(R["N0"]), list(R["M0"]))), ("最终", (cjob["Nk"], cjob["Mj"]))):
        tr = {}; r = ref.run(I, nk_, mj_, H, False, rlp=R, planf=planf, trace=tr)
        views[lab] = dict(Nk=nk_, Mj=mj_, gaps_H=r["gaps"], resets_H=r["resets"],
                          plan_vs_exec=plan_vs_exec(I, tr, H) if planf is not None else [], occ=occupancy(I, tr, nk_, H))
    lpd, lpm, lpN = occupancy_lp(I, name, kind, H)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, f"{name}__{tag}__occ.csv"), "w", newline="") as f:
        wr = csv.writer(f); wr.writerow(["t", "type", "lp_duty", "lp_maint", "lp_idle", "sim_start_duty", "sim_start_maint", "sim_start_idle", "sim_final_duty", "sim_final_maint", "sim_final_idle"])
        for k in range(len(I.T)):
            for t in range(H):
                wr.writerow([t, I.tn[k], round(lpd[k, t], 4), round(lpm[k, t], 4), round(lpN[k] - lpd[k, t] - lpm[k, t], 4)]
                            + [int(views[v]["occ"][x][k, t]) for v in ("起点", "最终") for x in range(3)])
    # 曲线的分段平均（每段 H/8 格）
    seg = []
    for k in range(len(I.T)):
        for a in range(0, H, H // 8):
            b = min(H, a + H // 8)
            seg.append(dict(type=I.tn[k], t0=a, t1=b, lp_maint=float(lpm[k, a:b].mean()), lp_idle=float((lpN[k] - lpd[k, a:b] - lpm[k, a:b]).mean()),
                            start_maint=float(views["起点"]["occ"][1][k, a:b].mean()), final_maint=float(views["最终"]["occ"][1][k, a:b].mean()),
                            start_idle=float(views["起点"]["occ"][2][k, a:b].mean()), final_idle=float(views["最终"]["occ"][2][k, a:b].mean())))
    res["occ_segments"] = seg; res["lp_N"] = [float(x) for x in lpN]
    for v in views.values(): v.pop("occ")
    res["views"] = views
    res["first_gap_window"] = first_gap_window(I, R, planf)
    with open(os.path.join(OUT, f"{name}__{tag}.json.tmp"), "w") as f: json.dump(res, f, ensure_ascii=False, indent=1)
    os.replace(os.path.join(OUT, f"{name}__{tag}.json.tmp"), os.path.join(OUT, f"{name}__{tag}.json"))
    return res, (lpm, views, I, H)


def plot(name, tags):
    """每型一行：LP 维修中台数 vs 仿真（起点、最终台数），各变体一列。只读 results/T3c 的 csv，可单独重画。"""
    import glob
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for fp in glob.glob("/usr/share/fonts/opentype/noto/NotoSerifCJK-Regular.ttc"): font_manager.fontManager.addfont(fp)
    plt.rcParams["font.family"] = ["Noto Serif CJK JP", "DejaVu Sans"]
    os.makedirs(FIG, exist_ok=True)
    rows = {}
    for tag in tags:
        with open(os.path.join(OUT, f"{name}__{tag}__occ.csv")) as f:
            for r in csv.DictReader(f): rows.setdefault((tag, r["type"]), []).append(r)
    types = list(dict.fromkeys(t for (_, t) in rows))
    fig, ax = plt.subplots(len(types), len(tags), figsize=(6 * len(tags), 1.8 * len(types)), squeeze=False, sharex=True)
    for c, tag in enumerate(tags):
        for k, ty in enumerate(types):
            rr = rows[(tag, ty)]; a = ax[k][c]
            a.plot([float(r["lp_maint"]) for r in rr], color="k", lw=1.2, label="LP")
            a.plot([float(r["sim_start_maint"]) for r in rr], lw=0.7, alpha=0.8, label="仿真（起点台数）")
            a.plot([float(r["sim_final_maint"]) for r in rr], lw=0.7, alpha=0.8, label="仿真（最终台数）")
            a.set_ylabel(ty, fontsize=8)
        ax[0][c].set_title(f"{name} {tag}：维修中台数", fontsize=9)
    ax[0][0].legend(fontsize=7, loc="upper right")
    fig.supxlabel("时刻（格）", fontsize=8); fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"{name}.png"), dpi=110); plt.close(fig)


if __name__ == "__main__":
    name = sys.argv[1]
    tags = sys.argv[2].split(",") if len(sys.argv) > 2 else ["fresh_rlp_u", "fresh_noreset_u"]
    if len(sys.argv) > 3 and sys.argv[3] == "plot": plot(name, tags); sys.exit()
    for tag in tags:
        res, _ = diag(name, tag)
        print(name, tag, "Python 重放与 C++ 台数相同" if res["same_as_cpp"] else f"不同：py {res['py_Nk']} {res['py_Mj']} c {res['c_Nk']} {res['c_Mj']}", flush=True)
    plot(name, tags)
