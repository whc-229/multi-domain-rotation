"""T13 RH-MIP：滚动时域整数规划（对手，模型预测式闭环）。只用实例本身的信息（不用我们的 LP 规则或计划）。
每一步（起点 ff = j·Δ，时域 [ff, H)，H = ff + L_w）：
  1. 状态：用 Python 参考仿真（与 rost5 逐位一致）按已提交的计划（planx，终点 ff）+ PB 从 0 推进到 ff，取此刻的真实状态
     （各机器人的型、相位、就绪虚拟时刻，各工具的空闲时刻，各岗位下一次到点时刻）；
  2. 窗口模型：freshlb 的时间展开网络（机器人按 (型, 相位, 虚拟时刻)、工具按 (类, 时刻) 聚合），只建 [ff, H) 内的节点与弧，全部取整数：
     - 供给 = 状态里的机器人与工具；第 0 步没有状态，供给是机队变量 N_k、M_j（代价 1，即定机队）；之后只能加台（代价 1，在 ff 加入）；
     - 覆盖：每格每站 K_s 台在岗，已在岗的班次扣掉；缺岗松弛变量，代价 GAP = 100（保证模型总有可行解）；
     - 提前大修弧只在值岗开始 D < W = 2·S_max 时有；工具卸载按 q_k·最短卸载（与 freshlb 相同，比仿真早还，是松弛）；
  3. 热启动：同一个仿真按"已提交计划 + PB"推进到 H 的排班换算成模型变量（缺岗记在松弛上），交给 HiGHS 作为初始解；
  4. HiGHS（单线程，mip_rel_gap = 0）限时求解，取最好的整数解；把 D ∈ [ff, ff + Δ) 的值岗弧提交进计划；j + 1。
  滚到 4·S_max 为止。计划之后按 PB。机队 = 第 0 步的 N、M 加之后的加台。
最后用 rost5 size 定台数（起点 = 上面的机队，规则 = 按计划派 + PB 兜底，放大起点循环）+ 60 万格复核，与其他方法相同的判定。
技术选择（写进报告）：L_w = S_max，Δ = ⌈S_max/2⌉，每步限时 TL 秒，总限时 TOT 秒（到了就停止滚动，之后按 PB）。
结果写进 <输出目录>/<实例>.json（原子写入，已存在就跳过）。
用法：python3 code/t13_rhmip.py <实例> [每步限时秒，默认 60] [总限时秒，默认 3600] [输出目录]"""
import os, sys, json, time, math, subprocess
from collections import defaultdict
import numpy as np
import highspy
from scipy.sparse import coo_matrix
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, U, ToolInst, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import rost5_ref as ref  # noqa: E402
import t13_methods as TM  # noqa: E402

R5 = os.environ.get("ROST5_BIN", os.path.join(ROOT, "code", "cpp", "rost5"))   # 环境变量只用于安装前试算
GAP = 100.0


class Model:
    """窗口 [ff, H) 的整数模型。变量键：("w",k,i,v) 等待、("d",k,i,v,s,r) 值岗、("tw",j,v) 工具等待、("td",k,s,D,j) 工具占用、
    ("N",k)/("M",j) 第 0 步的机队、("E",k)/("F",j) 之后的加台、("u",s,t) 缺岗松弛。"""

    def __init__(self, I, TI, ff, H, W, state):
        self.I, self.TI, self.ff, self.H = I, TI, ff, H
        nk, nt = len(I.T), TI.nt
        v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]; c0 = max(v0)
        maxld = max(I.step(k, i, s)[0] for k in range(nk) for i in range(I.L[k]) for s in I.elig[k])
        self.vlo = vlo = 0 if state is None else max(0, ff - maxld)     # 机器人节点从 ff − 最大前置开始（值岗开始 D 仍须 ≥ ff）
        self.keys, cost = [], []; self.col = {}
        rows = {}; R, C, V = [], [], []
        rhs = defaultdict(float)

        def row(key):
            if key not in rows: rows[key] = len(rows)
            return rows[key]

        def var(key, c=0.0):
            self.col[key] = len(self.keys); self.keys.append(key); cost.append(c); return self.col[key]

        def add(r, a, v): R.append(r); C.append(a); V.append(v)
        # 覆盖需求：K_s − 已在岗（岗位下一次到点时刻 > t 的）
        busy = defaultdict(int)
        if state is not None:
            for q, f in enumerate(state["fr"]):
                for t in range(ff, min(f, H)): busy[(state["psite"][q], t)] += 1
        for s in range(len(I.sites)):
            for t in range(max(ff, c0 + I.lead(s)), H):
                rhs[row(("c", s, t))] = I.sites[s]["K"] - busy[(s, t)]
                add(row(("c", s, t)), var(("u", s, t), GAP), 1)
        duty = defaultdict(list)
        for k in range(nk):
            L = I.L[k]
            for i in range(L):
                for v in range(max(v0[k], vlo), H):
                    a = var(("w", k, i, v)); add(row(("n", k, i, v)), a, -1)
                    if v + 1 < H: add(row(("n", k, i, v + 1)), a, 1)
                    for s in I.elig[k]:
                        ld, nx = I.step(k, i, s); D = v + ld
                        if D < ff: continue
                        opts = [(nx, (i + 1) % L, 0)] + ([(I.step(k, L - 1, s)[1], 0, 1)] if i != L - 1 and D < W else [])
                        for (dv, ni, r) in opts:
                            a = var(("d", k, i, v, s, r)); add(row(("n", k, i, v)), a, -1)
                            if v + dv < H: add(row(("n", k, ni, v + dv)), a, 1)
                            for q in range(I.d[k]):
                                if ("c", s, D + q) in rows: add(rows[("c", s, D + q)], a, 1)
                            if D < H and nt: duty[(k, s, D)].append(a)
        tlo = 0 if state is None else ff
        for j in range(nt):
            for v in range(tlo, H):
                a = var(("tw", j, v)); add(row(("t", j, v)), a, -1)
                if v + 1 < H: add(row(("t", j, v + 1)), a, 1)
        for (k, s, D), cols in duty.items():
            Rk = TI.R[k]; pre = TI.pre(k, s); unin_lb = Rk["q"] * min(TI.TT[jj]["uninstall"] for jj in Rk["C"])
            qr = row(("q", k, s, D))
            for a in cols: add(qr, a, -Rk["q"])
            for j in Rk["C"]:
                a0, b0 = max(tlo, D - pre), D + I.d[k] + TI.post(k, s, j, unin_lb)
                a = var(("td", k, s, D, j)); add(row(("t", j, a0)), a, -1)
                if b0 < H: add(row(("t", j, b0)), a, 1)
                add(qr, a, 1)
        # 供给
        if state is None:
            for k in range(nk): add(row(("n", k, 0, v0[k])), var(("N", k), 1.0), 1)
            for j in range(nt): add(row(("t", j, 0)), var(("M", j), 1.0), 1)
        else:
            for r, k in enumerate(state["rk"]):
                v = max(state["rv"][r], vlo, v0[k])
                if v < H: rhs[row(("n", k, state["rph"][r], v))] -= 1
            for j in range(nt):
                for f in state["tf"][j]:
                    v = max(f, ff)
                    if v < H: rhs[row(("t", j, v))] -= 1
            for k in range(nk):
                v = max(ff, v0[k])
                if v < H: add(row(("n", k, 0, v)), var(("E", k), 1.0), 1)
            for j in range(nt): add(row(("t", j, ff)), var(("F", j), 1.0), 1)
        self.rows = rows; self.cost = np.array(cost)
        self.b = np.array([rhs[r] for r in range(len(rows))])
        self.A = coo_matrix((V, (R, C)), shape=(len(rows), len(self.keys))).tocsc()

    def warm(self, disp, state, gaps_cover, fleet):
        """把仿真排班（[ff, H) 内的派工）换算成变量值。disp：派工记录（rost5_ref 的 trace）；fleet：第 0 步的 (N, M)。"""
        I, TI, ff, H = self.I, self.TI, self.ff, self.H
        x = np.zeros(len(self.keys)); nk = len(I.T)
        v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]
        def inc(key, c=1):
            if key in self.col: x[self.col[key]] += c
            else: raise KeyError(key)
        # 机器人：逐台走路径
        if state is None:
            rk = [k for k in range(nk) for _ in range(fleet[0][k])]; pos = {r: (rk[r], 0, v0[rk[r]]) for r in range(len(rk))}
            for k in range(nk): inc(("N", k), fleet[0][k])
        else:
            pos = {r: (k, state["rph"][r], max(state["rv"][r], self.vlo, v0[k])) for r, k in enumerate(state["rk"])}
        byr = defaultdict(list)
        for e in disp:
            if ff <= e[0] < H: byr[e[3]].append(e)
        for r, (k, i, u) in pos.items():
            for e in sorted(byr.get(r, []), key=lambda e: e[0]):
                t, s, _, _, k2, i2, _, rvn, _, v, _, rs = e
                for vv in range(u, v): inc(("w", k, i, vv))
                inc(("d", k, i, v, s, rs)); i, u = (0 if rs else (i + 1) % I.L[k]), rvn
                if u >= H: break
            else:
                for vv in range(u, H): inc(("w", k, i, vv))
        # 工具：逐件走路径
        tpos = {}
        if state is None:
            for j in range(TI.nt):
                inc(("M", j), fleet[1][j])
                for idx in range(fleet[1][j]): tpos[(j, idx)] = 0
        else:
            for j in range(TI.nt):
                for idx, f in enumerate(state["tf"][j]): tpos[(j, idx)] = max(f, ff)
        tl = 0 if state is None else ff
        uses = defaultdict(list)
        for e in disp:
            if ff <= e[0] < H:
                for (j, idx) in e[10]: uses[(j, idx)].append(e)
        Rk = TI.R
        for (j, idx), u in tpos.items():
            for e in sorted(uses.get((j, idx), []), key=lambda e: e[0]):
                t, s, k = e[0], e[1], e[4]; pre = TI.pre(k, s)
                a0 = max(tl, t - pre); b0 = t + I.d[k] + TI.post(k, s, j, Rk[k]["q"] * min(TI.TT[jj]["uninstall"] for jj in Rk[k]["C"]))
                for vv in range(u, a0): inc(("tw", j, vv))
                inc(("td", k, s, t, j)); u = b0
                if u >= H: break
            else:
                for vv in range(u, H): inc(("tw", j, vv))
        for (s, t), g in gaps_cover.items():
            if ("u", s, t) in self.col: x[self.col[("u", s, t)]] = g
        return x

    def solve(self, tl, x0=None):
        """HiGHS 限时求整数解（可带初始解）。返回 (状态, 目标, x)。"""
        n = len(self.keys)
        h = highspy.Highs(); h.setOptionValue("output_flag", False); h.setOptionValue("threads", 1)
        h.setOptionValue("mip_rel_gap", 0.0); h.setOptionValue("time_limit", float(tl))
        lp = highspy.HighsLp(); lp.num_col_ = n; lp.num_row_ = len(self.b)
        lp.col_cost_ = self.cost; lp.col_lower_ = np.zeros(n); lp.col_upper_ = np.full(n, highspy.kHighsInf)
        lp.row_lower_ = self.b; lp.row_upper_ = self.b
        lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
        lp.a_matrix_.start_ = self.A.indptr; lp.a_matrix_.index_ = self.A.indices; lp.a_matrix_.value_ = self.A.data
        lp.a_matrix_.num_col_ = n; lp.a_matrix_.num_row_ = len(self.b)
        lp.integrality_ = [highspy.HighsVarType.kInteger] * n
        h.passModel(lp)
        if x0 is not None:
            sol = highspy.HighsSolution(); sol.col_value = list(map(float, x0)); sol.value_valid = True; h.setSolution(sol)
        h.run()
        st = h.modelStatusToString(h.getModelStatus()); s = h.getSolution()
        if not s.value_valid: return st, None, None
        return st, float(h.getInfo().objective_function_value), np.round(np.array(s.col_value)).astype(int)


def simulate(I, Nk, Mj, plan, end, W, H, snap_at):
    """Python 参考仿真：按计划（终点 end）+ PB 推进到 H。返回 (trace, 结果)。"""
    pol = dict(w=[float(x) for x in TM.W_PB] + [0.0] * 6, plan_until=W, planx=dict(end=end, e=plan))
    tr = {"snap_at": snap_at}
    r = ref.run(I, Nk, Mj, H, False, trace=tr, **pol)
    return tr, r


def rolling(I, S, txt, tl, tot):
    """滚动求解，返回 (计划 dict, 终点, 机队 N、M, 逐步记录)。"""
    TI = ToolInst(I); W = 2 * S; Tplan = 4 * S; Lw = S; step = int(math.ceil(S / 2)); t0 = time.time()
    pb = json.loads(subprocess.run([R5, "size", txt, "365", "60000", "60000"] + TM.opts("pb", I.name if hasattr(I, "name") else "", S),
                                   capture_output=True, text=True, check=True).stdout)
    if pb["R"] < 0: return {}, 0, None, None, [dict(step=-1, status="PB 定不出台数，无热启动")]
    N, M = list(pb["Nk"]), list(pb["Mj"]); plan = defaultdict(list); log = []; j = 0
    while j * step < Tplan:
        ff = j * step; H = ff + Lw
        if time.time() - t0 > tot: log.append(dict(step=j, status="总限时到，之后按 PB")); break
        tr, rs = simulate(I, N, M, dict(plan), ff, W, H, ff if j > 0 else None)
        st0 = None
        if j > 0:
            sn = tr["snap"]; st0 = dict(rk=[k for k in range(len(I.T)) for _ in range(N[k])], rph=sn["rph"], rv=sn["rv"], tf=sn["tf"], fr=sn["fr"], psite=sn["psite"])
        Mo = Model(I, TI, ff, H, W, st0)
        # 热启动的缺岗：覆盖行的需求（已扣在岗）减去仿真在 [ff, H) 内派出的班次的覆盖
        covw = defaultdict(int)
        for e in tr["disp"]:
            if ff <= e[0] < H:
                for q in range(I.d[e[4]]): covw[(e[1], e[0] + q)] += 1
        gaps = {(key[1], key[2]): int(round(Mo.b[r])) - covw[(key[1], key[2])] for key, r in Mo.rows.items() if key[0] == "c"}
        x0 = Mo.warm(tr["disp"], st0, gaps, (N, M))
        res0 = Mo.A @ x0 - Mo.b
        t1 = time.time(); st, obj, x = Mo.solve(min(tl, max(1.0, tot - (time.time() - t0))), x0)
        warm_obj = float(Mo.cost @ x0)
        log.append(dict(step=j, from_=ff, H=H, n=len(Mo.keys), warm_obj=warm_obj, warm_resid=float(np.abs(res0).max()) if len(res0) else 0.0,
                        status=st, obj=obj, sec=round(time.time() - t1, 1)))
        if x is None: break
        if j == 0:
            N = [int(x[Mo.col[("N", k)]]) for k in range(len(I.T))]; M = [int(x[Mo.col[("M", q)]]) for q in range(TI.nt)]
        else:
            N = [N[k] + int(x[Mo.col[("E", k)]]) for k in range(len(I.T))]; M = [M[q] + int(x[Mo.col[("F", q)]]) for q in range(TI.nt)]
        for key, a in Mo.col.items():
            if key[0] == "d" and x[a] > 0:
                _, k, i, v, s, r = key; D = v + I.step(k, i, s)[0]
                if ff <= D < ff + step: plan[(s, D)] += [(k, i, r)] * int(x[a])
        j += 1
    end = min(j * step, Tplan)
    return {k: v for k, v in plan.items() if k[1] < end}, end, N, M, log


def main():
    name = sys.argv[1]; tl = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0; tot = float(sys.argv[3]) if len(sys.argv) > 3 else 3600.0
    out = sys.argv[4] if len(sys.argv) > 4 else os.path.join(ROOT, "results", "T13", "rhmip")
    os.makedirs(out, exist_ok=True); path = os.path.join(out, f"{name}.json")
    if os.path.exists(path): print("已完成，跳过", name); return
    I = load_train(name); I.name = name; S = S_max(I); W = 2 * S; _, txt = inst_paths(name, "train")
    t0 = time.time()
    plan, end, N, Mv, log = rolling(I, S, txt, tl, tot)
    sec_mip = time.time() - t0
    pf = os.path.join(out, f"{name}.planx")
    with open(pf + ".tmp", "w") as f:
        f.write(f"end {end}\n")
        for (s, D) in sorted(plan, key=lambda x: (x[1], x[0])):
            for (k, i, r) in plan[(s, D)]: f.write(f"{s} {D} {k} {i} {r}\n")
    os.replace(pf + ".tmp", pf)
    opts = TM.opts("pb", name, S) + [f"planx={pf}", f"planW={W}"]
    if N is not None: opts += ["N0=" + ",".join(map(str, N)), "M0=" + ",".join(map(str, Mv))]
    r = json.loads(subprocess.run([R5, "size", txt, "365", "60000", "600000"] + opts, capture_output=True, text=True, check=True).stdout)
    rec = dict(inst=name, m=I.meta["m"], family=I.meta["family"], method="rhmip", S=S, W=W, tl_step=tl, tl_total=tot,
               mip_N=N, mip_M=Mv, mip_R=(sum(N) + sum(Mv)) if N else None, plan_end=end, plan_duties=sum(len(v) for v in plan.values()), steps=log,
               sec_mip=round(sec_mip, 1), wall=round(time.time() - t0, 1), cmd_opts=opts, **r)
    rec["proven"] = int(rec.get("verify_cert") == 1 and rec.get("exact") == 1)
    rec["window_ok"] = int(rec.get("last_reset", -1) < W)
    with open(path + ".tmp", "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(path + ".tmp", path)
    print(json.dumps({k: rec[k] for k in ("inst", "mip_R", "plan_end", "R", "status", "verify_gaps", "sec_mip", "wall")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
