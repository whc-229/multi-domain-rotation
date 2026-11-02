"""T15 小检查 2（11-02 审核）：只算机器人时，2·S_max 的整数最优解当启动期，之后接规则，能否用同样的资源永远可行。
对象：T15 里 MIP 解到最优的 2 个实例（m02-K1-05、m02-K12-03），不提前大修。
1. 重新求解只算机器人的完整整数模型（freshlb.fresh_lb(full_int=True, tools=False)，H = 2·S_max，HiGHS 单线程，与 T15 的 mip 相同）。
   freshlb 在 prior/ 下、不返回解向量，这里把它调用的 milp 换成记下解向量的 HiGHS 包装（与 t15_job 换 HiGHS 的做法相同），
   再按 freshlb 建模时的列顺序重建"列 → 值岗弧 (机型 k, 相位 i, 站 s, 值岗开始 D)"，核对每格每站在岗数 = K_s、机队 = 目标值。
   与仿真对齐（技术选择）：freshlb 允许在岗位开放前（D < c0 + lead_s）值岗，这些值岗不受覆盖约束、却推进机器人的相位与时间，
   仿真里不可能发生（岗位从 c0 + lead_s 才开放）。包装里把这些列的上界设为 0：仿真做不到的排班本来就不在可行集里，所以仍是下界（只会更紧）。
   另把不加这条限制的 T15 mip 结果一并记下作对照。
2. 把值岗开始 D < H 的弧写成 rost5 的计划文件（planx：首行 "end H"，每行 "s D k i 0" 一台），t < H 时按计划派，之后按规则。
   按计划派加 planxlate（11-02 加的 rost5 选项）：同型同相位的就绪候选里取就绪最晚的一台。整数解里同型同相位的机器人只差在就绪时刻，
   取最早就绪的会把后面计划要用的那台提前用掉。不加时的计划期缺岗数也记下作对照。
   规则两种：PB（T13 的手工优先级规则）与松弛规则 fresh_noreset_u（T15 只算机器人的上界用的同一个规则文件）。
3. 用 rost5 eval（notools）判定：
   - Tmax = H 不判环：计划期内的缺岗数（应为 0，说明仿真照着整数解走了）；
   - 60 万格判环：缺岗、判环结果、周期；可行 = 严格证明或 0 缺岗；
   - 不可行时：gaplog 记下前 20 个缺岗 (时刻, 站)；再逐型加 1 台（计划不变），看加哪一型能变成可行；
   - 对照：同样的台数、不用计划，只用规则。
结果写进 results/T15/plan/<实例>.json（原子写入，已存在就跳过），计划文件 <实例>.planx 放在同一目录。
用法：python3 code/t15b_plan.py <实例> [MIP 限时秒，默认 7200]"""
import os, sys, json, time, math, subprocess, tempfile, types
import numpy as np
import highspy
from scipy.sparse import vstack
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, U, load_train, inst_paths  # noqa: E402
import freshlb  # noqa: E402
from freshlb import S_max  # noqa: E402
import t13_methods as TM  # noqa: E402
import t8_rules  # noqa: E402

R5 = os.environ.get("ROST5_BIN", os.path.join(ROOT, "code", "cpp", "rost5"))
OUT = os.environ.get("T15B_OUT", os.path.join(ROOT, "results", "T15", "plan"))   # 环境变量只用于开发期试算


def highs_milp(time_limit, zero_cols, got):
    """与 t4_job._milp_highs 相同的 HiGHS 包装（单线程、mip_rel_gap = 0），另把 zero_cols 的上界设为 0，并把结果记进 got["r"]。"""
    def milp(c, constraints=None, integrality=None, bounds=None, options=None):
        """scipy.optimize.milp 兼容的调用方式（只支持 freshlb 用到的参数）。"""
        A = vstack([cc.A for cc in constraints]).tocsc()
        lo = np.concatenate([np.atleast_1d(cc.lb) * np.ones(cc.A.shape[0]) for cc in constraints])
        hi = np.concatenate([np.atleast_1d(cc.ub) * np.ones(cc.A.shape[0]) for cc in constraints])
        n, m = A.shape[1], A.shape[0]
        h = highspy.Highs(); h.setOptionValue("output_flag", False); h.setOptionValue("threads", 1)
        h.setOptionValue("time_limit", float(time_limit)); h.setOptionValue("mip_rel_gap", 0.0)
        ub = np.full(n, highspy.kHighsInf); ub[list(zero_cols)] = 0.0
        lp = highspy.HighsLp(); lp.num_col_ = n; lp.num_row_ = m
        lp.col_cost_ = np.asarray(c, float); lp.col_lower_ = np.zeros(n); lp.col_upper_ = ub
        lp.row_lower_ = lo; lp.row_upper_ = hi
        lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
        lp.a_matrix_.start_ = A.indptr; lp.a_matrix_.index_ = A.indices; lp.a_matrix_.value_ = A.data
        lp.a_matrix_.num_col_ = n; lp.a_matrix_.num_row_ = m
        lp.integrality_ = [highspy.HighsVarType.kInteger if v else highspy.HighsVarType.kContinuous for v in integrality]
        h.passModel(lp); h.run()
        info = h.getInfo(); st = h.modelStatusToString(h.getModelStatus())
        has = info.primal_solution_status >= 2
        r = types.SimpleNamespace(x=np.array(h.getSolution().col_value) if has else None, fun=info.objective_function_value if has else None,
                                  status=st, mip_dual_bound=info.mip_dual_bound, mip_gap=info.mip_gap)
        got["r"] = r; return r
    return milp


def solve_capture(I, H, tl, zero_cols):
    """求解只算机器人的完整整数模型（zero_cols 的列固定为 0），返回 (freshlb 的结果, 解向量 x)。"""
    got = {}
    milp = highs_milp(tl, zero_cols, got)
    old = freshlb.milp; freshlb.milp = milp
    try: res = freshlb.fresh_lb(I, H, tools=False, reset=False, full_int=True, time_limit=tl)
    finally: freshlb.milp = old
    return res, got["r"].x


def duty_arcs(I, H):
    """按 freshlb.fresh_lb（tools=False、reset=False）的列顺序重建值岗弧：[(列, k, i, s, D)]，以及机队变量的起始列。"""
    nk = len(I.T); v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]
    col = 0; arcs = []
    for k in range(nk):
        for i in range(I.L[k]):
            for v in range(v0[k], H):
                col += 1                                  # 等待弧
                for s in I.elig[k]:
                    ld, _ = I.step(k, i, s)
                    arcs.append((col, k, i, s, v + ld)); col += 1
    return arcs, col, v0


def ev(txt, T, cert, Nk, opts, gaplog=None):
    """rost5 eval（只算机器人），返回 JSON。"""
    a = [R5, "eval", txt, str(T), str(cert)] + [str(x) for x in Nk] + ["|", "|", "|"] + opts + ["notools"] + ([f"gaplog={gaplog}"] if gaplog else [])
    return json.loads(subprocess.run(a, capture_output=True, text=True, check=True).stdout)


def verdict(e):
    """eval 结果是否长期可行（严格证明或 0 缺岗）。"""
    return int(e["gaps"] == 0)


def main():
    name = sys.argv[1]; tl = float(sys.argv[2]) if len(sys.argv) > 2 else 7200.0
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f"{name}.json")
    if os.path.exists(path): print("已完成，跳过", name); return
    I = load_train(name); S = S_max(I); H = 2 * S; _, txt = inst_paths(name, "train")
    t0 = time.time()
    arcs, iN, v0 = duty_arcs(I, H); c0 = max(v0)
    early = [c for (c, k, i, s, D) in arcs if D < c0 + I.lead(s)]          # 岗位开放前的值岗弧，固定为 0
    res, x = solve_capture(I, H, tl, early)
    t15 = json.load(open(os.path.join(ROOT, "results", "T15", "jobs", f"{name}__mip.json")))
    rec = dict(inst=name, S=S, H=H, tl=tl, n_early_cols=len(early), mip_status=res.get("status"), mip_obj=res.get("obj"), mip_bound=res.get("bound"),
               mip_sec=res.get("sec"), t15_mip=dict(status=t15["status"], obj=t15["obj"], bound=t15["bound"]))
    if x is None: rec["error"] = "MIP 没有整数解"
    else:
        Nk = [int(round(v)) for v in x[iN:iN + len(I.T)]]
        assert len(x) == iN + len(I.T), (len(x), iN)
        # 核对：每格每站在岗数 = K_s（从 c0 + lead_s 起），机队 = 目标值，解是整数，岗位开放前没有值岗
        cover = {}
        plan = []
        for (c, k, i, s, D) in arcs:
            n = int(round(x[c]))
            if n <= 0: continue
            for q in range(I.d[k]):
                if D + q < H: cover[(s, D + q)] = cover.get((s, D + q), 0) + n
            if D < H: plan += [(s, D, k, i)] * n
        bad_cov = sum(1 for s in range(len(I.sites)) for t in range(c0 + I.lead(s), H) if cover.get((s, t), 0) != I.sites[s]["K"])
        rec.update(Nk=Nk, R=sum(Nk), n_plan=len(plan), check_cover_bad=bad_cov, check_int=float(np.max(np.abs(x - np.round(x)))),
                   check_obj=int(sum(Nk) == round(res["obj"])), check_early=int(sum(1 for (s, D, k, i) in plan if D < c0 + I.lead(s))))
        pf = os.path.join(OUT, f"{name}.planx")
        with open(pf + ".tmp", "w") as f:
            f.write(f"end {H}\n")
            for s, D, k, i in sorted(plan, key=lambda e: (e[1], e[0], e[2], e[3])): f.write(f"{s} {D} {k} {i} 0\n")
        os.replace(pf + ".tmp", pf)
        rules = {"pb": [o for o in TM.opts("pb", name, S) if o != "facloop"], "rlp": t8_rules.rule_opts(name, "rlp")}
        rec["rules"] = {}
        with tempfile.TemporaryDirectory(prefix=f"t15b_{name}_") as tmpd:
            for rn, ro in rules.items():
                o = ro + [f"planx={pf}", "planxlate"]
                r = dict(plan_gaps=ev(txt, H, 0, Nk, o)["gaps"], plan_gaps_first=ev(txt, H, 0, Nk, ro + [f"planx={pf}"])["gaps"])
                e = ev(txt, 600000, 1, Nk, o); r["eval"] = {k: e.get(k) for k in ("gaps", "cert", "period", "exact", "t_end")}; r["feasible"] = verdict(e)
                r["rule_only"] = {k: v for k, v in ev(txt, 600000, 1, Nk, ro).items() if k in ("gaps", "cert", "period", "exact")}
                if not r["feasible"]:
                    gl = os.path.join(tmpd, f"{rn}.gaps"); ev(txt, 600000, 1, Nk, o, gaplog=gl)
                    G = [tuple(map(int, ln.split()[:2])) for ln in open(gl)]
                    r["first_gaps"] = G[:20]; r["gap_sites"] = {str(s): sum(1 for _, x_ in G if x_ == s) for s in sorted({g[1] for g in G})}
                    r["add_one"] = {}
                    for k in range(len(I.T)):
                        N2 = list(Nk); N2[k] += 1
                        e2 = ev(txt, 600000, 1, N2, o)
                        r["add_one"][str(k)] = dict(gaps=e2["gaps"], cert=e2["cert"], feasible=verdict(e2), plan_gaps=ev(txt, H, 0, N2, o)["gaps"])
                rec["rules"][rn] = r
                print(name, rn, json.dumps(r, ensure_ascii=False), flush=True)
    rec["wall"] = round(time.time() - t0, 1)
    with open(path + ".tmp", "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(path + ".tmp", path)
    print(json.dumps({k: v for k, v in rec.items() if k != "rules"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
