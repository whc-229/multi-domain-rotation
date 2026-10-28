"""T4 单个作业（m ∈ {2, 3} 的训练实例）。结果写进 results/T4/jobs/<名>.json（原子写入，已存在就跳过）。
- lb：时间展开 ILP 下界 = freshlb.fresh_lb(full_int=True)，H = mult·S_max。建模用 freshlb 原代码，只把 scipy 的 milp 换成 highspy
  （单线程、mip_rel_gap = 0、限时）。下界取 ⌈对偶界 − 1e-6⌉；解到最优时对偶界 = 最优值。
- ub：code/fresh_tail_mdrr.py 的周期尾整数模型，参数 (T, P)；最优值（或限时内的最好可行解）是无限期全新启动问题的上界。
- esub：上界的替代来源（HiGHS / CP-SAT 在周期尾模型上找不到可行解，见报告）：更长的 ES（30 代），
  rost5 es <实例> certs 365 30 4 12 <种子> 60000 600000 <curdef | lpg+planf>；解由仿真器的严格证书与精确重放判定。
- es100（T4b）：同 esub，但 100 代、种子 0–4，结果写进 results/T4b/jobs/<实例>__es100_<配置>__s<种子>.json。
  同一种子下 ES 的随机数流与代数无关，所以前 30 代与 esub 相同（汇总时核对），逐代最好值（hist）给出 10 / 30 / 100 代曲线。
用法：python3 code/t4_job.py lb <实例> <reset 0/1> <mult> <限时秒>
      python3 code/t4_job.py esub <实例> <curdef|lpgplan> <种子>
      python3 code/t4_job.py es100 <实例> <curdef|lpgplan> <种子>
      python3 code/t4_job.py ub <实例> <reset 0/1> <T> <P> <限时秒> [mip_lp_solver]"""
import os, sys, json, time, types
import numpy as np
from scipy.sparse import vstack
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train  # noqa: E402
import highspy  # noqa: E402
import freshlb  # noqa: E402

OUT = os.environ.get("T4_OUT", os.path.join(ROOT, "results", "T4", "jobs"))   # 环境变量只用于开发期试算


def _milp_highs(time_limit):
    """返回一个与 scipy.optimize.milp 调用方式兼容的函数（只支持 freshlb 用到的参数）。"""
    def milp(c, constraints=None, integrality=None, bounds=None, options=None):
        A = vstack([cc.A for cc in constraints]).tocsc()
        lo = np.concatenate([np.atleast_1d(cc.lb) * np.ones(cc.A.shape[0]) for cc in constraints])
        hi = np.concatenate([np.atleast_1d(cc.ub) * np.ones(cc.A.shape[0]) for cc in constraints])
        n, m = A.shape[1], A.shape[0]
        h = highspy.Highs(); h.setOptionValue("output_flag", False); h.setOptionValue("threads", 1)
        h.setOptionValue("time_limit", float(time_limit)); h.setOptionValue("mip_rel_gap", 0.0)
        lp = highspy.HighsLp(); lp.num_col_ = n; lp.num_row_ = m
        lp.col_cost_ = np.asarray(c, float); lp.col_lower_ = np.zeros(n); lp.col_upper_ = np.full(n, highspy.kHighsInf)
        lp.row_lower_ = lo; lp.row_upper_ = hi
        lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
        lp.a_matrix_.start_ = A.indptr; lp.a_matrix_.index_ = A.indices; lp.a_matrix_.value_ = A.data
        lp.a_matrix_.num_col_ = n; lp.a_matrix_.num_row_ = m
        lp.integrality_ = [highspy.HighsVarType.kInteger if v else highspy.HighsVarType.kContinuous for v in integrality]
        h.passModel(lp); h.run()
        info = h.getInfo(); st = h.modelStatusToString(h.getModelStatus())
        has = info.primal_solution_status >= 2
        return types.SimpleNamespace(x=np.array(h.getSolution().col_value) if has else None,
                                     fun=info.objective_function_value if has else None, status=st,
                                     mip_dual_bound=info.mip_dual_bound, mip_gap=info.mip_gap)
    return milp


def lb_job(I, reset, mult, tl):
    H = int(mult * freshlb.S_max(I))
    old = freshlb.milp; freshlb.milp = _milp_highs(tl)
    try:
        t0 = time.time()
        r = freshlb.fresh_lb(I, H, tools=True, reset=reset, full_int=True, time_limit=tl)
        sec = time.time() - t0
    finally:
        freshlb.milp = old
    b = r.get("bound")
    return dict(H=H, obj=r.get("obj"), bound=b, LB=int(np.ceil(b - 1e-6)) if b is not None and np.isfinite(b) else None,
                status=r.get("status"), N=r.get("N"), M=r.get("M"), n=r.get("n"), sec=round(sec, 1))


OUT_B = os.path.join(ROOT, "results", "T4b", "jobs")


def esub_job(name, cfg, seed, gens=30, out=None, prefix="esub"):
    """长 ES 上界（默认 30 代；T4b 用 gens=100、prefix="es100"、out=OUT_B）。规则文件 / lpg 由 T3 的 ensure_rules 生成（加锁）。"""
    out = out or OUT
    import subprocess
    import t3_job
    from mdrr import inst_paths
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"{name}__{prefix}_{cfg}__s{seed}.json")
    if os.path.exists(path): print("已完成，跳过", path); return
    if not t3_job.deps_ok(name): print("T2 依赖未齐"); sys.exit(3)
    import fcntl
    lk = open(path + ".lock", "w")      # 两个调度器从两头跑同一张作业表时，保证同一作业只有一个进程在写
    try: fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError: print("另一进程正在跑这个作业，退出", path); sys.exit(4)
    if os.path.exists(path): print("已完成，跳过", path); return
    t3_job.ensure_rules(name)
    _, txt = inst_paths(name, "train"); R = t3_job.RDIR
    opts = ["curdef"] if cfg == "curdef" else ["lpg=" + os.path.join(R, f"{name}.lpg.txt"), "planf=" + os.path.join(R, f"{name}__fresh_rlp.planf")]
    args = [t3_job.R5, "es", txt, "certs", "365", str(gens), "4", "12", str(seed), "60000", "600000"] + opts
    t0 = time.time(); r = subprocess.run(args, capture_output=True, text=True)
    d = json.loads(r.stdout)
    I = load_train(name)
    d.update(inst=name, m=I.meta["m"], family=I.meta["family"], kind=f"{prefix}_{cfg}", gens=gens, seed=seed, wall=round(time.time() - t0, 1), cmd=" ".join(args))
    d["proven"] = int(d.get("verify_cert") == 1 and d.get("exact") == 1)
    with open(path + ".tmp", "w") as f: f.write(json.dumps(d) + "\n")
    os.replace(path + ".tmp", path)
    print(json.dumps({k: d.get(k) for k in ("inst", "kind", "seed", "R", "status", "verify_cert", "wall")}))


def main():
    if sys.argv[1] == "esub": return esub_job(sys.argv[2], sys.argv[3], int(sys.argv[4]))
    if sys.argv[1] == "es100": return esub_job(sys.argv[2], sys.argv[3], int(sys.argv[4]), gens=100, out=OUT_B, prefix="es100")
    kind, name, rs = sys.argv[1], sys.argv[2], int(sys.argv[3])
    os.makedirs(OUT, exist_ok=True)
    if kind == "lb":
        mult, tl = float(sys.argv[4]), float(sys.argv[5]); jn = f"{name}__lb__r{rs}__m{mult:g}"
    else:
        T, P, tl = int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]); jn = f"{name}__ub__r{rs}__T{T}__P{P}"
        lps = sys.argv[7] if len(sys.argv) > 7 else None
    path = os.path.join(OUT, jn + ".json")
    if os.path.exists(path): print("已完成，跳过", jn); return
    I = load_train(name)
    rec = dict(inst=name, m=I.meta["m"], family=I.meta["family"], kind=kind, reset=rs, S=freshlb.S_max(I))
    t0 = time.time()
    if kind == "lb":
        rec.update(lb_job(I, bool(rs), mult, tl)); rec["mult"] = mult
    else:
        from fresh_tail_mdrr import solve
        rec.update(solve(I, T, P, bool(rs), tl, opts={"mip_lp_solver": lps} if lps else None))
    rec["wall"] = round(time.time() - t0, 1)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: v for k, v in rec.items() if k not in ("duties", "z")}))


if __name__ == "__main__":
    main()
