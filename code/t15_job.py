"""T15 单个作业：下界松在哪里（m ∈ {2, 3} 的训练实例，不提前大修场景，只看机器人）。结果写进 results/T15/jobs/<名>.json（原子写入，已存在就跳过）。
"只看机器人"的实现（技术选择）：rost5 的 notools 选项在读入实例后去掉全部工具（工具检查恒成立、不计入资源，定台数只调机器人台数）；
LP / MIP 用 freshlb.fresh_lb(tools=False)。
种类：
- lp：只算机器人的全新启动 LP，H = 4·S_max 与 2·S_max 两个（不提前大修）；下界 = ⌈LP(4S)⌉；
- mip：只算机器人的完整整数模型（所有流量取整数，freshlb full_int=True），H = 2·S_max，HiGHS 单线程、限时（默认 7200 秒），记最优值或对偶界；
- ub_rlp：fresh_noreset_u 规则（T3b 的规则文件）+ notools，rost5 rlp 定台数 + 60 万格复核；
- ub_fd <种子>：T11 的 RDD-FD（默认规则 rlp，main 设置：T_eval = 6·S_max、1000 代、互换 0.5）+ notools；
- ub_es <种子>：ES curdef（certs 模式，100 代，(4 + 12)，与 T4b 相同）+ notools。
用法：python3 code/t15_job.py lp|mip|ub_rlp <实例> [限时]  |  python3 code/t15_job.py ub_fd|ub_es <实例> <种子>"""
import os, sys, json, time, subprocess, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
import freshlb  # noqa: E402
from freshlb import S_max  # noqa: E402
from t4_job import _milp_highs  # noqa: E402
import t8_rules  # noqa: E402

R5 = os.environ.get("ROST5_BIN", os.path.join(ROOT, "code", "cpp", "rost5"))
OUT = os.path.join(ROOT, "results", "T15", "jobs")


def main():
    kind, name = sys.argv[1], sys.argv[2]
    arg = sys.argv[3] if len(sys.argv) > 3 else None
    os.makedirs(OUT, exist_ok=True)
    jn = f"{name}__{kind}" + (f"__s{arg}" if kind in ("ub_fd", "ub_es") else "")
    path = os.path.join(OUT, jn + ".json")
    if os.path.exists(path): print("已完成，跳过", jn); return
    I = load_train(name); S = S_max(I); _, txt = inst_paths(name, "train")
    rec = dict(inst=name, m=I.meta["m"], family=I.meta["family"], kind=kind, S=S)
    t0 = time.time()
    if kind == "lp":
        for mult in (4, 2):
            r = freshlb.fresh_lb(I, mult * S, tools=False, reset=False, integer=False, time_limit=7200)
            rec[f"LP{mult}"] = r["obj"]; rec[f"N{mult}"] = r.get("N"); rec[f"sec{mult}"] = r.get("sec")
        rec["LB"] = int(math.ceil(rec["LP4"] - 1e-6)); rec["LB2"] = int(math.ceil(rec["LP2"] - 1e-6))
    elif kind == "mip":
        tl = float(arg) if arg else 7200.0
        old = freshlb.milp; freshlb.milp = _milp_highs(tl)
        try: r = freshlb.fresh_lb(I, 2 * S, tools=False, reset=False, full_int=True, time_limit=tl)
        finally: freshlb.milp = old
        b = r.get("bound")
        rec.update(H=2 * S, tl=tl, obj=r.get("obj"), bound=b, LB_int=int(np.ceil(b - 1e-6)) if b is not None and np.isfinite(b) else None,
                   status=r.get("status"), N=r.get("N"), n=r.get("n"), sec=r.get("sec"))
    elif kind == "ub_rlp":
        args = [R5, "rlp", txt, os.path.join(ROOT, "results", "T3", "rlp", f"{name}__fresh_noreset_u.rlp"), "365", "60000", "600000", "notools"]
        rec.update(json.loads(subprocess.run(args, capture_output=True, text=True, check=True).stdout)); rec["cmd"] = " ".join(args)
    elif kind == "ub_fd":
        args = [R5, "rddfd", txt, "365", "60000", "600000", str(6 * S), "1000", "36000", "4", "12", arg, "0.1", "0.4"] + t8_rules.rule_opts(name, "rlp") + \
               [f"W={2 * S}", "stall=30", "swap=0.5", "fdsel=best", "notools"]
        d = json.loads(subprocess.run(args, capture_output=True, text=True, check=True).stdout); d.pop("genes", None); d.pop("cdf", None); d.pop("dualw", None)
        rec.update(d); rec["cmd"] = " ".join(args)
    elif kind == "ub_es":
        args = [R5, "es", txt, "certs", "365", "100", "4", "12", arg, "60000", "600000", "curdef", "notools"]
        rec.update(json.loads(subprocess.run(args, capture_output=True, text=True, check=True).stdout)); rec["cmd"] = " ".join(args)
    else: raise ValueError(kind)
    if kind.startswith("ub"):
        rec["proven"] = int(rec.get("verify_cert") == 1 and rec.get("exact") == 1)
        rec["feasible"] = int(rec.get("R", -1) >= 0 and (rec["proven"] == 1 or rec.get("verify_gaps") == 0))
        rec["tools_in_sol"] = len(rec.get("Mj") or [])
    rec["wall"] = round(time.time() - t0, 1)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: rec.get(k) for k in ("inst", "kind", "LB", "LB2", "obj", "bound", "R", "status", "feasible", "wall")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
