"""T20 队尾重跑表（总调度第 6 段，TASKS"无人值守安排"中 (e) 的上线规则）。
(e) 的状态记在 results/T20/e_status.json：{"status": "installed" | "rejected", "md5": …, "time": …}；没有这个文件 = 还没定。
- 已装上（installed）：扫 results/T20/opp/jobs 与 results/T20/ours/jobs 里搜索型方法的结果，
  凡是"因定台数失败记无解"且不是用装了 (e) 的程序跑的（rost5 输出里没有 e_fixed 字段），排进队尾重跑。
  "因定台数失败"（11-09 审核：只重跑装了 (e) 会触发的作业，别的原因记"—"的不重跑）= R < 0、起点定不出（base_R < 0）、
  原有后备也没定出（f0_src < 0），也就是 fd / rddfd 里 (e) 的触发条件。
  例：shorteval m06-K12-13 s1 起点定出 109 个，失败在长期复核，不进重跑表。
  汇总时（11-09 审核）：重跑里 (e) 实际没触发（e_fixed = −1）的仍用旧结果，并在备查表注明。
- 后备上限（11-09 下午审核 (b)：后备最多用预算的一半）装上后（e_status.json 里 "f0cap": "installed"），另把
  "装了 (e)、没装上限时跑的（输出有 e_fixed、没有 f0_capped），后备没定出、(e) 没触发、gens = 0"的作业排进来，规则同上
  （新结果写 _e，旧结果备查；重跑中没走到 (e) 的用旧结果）。
  重跑的结果写新阶段目录 opp_e / ours_e（不覆盖旧结果）；ours_e 共用 ours 的预处理缓存（见 code/t20_job.py）。
- 没装或不装：写空表。
用法：python3 code/t20_requeue_e.py <输出作业表>     生成 q6
      python3 code/t20_requeue_e.py --final-check     退出码 0 = 可以收尾（(e) 已定，且要重跑的都有结果）"""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402
from t20_job import TLIM  # noqa: E402

SEARCH = {"opp": ("esdw", "rk", "d5", "rddhand"), "ours": ("ours", "ours_w0.25", "ours_w1", "ours_w2", "noconn", "shorteval")}
T20 = os.path.join(ROOT, "results", "T20")


def e_status():
    """(e) 的状态：installed / rejected / None（未定）。"""
    p = os.path.join(T20, "e_status.json")
    return json.load(open(p))["status"] if os.path.exists(p) else None


def f0cap_installed():
    """后备上限是否已装上。"""
    p = os.path.join(T20, "e_status.json")
    return os.path.exists(p) and json.load(open(p)).get("f0cap") == "installed"


def need_redo(r, cap):
    """一个结果是否要重跑：(i) 装 (e) 前跑的、因定台数失败记无解；(ii) 上限已装时，装上限前跑的、后备用完预算 (e) 没触发。"""
    if r.get("R", -1) >= 0 or r.get("base_R", 0) >= 0 or r.get("f0_src", 0) >= 0: return False
    if "e_fixed" not in r: return True
    return cap and "f0_capped" not in r and r["e_fixed"] < 0 and r.get("gens", -1) == 0


def to_redo():
    """要重跑的搜索型作业（见 need_redo）：[(阶段, 方法, 实例, 种子)]。"""
    out = []; cap = f0cap_installed()
    for stage, meths in SEARCH.items():
        for f in sorted(glob.glob(os.path.join(T20, stage, "jobs", "*.json"))):
            meth, inst, s = os.path.basename(f)[:-5].split("__")
            if meth not in meths: continue
            try: r = json.load(open(f))
            except Exception: continue
            if need_redo(r, cap): out.append((stage, meth, inst, int(s[1:])))
    return out


def main():
    """生成 q6，或做收尾检查。"""
    st = e_status()
    if sys.argv[1] == "--final-check":
        if st is None: sys.exit(1)
        if st != "installed": sys.exit(0)
        miss = [x for x in to_redo() if not os.path.exists(os.path.join(T20, x[0] + "_e", "jobs", f"{x[1]}__{x[2]}__s{x[3]}.json"))]
        sys.exit(1 if miss else 0)
    lines = []
    if st == "installed":
        for stage, meth, inst, s in to_redo():
            TL = TLIM[int(inst[1:3])]
            lines.append(f"{2 * TL + 1800} job {stage}_e {meth} {inst} {s}")
    with open(sys.argv[1], "w") as f: f.write("\n".join(lines) + ("\n" if lines else ""))
    print("q6", st, len(lines), "行")


if __name__ == "__main__":
    main()
