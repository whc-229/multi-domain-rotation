"""T20 作业表（code/t20_run.sh 的输入）。每行 "<超时秒> job <阶段> <方法> <实例> <种子>" 或 "<超时秒> lb <实例> <种类> <倍数>"。
- lb（T20-1）：160 个测试实例 × (rate + LP_fr、LP_fr^Rw 各 6 个 H)；超时 = 2 小时 + 10 分钟（单个 LP 限时 1 小时，另有索引核对与对偶证书）；
- opp（T20-2）：7 个对手——pb（1 次）、gphh（5 条规则）、esdw / rk / d5 / rddhand（种子 0–4）、rhmip（1 次）；
- ours（T20-3，等审核方说"开始"再派）：ours（种子 0–4，160 个）；消融（种子 0–2）noconn、shorteval（m = 4–10），conn_only（1 次，160 个）；
  窗口扫描 ours_w0.25 / ours_w1 / ours_w2（m = 4、6，种子 0–2）。
作业超时（防卡死）= 2·预算 + 30 分钟。长作业排前面（按预算从大到小），同预算内按方法、实例、种子。
可选第 3 个参数：只取这些 m（逗号分隔，例如 2,3,4,6），只作用于 job 行（lb 行总是全部 160 个实例）。
用法：python3 code/t20_make_jobs.py <lb|opp|lbopp|ours> <输出文件> [m 列表]"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import test_names  # noqa: E402
from t20_job import TLIM  # noqa: E402

MULTS = (0.5, 1, 1.5, 2, 3, 4)


def m_of(n):
    """实例名里的 m。"""
    return int(n[1:3])


def lb_jobs():
    """T20-1：大实例在前。"""
    out = []
    for n in sorted(test_names(), key=lambda n: (-m_of(n), n)):
        out.append((-1, f"{2 * 3600 + 600} lb {n} rate 0"))
        for kind in ("fr", "frRw"):
            for mu in MULTS: out.append((-1, f"{2 * 3600 + 600} lb {n} {kind} {mu:g}"))
    return [x[1] for x in out]


def job_lines(stage, plan):
    """plan = [(方法, 实例列表, 种子列表)] → 按预算从大到小排好的作业行。"""
    rows = []
    for meth, insts, seeds in plan:
        for n in insts:
            TL = TLIM[m_of(n)]
            for s in seeds: rows.append((-TL, meth, n, s, f"{2 * TL + 1800} job {stage} {meth} {n} {s}"))
    rows.sort(key=lambda r: r[:4])
    return [r[4] for r in rows]


def main():
    """按种类写作业表。"""
    kind, out = sys.argv[1], sys.argv[2]
    ms = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else None
    keep = lambda L: [n for n in L if ms is None or m_of(n) in ms]
    allt = keep(test_names()); mid = keep(test_names([4, 6, 8, 10])); small = keep(test_names([4, 6]))
    lines = []
    if kind in ("lb", "lbopp"): lines += lb_jobs()
    if kind in ("opp", "lbopp"):
        lines += job_lines("opp", [("pb", allt, [0]), ("gphh", allt, range(5)), ("esdw", allt, range(5)), ("rk", allt, range(5)),
                                   ("d5", allt, range(5)), ("rddhand", allt, range(5)), ("rhmip", allt, [0])])
    if kind == "ours":
        lines += job_lines("ours", [("ours", allt, range(5)), ("conn_only", allt, [0]), ("noconn", mid, range(3)), ("shorteval", mid, range(3)),
                                    ("ours_w0.25", small, range(3)), ("ours_w1", small, range(3)), ("ours_w2", small, range(3))])
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w") as f: f.write("\n".join(lines) + "\n")
    print(kind, len(lines), "行 →", out)


if __name__ == "__main__":
    main()
