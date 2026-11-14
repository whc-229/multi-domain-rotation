"""T20 长期可行的两种依据按 (方法, m) 计数（11-14 审核：区分"判出循环（证明）"与"只有 60 万格复核"）。
作业集合同 code/t20_paper_tables.py（opp、ours；opp_e / ours_e 里 (e) 实际触发的替换同名作业；只算 (a) 口径可行的）。
字段（rost5 judge 的输出）：
- 搜索型方法：作业内 judge 重判最终解的结果 rec["judge"]（check_judge = 1 表示与搜索时的记录逐字段相同）；
- shorteval：最终解本身就是搜索后按长期口径 judge（或 judge + rescue）的结果，取作业记录的同名字段；
- 确定性方法（pb、gphh、conn_only、rhmip）：没有单独重判，取定台数输出的同名字段（size 的判定与 judge 是同一段代码）。
判出循环 = verify_cert = 1 且 exact = 1（60 万格复核中判出循环，且循环经精确核对）；只有 60 万格复核 = 可行但不满足前者（复核 0 缺岗、没判出循环）。
输出 results/T20/cert_counts.md（markdown 表）并打印。用法：python3 code/t20_cert_counts.py"""
import os, sys, json, glob
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402

T20 = os.path.join(ROOT, "results", "T20")
ORDER = ["ours", "pb", "gphh", "esdw", "rk", "d5", "rddhand", "rhmip", "conn_only", "noconn", "pbnoconn", "shorteval", "ours_w0.25", "ours_w1", "ours_w2"]
MS = [2, 3, 4, 6, 8, 10, 15, 20]


def collect():
    """与论文表格同一作业集合：(方法, 实例, 种子键) → 记录。"""
    jobs = {}
    for d in ("opp", "ours"):
        for f in glob.glob(os.path.join(T20, d, "jobs", "*.json")):
            meth, inst, key = os.path.basename(f)[:-5].split("__"); jobs[(meth, inst, key)] = json.load(open(f))
    for d in ("opp_e", "ours_e"):
        for f in glob.glob(os.path.join(T20, d, "jobs", "*.json")):
            meth, inst, key = os.path.basename(f)[:-5].split("__"); rec = json.load(open(f))
            if rec.get("e_fixed", -1) >= 0: jobs[(meth, inst, key)] = rec
    return jobs


def basis(rec):
    """一个可行作业的依据：("证明" | "复核", 是否来自 judge 重判)；不可行返回 None。"""
    if rec.get("status") not in ("cert", "uncert"): return None
    j = rec.get("judge")
    src = j if (j is not None and rec.get("check_judge") == 1) else rec
    return ("证明" if src.get("verify_cert") == 1 and src.get("exact") == 1 else "复核"), j is not None


def main():
    """计数并写表。"""
    cnt = defaultdict(lambda: defaultdict(int))
    for (meth, inst, key), rec in collect().items():
        b = basis(rec); m = int(inst[1:3])
        cnt[meth][(m, "总数")] += 1
        if b: cnt[meth][(m, b[0])] += 1; cnt[meth][(m, "重判")] += int(b[1])
    lines = ["| 方法 | 字段来源 | " + " | ".join(f"m = {m}" for m in MS) + " | 合计 |", "|---|---|" + "---|" * (len(MS) + 1)]
    for meth in [x for x in ORDER if x in cnt]:
        C = cnt[meth]; tot = [0, 0, 0]
        cells = []
        for m in MS:
            n = C.get((m, "总数"), 0)
            if not n: cells.append("—"); continue
            p, v = C.get((m, "证明"), 0), C.get((m, "复核"), 0); tot[0] += p; tot[1] += v; tot[2] += n
            cells.append(f"{p} / {v}（{n}）")
        nre = sum(C.get((m, "重判"), 0) for m in MS)
        src = "judge 重判" if nre else ("作业记录（长期 judge）" if meth == "shorteval" else "定台数输出")
        lines.append(f"| {meth} | {src} | " + " | ".join(cells) + f" | {tot[0]} / {tot[1]}（{tot[2]}） |")
    txt = "格内为 判出循环 / 只有 60 万格复核（运行数），只计可行的运行。\n\n" + "\n".join(lines) + "\n"
    open(os.path.join(T20, "cert_counts.md"), "w").write(txt); print(txt)


if __name__ == "__main__":
    main()
