"""T3 汇总：合并 results/T3/jobs/*.json → results/T3_rlp.jsonl，并按 (m, 族) 打印各方法离 ⌈LP_fr⌉、⌈LP_fr^R⌉ 的差距。
每行：实例、方法、场景、种子、资源数、机器人数、工具数、提前大修次数、是否严格证明、周期长度、60 万格缺岗数、用时、定台数起点倍数。
60 万格缺岗数：证明了的记 0；判不出循环的取复核 60 万格内的缺岗数；判出循环但环内有缺岗的取复核缺岗数（长期不可行）。
ES 每个 (实例, 配置) 有 3 个种子，表里先对种子取平均。
用法：python3 code/t3_summary.py"""
import os, sys, json, glob, math
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402

METHODS = ["fresh_rlp", "fresh_noreset", "fresh_rlp_c", "fresh_noreset_c", "es_curdef", "es_lpgplan", "steady"]


def load_bounds():
    b = {}
    p = os.path.join(ROOT, "results", "T2_bounds.jsonl")
    if os.path.exists(p):
        for l in open(p):
            d = json.loads(l)
            if d["mult"] == 4.0: b[d["inst"]] = d
    return b


def rows(jdir=None):
    """读入一个作业目录（默认 results/T3/jobs）的全部结果，每个作业一行。"""
    out = []
    jdir = jdir or os.path.join(ROOT, "results", "T3", "jobs")
    for p in sorted(glob.glob(os.path.join(jdir, "*.json"))):
        d = json.loads(open(p).read())
        ok = d.get("R", -1) is not None and d.get("R", -1) >= 0
        vc = d.get("verify_cert")
        gaps600 = 0 if (vc == 1 and d.get("exact") == 1) else d.get("verify_gaps")
        out.append(dict(inst=d["inst"], m=d["m"], family=d["family"], method=d["kind"],
                        scenario="steady" if d["kind"] == "steady" else "fresh", seed=d.get("seed"),
                        R=d["R"] if ok else None, robots=d.get("robots", sum(d.get("Nk", [])) if ok else None),
                        tools=d.get("tools", sum(d.get("Mj", [])) if ok else None),
                        resets=d.get("resets"), proven=d.get("proven"), period=d.get("verify_period"),
                        gaps600k=gaps600, status=d.get("status"), start_factor=d.get("start_factor"),
                        sec=d.get("sec", d.get("wall")), wall=d.get("wall")))
    return out


def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def table(rs, B):
    g = defaultdict(lambda: defaultdict(list))
    for r in rs:
        key = (r["m"], r["family"])
        g[key][(r["method"], r["inst"])].append(r)
    lines = ["| m | 族 | 方法 | 完成/证明/失败 | 平均资源 | 离 ⌈LP_fr⌉ | % | 离 ⌈LP_fr^R⌉ | % | 平均用时 s |", "|---|---|---|---|---|---|---|---|---|---|"]
    for (m, fam) in sorted(g):
        for meth in METHODS:
            per = [v for (mm, n), v in g[(m, fam)].items() if mm == meth]
            if not per: continue
            Rm, d1, p1, d2, p2, secs = [], [], [], [], [], []
            ndone = nprov = nfail = 0
            for runs in per:
                for r in runs:
                    ndone += 1; nprov += bool(r["proven"]); nfail += (r["R"] is None or (r["gaps600k"] or 0) > 0)
                    secs.append(r["wall"])
                good = [r["R"] for r in runs if r["R"] is not None and (r["gaps600k"] or 0) == 0]
                if not good: continue
                R = sum(good) / len(good); Rm.append(R)
                b = B.get(runs[0]["inst"])
                if b is None: continue
                if meth == "steady":
                    ref = b.get("R_rate_int")
                    if ref: d1.append(R - ref); p1.append(100 * (R - ref) / ref)
                    continue
                if b.get("LP_fr") is not None:
                    c = math.ceil(b["LP_fr"] - 1e-6); d1.append(R - c); p1.append(100 * (R - c) / c)
                if b.get("LP_frR") is not None:
                    c = math.ceil(b["LP_frR"] - 1e-6); d2.append(R - c); p2.append(100 * (R - c) / c)
            fm = lambda x, nd=2: "—" if x is None else f"{x:.{nd}f}"
            label = meth + ("（对 R_rate_int）" if meth == "steady" else "")
            lines.append(f"| {m} | {fam} | {label} | {ndone}/{nprov}/{nfail} | {fm(mean(Rm),1)} | {fm(mean(d1))} | {fm(mean(p1),1)} | {fm(mean(d2))} | {fm(mean(p2),1)} | {fm(mean(secs),0)} |")
    return "\n".join(lines)


if __name__ == "__main__":
    rs = rows(); B = load_bounds()
    outp = os.path.join(ROOT, "results", "T3_rlp.jsonl")
    with open(outp + ".tmp", "w") as f:
        for r in rs: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(outp + ".tmp", outp)
    print(table(rs, B))
