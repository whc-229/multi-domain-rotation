"""T2 汇总：合并 results/T2/jobs/*.json → results/T2_bounds.jsonl（每个 (实例, H 倍数) 一行），并按 (m, 族) 打印 Markdown 表。
每行：实例、m、族、H 倍数、H、R_rate、R_rate_int、LP_fr、LP_fr^R、用时、是否解到最优、对偶证书最大违反量；
4S 行另有 "收敛" 字段 = LP(4S) − LP(2S)（两者都解出时）。
用法：python3 code/t2_summary.py [结果目录] [输出 jsonl]"""
import os, sys, json, glob, math
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402


def load(jdir):
    rec = {}
    for p in glob.glob(os.path.join(jdir, "*.json")):
        d = json.loads(open(p).read()); rec[(d["inst"], d["kind"], d["mult"])] = d
    return rec


def merge(rec):
    insts = sorted({k[0] for k in rec})
    rows = []
    for n in insts:
        rate = rec.get((n, "rate", 0.0))
        for mult in (2.0, 4.0):
            fr, frR = rec.get((n, "fr", mult)), rec.get((n, "frR", mult))
            if fr is None and frR is None: continue
            base = fr or frR
            row = dict(inst=n, m=base["m"], family=base["family"], mult=mult, H=base["H"], S=base["S"],
                       R_rate=rate and rate["R_rate"], R_rate_int=rate and rate["R_rate_int"])
            for tag, d in (("fr", fr), ("frR", frR)):
                if d is None: row[f"LP_{tag}"] = None; row[f"{tag}_done"] = False; continue
                ok = d.get("obj") is not None and d.get("lp_status") == "Optimal"
                row[f"LP_{tag}"] = d.get("obj_full") if ok else None
                row[f"{tag}_done"] = ok
                row[f"{tag}_sec"] = d.get("wall")
                row[f"{tag}_cert_viol"] = max(d.get("path_viol", math.nan), d.get("tool_viol", math.nan), d.get("wait_viol", math.nan), d.get("M_viol", math.nan)) if ok else None
                row[f"{tag}_index_mismatch"] = d.get("index_mismatch")
                row[f"{tag}_N"] = d.get("N"); row[f"{tag}_M"] = d.get("M")
                row[f"{tag}_resets_flow"] = d.get("resets_flow")
            rows.append(row)
    by = {(r["inst"], r["mult"]): r for r in rows}
    for r in rows:
        if r["mult"] == 4.0 and (r["inst"], 2.0) in by:
            o = by[(r["inst"], 2.0)]
            for tag in ("fr", "frR"):
                a, b = r[f"LP_{tag}"], o[f"LP_{tag}"]
                r[f"conv_{tag}"] = (a - b) if (a is not None and b is not None) else None
    return rows


def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def f(x, nd=2):
    return "—" if x is None else f"{x:.{nd}f}"


def tables(rows, mult=4.0):
    out = []
    g = defaultdict(list)
    for r in rows:
        if r["mult"] == mult: g[(r["m"], r["family"])].append(r)
    out.append(f"### 下界（H = {mult:g}·S_max，各组 6 个训练实例的平均）\n")
    out.append("| m | 族 | 完成 fr/frR | R_rate | R_rate_int | LP_fr | LP_fr^R | 启动代价 LP_fr−R_rate | 启动代价 % | 提前大修收益 LP_fr−LP_fr^R | 收益 % | ⌈LP_fr⌉−R_rate_int | ⌈LP_fr⌉−⌈LP_fr^R⌉ |")
    out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for (m, fam) in sorted(g):
        rs = g[(m, fam)]
        both = [r for r in rs if r["LP_fr"] is not None and r["LP_frR"] is not None and r["R_rate"] is not None]
        sc = [r["LP_fr"] - r["R_rate"] for r in both]; scp = [100 * (r["LP_fr"] - r["R_rate"]) / r["R_rate"] for r in both]
        gn = [r["LP_fr"] - r["LP_frR"] for r in both]; gnp = [100 * (r["LP_fr"] - r["LP_frR"]) / r["LP_fr"] for r in both]
        ci = [math.ceil(r["LP_fr"] - 1e-6) - r["R_rate_int"] for r in both]
        cg = [math.ceil(r["LP_fr"] - 1e-6) - math.ceil(r["LP_frR"] - 1e-6) for r in both]
        out.append(f"| {m} | {fam} | {sum(r['fr_done'] for r in rs)}/{sum(r['frR_done'] for r in rs)} | {f(mean([r['R_rate'] for r in rs]))} | {f(mean([r['R_rate_int'] for r in rs]),1)} | "
                   f"{f(mean([r['LP_fr'] for r in rs]))} | {f(mean([r['LP_frR'] for r in rs]))} | {f(mean(sc))} | {f(mean(scp),1)} | {f(mean(gn))} | {f(mean(gnp),1)} | {f(mean(ci),2)} | {f(mean(cg),2)} |")
    return "\n".join(out)


def time_table(rows):
    out = ["### 用时（秒，墙钟；每格为该组 fr / frR 的平均，括号内最大）\n",
           "| m | 族 | 2S fr | 2S frR | 4S fr | 4S frR |", "|---|---|---|---|---|---|"]
    g = defaultdict(list)
    for r in rows: g[(r["m"], r["family"])].append(r)
    for (m, fam) in sorted(g):
        cells = []
        for mult in (2.0, 4.0):
            for tag in ("fr", "frR"):
                xs = [r.get(f"{tag}_sec") for r in g[(m, fam)] if r["mult"] == mult and r.get(f"{tag}_sec") is not None]
                cells.append(f"{mean(xs):.0f} ({max(xs):.0f})" if xs else "—")
        out.append(f"| {m} | {fam} | " + " | ".join(cells) + " |")
    return "\n".join(out)


def conv_table(rows):
    out = ["### H 收敛：LP(4S) − LP(2S)（各组最小 / 最大）\n", "| m | 族 | fr | frR |", "|---|---|---|---|"]
    g = defaultdict(list)
    for r in rows:
        if r["mult"] == 4.0: g[(r["m"], r["family"])].append(r)
    for (m, fam) in sorted(g):
        cells = []
        for tag in ("fr", "frR"):
            xs = [r.get(f"conv_{tag}") for r in g[(m, fam)] if r.get(f"conv_{tag}") is not None]
            cells.append(f"{min(xs):+.3f} / {max(xs):+.3f}" if xs else "—")
        out.append(f"| {m} | {fam} | " + " | ".join(cells) + " |")
    return "\n".join(out)


def check_table(rows):
    v = [r.get(f"{t}_cert_viol") for r in rows for t in ("fr", "frR") if r.get(f"{t}_cert_viol") is not None]
    im = [r.get(f"{t}_index_mismatch") for r in rows for t in ("fr", "frR") if r.get(f"{t}_index_mismatch") is not None]
    return (f"- 对偶证书：{len(v)} 个 LP，最大违反量 {max(v):.1e}；超过 1e-6 的 {sum(x > 1e-6 for x in v)} 个。\n"
            f"- 索引重建比对：{len(im)} 个 LP，不一致非零元合计 {sum(im)}。") if v else "- 尚无结果"


if __name__ == "__main__":
    jdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results", "T2", "jobs")
    outp = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "results", "T2_bounds.jsonl")
    rows = merge(load(jdir))
    with open(outp + ".tmp", "w") as fo:
        for r in rows: fo.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(outp + ".tmp", outp)
    print(tables(rows, 4.0)); print(); print(tables(rows, 2.0)); print(); print(conv_table(rows)); print(); print(time_table(rows)); print(); print(check_table(rows))
