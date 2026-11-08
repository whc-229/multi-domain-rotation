"""T20-1 汇总（测试实例 160 个）：下界与时间界的 LP。
- 每个 (实例, 种类, 倍数) 取正式结果；正式结果没有解（超时）时取 __retry1（与 T17 汇总相同的取法）。
- LP 值 = 原始目标与证书对偶目标中较小者；证书违反量 = max(路径、工具、等待、M)（同 code/t17_summary.py）。
- 速率下界取 R_rate_g（网格口径 LP）与 R_rate_g_int（整数）。
- H₀ = LP(H) ≥ 99.5%·LP(4S) 的最小倍数（只在 4S 有解时算），按 m 取中位数。
输出 results/T20/lb_summary.jsonl（逐实例）与 markdown 表（打印）。
用法：python3 code/t20_lb_summary.py"""
import os, sys, json, math, statistics as st
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, test_names  # noqa: E402
from t17_summary import pick, lp_value, MULTS  # noqa: E402

J = os.path.join(ROOT, "results", "T20", "lb", "jobs")
MS = (2, 3, 4, 6, 8, 10, 15, 20)


def inst_row(n):
    """一个实例：速率下界、两种 LP 在各倍数下的值、来源（正式 / retry1）、证书违反量、H₀。"""
    rate = json.load(open(os.path.join(J, f"{n}__rate__m0.json")))
    row = dict(inst=n, m=int(n[1:3]), family=n.split("-")[1], R_rate_g=rate["R_rate_g"], R_rate_g_int=rate["R_rate_g_int"])
    for kind in ("fr", "frRw"):
        vals, src, viol = {}, {}, 0.0
        for mu in MULTS:
            d, b = pick(os.path.join(J, f"{n}__{kind}__m{mu:g}"))
            if d is None: vals[mu] = None; src[mu] = "无解"; continue
            v, vi = lp_value(d); vals[mu] = v; src[mu] = "retry1" if b.endswith("__retry1") else "正式"; viol = max(viol, vi)
        top = vals[4.0]
        h0 = None if top is None else next((mu for mu in MULTS if vals[mu] is not None and vals[mu] >= 0.995 * top), None)
        row[kind] = {f"{mu:g}": vals[mu] for mu in MULTS}; row[kind + "_src"] = {f"{mu:g}": src[mu] for mu in MULTS}
        row[kind + "_viol"] = viol; row[kind + "_H0"] = h0
    return row


def mean(xs):
    """均值（忽略 None）；全空返回 None。"""
    xs = [x for x in xs if x is not None]
    return st.mean(xs) if xs else None


def fmt(x, nd=1):
    """表格里的数字。"""
    return "—" if x is None else f"{x:.{nd}f}"


def main():
    """逐实例写 jsonl，打印四张表：速率下界、两种 LP 各倍数的平均值（共同实例）与缺失数、H₀ 中位数。"""
    rows = [inst_row(n) for n in test_names(list(MS))]
    with open(os.path.join(ROOT, "results", "T20", "lb_summary.jsonl"), "w") as f:
        for r in rows: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    by = defaultdict(list)
    for r in rows: by[r["m"]].append(r)
    print("| m | 实例 | R_rate_g | R_rate_g_int（整数） | 最大证书违反量 |")
    print("|---|---|---|---|---|")
    for m in MS:
        R = by[m]
        print(f"| {m} | {len(R)} | {fmt(mean(r['R_rate_g'] for r in R))} | {fmt(mean(r['R_rate_g_int'] for r in R))} | "
              f"{max(max(r['fr_viol'], r['frRw_viol']) for r in R):.1e} |")
    for kind, name in (("fr", "LP_fr(H)"), ("frRw", "LP_fr^Rw(H)")):
        print(f"\n{name}：只在 6 个倍数都有解的实例上取平均（不同倍数用同一组实例才可比）；最后两列为无解、取 retry1 的 (实例, 倍数) 个数")
        print("| m | 实例 | " + " | ".join(f"{mu:g}S" for mu in MULTS) + " | 无解 | retry1 |"); print("|---|---|" + "---|" * len(MULTS) + "---|---|")
        for m in MS:
            R = by[m]; C = [r for r in R if all(r[kind][f"{mu:g}"] is not None for mu in MULTS)]
            miss = sum(v == "无解" for r in R for v in r[kind + "_src"].values()); rt = sum(v == "retry1" for r in R for v in r[kind + "_src"].values())
            print(f"| {m} | {len(C)} | " + " | ".join(fmt(mean(r[kind][f'{mu:g}'] for r in C)) for mu in MULTS) + f" | {miss} | {rt} |")
    print("\nH₀（LP(H) ≥ 99.5%·LP(4S) 的最小倍数）中位数，括号为 4S 有解的实例数")
    print("| m | LP_fr | LP_fr^Rw |"); print("|---|---|---|")
    for m in MS:
        cs = []
        for kind in ("fr", "frRw"):
            h = [r[kind + "_H0"] for r in by[m] if r[kind + "_H0"] is not None]
            cs.append(f"{st.median(h):g} ({len(h)})" if h else "— (0)")
        print(f"| {m} | {cs[0]} | {cs[1]} |")


if __name__ == "__main__":
    main()
