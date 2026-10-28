"""T3b 汇总：不缩放的 R_LP（_u）与原版（按 N_LP 缩放）、_c（按 ⌈N_LP⌉ 缩放）逐实例对比。
输入：results/T3/jobs（原版、_c）、results/T3b/jobs（_u）、results/T2_bounds.jsonl（4S 的 LP 值）。
输出：results/T3b_rlp_u.jsonl（_u 的逐作业结果，字段同 T3_rlp.jsonl）和 markdown 表：
  1) 按 (m, 族) 的平均资源与离 ⌈LP⌉ 的差距（三种缩放 × {fresh_rlp, fresh_noreset}，以及 steady 对 steady_u）；
  2) 逐实例胜负：_u 对原版、_u 对 _c、fresh_rlp_u 对 fresh_noreset_u（更少 / 持平 / 更多）；
  3) 可行性：严格证明数、只靠 60 万格复核的个数、失败数。
资源只计可行的结果（严格证明或 60 万格复核 0 缺岗）。
用法：python3 code/t3b_summary.py"""
import os, sys, json, math
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402
import t3_summary  # noqa: E402

BASES = ["fresh_rlp", "fresh_noreset", "steady"]
SUF = ["", "_c", "_u"]


def feasible_R(r):
    """可行时返回资源数，否则 None。"""
    return r["R"] if r["R"] is not None and (r["gaps600k"] or 0) == 0 else None


def load():
    """返回 ({(实例, 方法): 行}, _u 的行列表)。"""
    a = t3_summary.rows(); u = t3_summary.rows(os.path.join(ROOT, "results", "T3b", "jobs"))
    return {(r["inst"], r["method"]): r for r in a + u}, u


def mean(xs):
    """忽略 None 的平均。"""
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def wlt(pairs):
    """pairs: [(a, b)]，返回 "a 更少 / 持平 / a 更多" 的计数和 a−b 的平均。"""
    pairs = [(x, y) for x, y in pairs if x is not None and y is not None]
    return (sum(x < y for x, y in pairs), sum(x == y for x, y in pairs), sum(x > y for x, y in pairs),
            mean([x - y for x, y in pairs]), len(pairs))


def tables(D, U, B):
    """生成三张表。"""
    insts = sorted({r["inst"] for r in U}, key=lambda n: (int(n[1:3]), n))
    groups = defaultdict(list)
    for n in insts: groups[(int(n[1:3]), n.split("-")[1])].append(n)
    fm = lambda x, nd=2: "—" if x is None else f"{x:.{nd}f}"
    ceil = lambda v: None if v is None else math.ceil(v - 1e-6)
    L = ["### 平均资源（离 4S 下界的差，%）\n",
         "| m | 族 | fresh_rlp | _c | _u | fresh_noreset | _c | _u | ⌈LP_fr^R⌉ | ⌈LP_fr⌉ | steady | steady_u |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for key in sorted(groups):
        ns = groups[key]; cells = []
        for base, lpk in (("fresh_rlp", "LP_frR"), ("fresh_noreset", "LP_fr")):
            for suf in SUF:
                Rs = [feasible_R(D[(n, base + suf)]) if (n, base + suf) in D else None for n in ns]
                gp = [100 * (R - ceil(B[n][lpk])) / ceil(B[n][lpk]) for n, R in zip(ns, Rs) if R is not None and n in B and B[n].get(lpk) is not None]
                cells.append(f"{fm(mean(Rs), 1)}（{fm(mean(gp), 1)}%）")
        lb = [fm(mean([ceil(B[n]["LP_frR"]) for n in ns]), 1), fm(mean([ceil(B[n]["LP_fr"]) for n in ns]), 1)]
        st = [fm(mean([feasible_R(D[(n, m)]) for n in ns if (n, m) in D]), 1) for m in ("steady", "steady_u")]
        L.append(f"| {key[0]} | {key[1]} | " + " | ".join(cells + lb + st) + " |")
    L += ["", "### 逐实例对比（a 比 b 少 / 持平 / 多，括号内 a − b 的平均）\n",
          "| m | 族 | rlp_u 对 rlp | rlp_u 对 rlp_c | noreset_u 对 noreset | noreset_u 对 noreset_c | **rlp_u 对 noreset_u** | rlp 对 noreset（原版） | steady_u 对 steady |",
          "|---|---|---|---|---|---|---|---|---|"]
    comps = [("fresh_rlp_u", "fresh_rlp"), ("fresh_rlp_u", "fresh_rlp_c"), ("fresh_noreset_u", "fresh_noreset"),
             ("fresh_noreset_u", "fresh_noreset_c"), ("fresh_rlp_u", "fresh_noreset_u"), ("fresh_rlp", "fresh_noreset"), ("steady_u", "steady")]
    tot = defaultdict(list)
    for key in sorted(groups):
        cells = []
        for a, b in comps:
            pr = [(feasible_R(D[(n, a)]) if (n, a) in D else None, feasible_R(D[(n, b)]) if (n, b) in D else None) for n in groups[key]]
            tot[(a, b)] += pr
            w, t, l, d, _ = wlt(pr); cells.append(f"{w}/{t}/{l}（{fm(d)}）")
        L.append(f"| {key[0]} | {key[1]} | " + " | ".join(cells) + " |")
    cells = []
    for a, b in comps:
        w, t, l, d, _ = wlt(tot[(a, b)]); cells.append(f"{w}/{t}/{l}（{fm(d)}）")
    L.append("| 合计 | | " + " | ".join(cells) + " |")
    L += ["", "### 可行性（_u）\n", "| 方法 | 作业数 | 严格证明 | 只靠 60 万格复核 | 不可行或失败 | 起点倍数 > 1 |", "|---|---|---|---|---|---|"]
    for m in ("fresh_rlp_u", "fresh_noreset_u", "steady_u"):
        rs = [r for r in U if r["method"] == m]
        pv = sum(bool(r["proven"]) for r in rs); ok = sum(feasible_R(r) is not None for r in rs)
        L.append(f"| {m} | {len(rs)} | {pv} | {ok - pv} | {len(rs) - ok} | {sum((r['start_factor'] or 100) > 100 for r in rs)} |")
    return "\n".join(L)


if __name__ == "__main__":
    D, U = load(); B = t3_summary.load_bounds()
    outp = os.path.join(ROOT, "results", "T3b_rlp_u.jsonl")
    with open(outp + ".tmp", "w") as f:
        for r in U: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(outp + ".tmp", outp)
    print(tables(D, U, B))
