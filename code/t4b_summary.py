"""T4b 汇总：长预算 ES（100 代 × 5 个种子）能把 m ∈ {2, 3} 的上界压到多低。
输入：results/T4b/jobs（es100），results/T4/jobs（esub，30 代），results/T3/jobs（es，10 代），results/T4_exact.jsonl（T4 的下界与上界）。
输出：results/T4b_es100.jsonl（逐作业）和 markdown 表：
  1) 核对：种子 0–2 的 es100 逐代最好值（hist）前 30 代与 esub 相同、前 10 代与 T3 的 es 相同（同一种子的随机数流与代数无关）；
  2) 10 / 30 / 100 代的最好值曲线：每个 (实例, 配置) 取种子 0–2 的最好可行值，再按 (m, 族) 平均；另列 100 代 5 个种子；
  3) 上界−下界差距的变化：T4 的最好上界与加入 es100 后的最好上界，各场景（不允许 / 允许提前大修）分开。
可行 = 严格证明或 60 万格复核 0 缺岗；不允许提前大修的场景只取提前大修次数为 0 的解（与 t4_summary 相同）。
hist 是搜索时的适应度（6 万格证书口径下的资源数），不是最终复核后的值；曲线用最终复核后的值，hist 只用于核对。
用法：python3 code/t4b_summary.py"""
import os, sys, json, glob
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402


def load(pattern):
    """读一组结果 json，返回 {(实例, 配置, 种子): 记录}。配置 = curdef / lpgplan。"""
    out = {}
    for p in glob.glob(os.path.join(ROOT, pattern)):
        d = json.loads(open(p).read()); out[(d["inst"], d["kind"].split("_", 1)[1], d["seed"])] = d
    return out


def feas(d, reset):
    """该解在场景 reset 下可行时返回资源数，否则 None。"""
    if d is None or d.get("R") is None or d["R"] < 0: return None
    if not (d.get("proven") == 1 or d.get("verify_gaps") == 0): return None
    if reset == 0 and d.get("resets", 0) != 0: return None
    return d["R"]


def mean(xs):
    """忽略 None 的平均。"""
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def best(vals):
    """忽略 None 的最小值。"""
    vals = [v for v in vals if v is not None]
    return min(vals) if vals else None


def main():
    E10, E30, E100 = load("results/T3/jobs/m0[23]-*__es_*__s*.json"), load("results/T4/jobs/m0[23]-*__esub_*__s*.json"), load("results/T4b/jobs/*__es100_*.json")
    T4 = {(r["inst"], r["reset"]): r for r in map(json.loads, open(os.path.join(ROOT, "results", "T4_exact.jsonl")))}
    # 1) 核对
    same = diff = 0
    for (n, c, s), d in E100.items():
        if s > 2: continue
        a, b = E30.get((n, c, s)), E10.get((n, c, s))
        ok = a is not None and b is not None and d["hist"][:31] == a["hist"] and d["hist"][:11] == b["hist"]
        same += ok; diff += (not ok)
    rows = []
    for (n, c, s), d in sorted(E100.items()):
        rows.append(dict(inst=n, m=d["m"], family=d["family"], cfg=c, seed=s, R=d.get("R"), resets=d.get("resets"), proven=d.get("proven"),
                         verify_gaps=d.get("verify_gaps"), hist10=d["hist"][10], hist30=d["hist"][30], hist100=d["hist"][-1], wall=d.get("wall")))
    with open(os.path.join(ROOT, "results", "T4b_es100.jsonl.tmp"), "w") as f:
        for r in rows: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(os.path.join(ROOT, "results", "T4b_es100.jsonl.tmp"), os.path.join(ROOT, "results", "T4b_es100.jsonl"))
    insts = sorted({k[0] for k in E100})
    fm = lambda x, nd=2: "—" if x is None else f"{x:.{nd}f}"
    L = [f"核对：种子 0–2 的 es100 前 30 代 hist 与 esub 相同、前 10 代与 T3 es 相同：{same} 个相同，{diff} 个不同。", "",
         "### 最好可行值随代数（每个 (实例, 配置) 取种子的最好值，再按 (m, 族) 平均；只计 100 代 5 个种子都跑完的实例）\n",
         "| m | 族 | 场景 | 配置 | 实例数 | 10 代（3 种子） | 30 代（3 种子） | 100 代（3 种子） | 100 代（5 种子） | 100 代 5 种子 比 30 代 少 | 有改进的实例 |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    full = [n for n in insts if all((n, c, s) in E100 for c in ("curdef", "lpgplan") for s in range(5))]
    for (m, fam) in sorted({(int(n[1:3]), n.split("-")[1]) for n in full}):
        ns = [n for n in full if int(n[1:3]) == m and n.split("-")[1] == fam]
        for rs in (0, 1):
            for c in ("curdef", "lpgplan"):
                v10 = [best(feas(E10.get((n, c, s)), rs) for s in range(3)) for n in ns]
                v30 = [best(feas(E30.get((n, c, s)), rs) for s in range(3)) for n in ns]
                v100 = [best(feas(E100.get((n, c, s)), rs) for s in range(3)) for n in ns]
                v105 = [best(feas(E100.get((n, c, s)), rs) for s in range(5)) for n in ns]
                dd = [a - b for a, b in zip(v30, v105) if a is not None and b is not None]
                if mean(v105) is None: continue        # 不允许提前大修时 lpgplan 的解都带提前大修，整行为空，不列
                L.append(f"| {m} | {fam} | {'允许' if rs else '不允许'} | {c} | {len(ns)} | {fm(mean(v10))} | {fm(mean(v30))} | {fm(mean(v100))} | {fm(mean(v105))} | "
                         f"{fm(mean(dd))} | {sum(x > 0 for x in dd)} / {len(dd)} |")
    L += ["", "### 上界−下界差距的变化（T4 的最好上界 → 加入 es100 后）\n",
          "| m | 族 | 场景 | 实例数 | 下界 | T4 最好上界 | 新最好上界 | 上界−下界 | % | 新 上界−下界 | % | 上界降低的实例 | 闭合的实例 |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    out = []
    for (m, fam) in sorted({(int(n[1:3]), n.split("-")[1]) for n in full}):
        ns = [n for n in full if int(n[1:3]) == m and n.split("-")[1] == fam]
        for rs in (0, 1):
            g = []
            for n in ns:
                t = T4[(n, rs)]; lb, ub = t["LB"], t["UB"]
                nb = best([ub] + [feas(E100.get((n, c, s)), rs) for c in ("curdef", "lpgplan") for s in range(5)])
                g.append((lb, ub, nb)); out.append(dict(inst=n, reset=rs, LB=lb, UB_T4=ub, UB_new=nb))
            ok = [x for x in g if None not in x]
            L.append(f"| {m} | {fam} | {'允许' if rs else '不允许'} | {len(ns)} | {fm(mean([x[0] for x in ok]))} | {fm(mean([x[1] for x in ok]))} | {fm(mean([x[2] for x in ok]))} | "
                     f"{fm(mean([x[1] - x[0] for x in ok]))} | {fm(mean([100 * (x[1] - x[0]) / x[0] for x in ok]), 1)} | "
                     f"{fm(mean([x[2] - x[0] for x in ok]))} | {fm(mean([100 * (x[2] - x[0]) / x[0] for x in ok]), 1)} | "
                     f"{sum(x[2] < x[1] for x in ok)} | {sum(x[2] == x[0] for x in ok)} |")
    with open(os.path.join(ROOT, "results", "T4b_gap.jsonl.tmp"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(os.path.join(ROOT, "results", "T4b_gap.jsonl.tmp"), os.path.join(ROOT, "results", "T4b_gap.jsonl"))
    # 用时与可行性
    w = defaultdict(list); fe = defaultdict(lambda: [0, 0, 0])
    for r in rows:
        w[(r["m"], r["family"], r["cfg"])].append(r["wall"])
        a = fe[(r["m"], r["family"], r["cfg"])]; a[0] += 1; a[1] += (r["proven"] == 1); a[2] += (feas(E100[(r["inst"], r["cfg"], r["seed"])], 1) is None)
    L += ["", "### 用时与可行性（es100 每个作业）\n", "| m | 族 | 配置 | 作业 | 严格证明 | 不可行 | 平均用时 s | 最长 s |", "|---|---|---|---|---|---|---|---|"]
    for k in sorted(w):
        L.append(f"| {k[0]} | {k[1]} | {k[2]} | {fe[k][0]} | {fe[k][1]} | {fe[k][2]} | {mean(w[k]):.0f} | {max(w[k]):.0f} |")
    print("\n".join(L))


if __name__ == "__main__":
    main()
