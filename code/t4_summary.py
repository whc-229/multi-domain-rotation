"""T4 汇总：m ∈ {2, 3} 训练实例上的上下界与差距拆分。输出 results/T4_gap.jsonl 和 markdown 表。

两个场景分开算（全新启动）：
- 不允许提前大修：上界只取 resets == 0 的已证明解；下界 = max(⌈LP_fr(4S)⌉, ILP 下界 r0)。R_LP 取 fresh_noreset。
- 允许提前大修：上界取全部已证明解；下界 = max(⌈LP_fr^R(4S)⌉, ILP 下界 r1)。R_LP 取 fresh_rlp。
已证明 = 仿真器判出循环（环内无缺岗）且精确重放通过（proven == 1）。稳态场景（steady）的解不是全新启动的解，不参与。
上界来源：T3 的 R_LP 各变体与 ES（10 代 × 3 种子），T4 的 esub（30 代 × 3 种子）。
差距拆分：R_LP − 下界 = (R_LP − 最好上界) + (最好上界 − 下界)；前一项是规则相对已知最好排班的损失，后一项是上下界之间尚未确定的部分。
用法：python3 code/t4_summary.py"""
import os, sys, json, glob, math
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402

FRESH_KINDS = ["fresh_rlp", "fresh_noreset", "fresh_rlp_c", "fresh_noreset_c", "es_curdef", "es_lpgplan", "esub_curdef", "esub_lpgplan"]


def load_jobs(pattern):
    """读取一组结果 json。"""
    return [json.loads(open(p).read()) for p in sorted(glob.glob(pattern))]


def lp_bounds():
    """T2 的 4S LP 值：{实例: {"fr": 值, "frR": 值}}，只取求解成功的。"""
    b = defaultdict(dict)
    for d in load_jobs(os.path.join(ROOT, "results", "T2", "jobs", "m0[23]-*__fr*__m4.json")):
        if d.get("obj") is not None: b[d["inst"]][d["kind"]] = d["obj"]
    return b


def ilp_bounds():
    """T4 的 ILP 下界：{实例: {0: (LB, 状态, 用时), 1: …}}。"""
    b = defaultdict(dict)
    for d in load_jobs(os.path.join(ROOT, "results", "T4", "jobs", "*__lb__*.json")):
        b[d["inst"]][d["reset"]] = (d.get("LB"), d.get("status"), d.get("wall"))
    return b


def solutions():
    """所有全新启动的已证明解：{实例: [(方法, R, resets)]}。"""
    s = defaultdict(list)
    for pat in [("T3", "m0[23]-*.json"), ("T4", "m0[23]-*__esub_*.json")]:
        for d in load_jobs(os.path.join(ROOT, "results", pat[0], "jobs", pat[1])):
            if d.get("kind") not in FRESH_KINDS or d.get("proven") != 1 or d.get("R") is None: continue
            s[d["inst"]].append((d["kind"], d["R"], d.get("resets", 0)))
    return s


def per_instance():
    """逐实例、逐场景算上下界与差距拆分。"""
    LP, IL, SOL = lp_bounds(), ilp_bounds(), solutions()
    out = []
    for inst in sorted(set(LP) | set(SOL)):
        m = int(inst[1:3]); fam = inst.split("-")[1]
        for rs, lpk, rlpk in [(0, "fr", "fresh_noreset"), (1, "frR", "fresh_rlp")]:
            cand = [(k, R) for (k, R, nres) in SOL[inst] if rs == 1 or nres == 0]
            ub = min((R for _, R in cand), default=None)
            ub_src = sorted({k for k, R in cand if R == ub})
            lp = LP[inst].get(lpk); lp_c = math.ceil(lp - 1e-6) if lp is not None else None
            il = IL[inst].get(rs, (None, None, None))
            lbs = [x for x in (lp_c, il[0]) if x is not None]
            lb = max(lbs) if lbs else None
            rlp = next((R for (k, R, _) in SOL[inst] if k == rlpk), None)
            es_best = min((R for (k, R) in cand if k.startswith("es")), default=None)
            out.append(dict(inst=inst, m=m, family=fam, reset=rs, LP=lp, LP_ceil=lp_c, ILP_LB=il[0], ILP_status=il[1],
                            LB=lb, UB=ub, UB_src=ub_src, R_LP=rlp, ES_best=es_best,
                            rlp_minus_ub=None if None in (rlp, ub) else rlp - ub,
                            ub_minus_lb=None if None in (ub, lb) else ub - lb))
    return out


def mean(xs):
    """忽略 None 的平均。"""
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def table(rows):
    """按 (m, 族, 场景) 汇总的 markdown 表。"""
    g = defaultdict(list)
    for r in rows: g[(r["m"], r["family"], r["reset"])].append(r)
    fm = lambda x, nd=2: "—" if x is None else f"{x:.{nd}f}"
    L = ["| m | 族 | 提前大修 | 实例数 | ⌈LP⌉ | ILP 下界 | 下界 | 最好上界 | R_LP | R_LP−上界 | 上界−下界 | 上界−下界=0 的实例 | ILP 证到最优 |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for key in sorted(g):
        rs = g[key]
        z = sum(1 for r in rs if r["ub_minus_lb"] == 0)
        opt = sum(1 for r in rs if r["ILP_status"] == "Optimal")
        L.append(f"| {key[0]} | {key[1]} | {'允许' if key[2] else '不允许'} | {len(rs)} | {fm(mean([r['LP_ceil'] for r in rs]))} | "
                 f"{fm(mean([r['ILP_LB'] for r in rs]))} | {fm(mean([r['LB'] for r in rs]))} | {fm(mean([r['UB'] for r in rs]))} | "
                 f"{fm(mean([r['R_LP'] for r in rs]))} | {fm(mean([r['rlp_minus_ub'] for r in rs]))} | {fm(mean([r['ub_minus_lb'] for r in rs]))} | "
                 f"{z} | {opt} |")
    return "\n".join(L)


def detail(rows):
    """逐实例明细表。"""
    L = ["| 实例 | 提前大修 | LP | ILP 下界（状态） | 下界 | 最好上界（来源） | R_LP |", "|---|---|---|---|---|---|---|"]
    for r in rows:
        lp_s = "—" if r["LP"] is None else f"{r['LP']:.2f}"
        L.append(f"| {r['inst']} | {'允许' if r['reset'] else '不允许'} | {lp_s} | "
                 f"{r['ILP_LB'] if r['ILP_LB'] is not None else '—'}（{r['ILP_status'] or '—'}） | {r['LB']} | {r['UB']}（{', '.join(r['UB_src'])}） | {r['R_LP']} |")
    return "\n".join(L)


if __name__ == "__main__":
    rows = per_instance()
    p = os.path.join(ROOT, "results", "T4_gap.jsonl")
    with open(p + ".tmp", "w") as f:
        for r in rows: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(p + ".tmp", p)
    print(table(rows)); print(); print(detail(rows))
