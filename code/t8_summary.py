"""T8 汇总：RDD 可行性试点（全新启动、不提前大修）。输出 results/T8_rdd.jsonl 和 markdown 表（打印）。
- 可行 = 严格证明（判出循环、环内无缺岗、精确重放通过）或 60 万格复核 0 缺岗；不可行的结果不进平均。
- 默认规则本身（base）：同一作业里、用同样的定台数流程算出的结果（与种子无关）。
- ES 饱和值：T4b 汇总的不提前大修最好上界（UB_new，含 T3 / T4 / T4b 的全部 ES 与 R_LP 解），只对 m ≤ 3 有。
- ρ 扫描（sweep）：m = 2 的 12 个实例 × 3 种规则 × ρ ∈ {0.02, 0.1, 0.4} × 种子 0–2；按 (实例, 规则, 种子) 配对比较各 ρ 的资源数。
用法：python3 code/t8_summary.py [sweep | full <ρ>]"""
import os, sys, json, glob
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402

J = os.path.join(ROOT, "results", "T8", "jobs")
RULES = ("rlp", "es", "hand")


def load():
    """全部 T8 结果：{(实例, 规则, ρ, 种子): 记录}。"""
    out = {}
    for f in glob.glob(os.path.join(J, "*.json")):
        d = json.load(open(f)); out[(d["inst"], d["rule"], float(d["rho"]), int(d["seed"]))] = d
    return out


def fR(d):
    """RDD 最好解可行时返回资源数。"""
    return d["R"] if d and d.get("R") is not None and d["R"] >= 0 and (d.get("proven") == 1 or d.get("verify_gaps") == 0) else None


def fB(d):
    """默认规则本身可行时返回资源数。"""
    return d["base_R"] if d and d.get("base_R") is not None and d["base_R"] >= 0 and (d.get("base_proven") == 1 or d.get("base_verify_gaps") == 0) else None


def mean(xs):
    """忽略 None 的平均。"""
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def fm(x, nd=2):
    """格式化。"""
    return "—" if x is None else f"{x:.{nd}f}"


def sweep(D):
    """ρ 扫描表：每个 (规则, ρ) 的平均资源数（三个 ρ 都可行的配对上）、对默认规则的平均减少、可行率。"""
    rhos = (0.02, 0.1, 0.4)
    keys = sorted({(n, r, s) for (n, r, rho, s) in D if n.startswith("m02-") and s <= 2})
    L = ["| 规则 | ρ | 作业 | 可行 | 平均资源数（三个 ρ 都可行的配对） | 比默认规则少（平均） | 比默认规则少的作业 | 平均用时 s |", "|---|---|---|---|---|---|---|---|"]
    best = {}
    allm = defaultdict(list)
    for rule in RULES:
        ks = [k for k in keys if k[1] == rule]
        common = [k for k in ks if all(fR(D.get((k[0], rule, rho, k[2]))) is not None for rho in rhos)]
        for rho in rhos:
            ds = [D.get((n, rule, rho, s)) for (n, _, s) in ks]
            ds = [d for d in ds if d]
            m = mean([fR(D[(n, rule, rho, s)]) for (n, _, s) in common])
            red = mean([fB(d) - fR(d) for d in ds if fB(d) is not None and fR(d) is not None])
            allm[rho] += [fR(D[(n, rule, rho, s)]) for (n, _, s) in common]
            L.append(f"| {rule} | {rho:g} | {len(ds)} | {sum(fR(d) is not None for d in ds)} | {fm(m)}（{len(common)} 组） | {fm(red)} | "
                     f"{sum(1 for d in ds if fB(d) is not None and fR(d) is not None and fR(d) < fB(d))} | {fm(mean([d.get('wall') for d in ds]), 0)} |")
    L.append("")
    L.append("三种规则合在一起（各自三个 ρ 都可行的配对）：" + "，".join(f"ρ = {rho:g}：{fm(mean(allm[rho]), 3)}" for rho in rhos))
    return L


def full(D, rho):
    """全部实例在选定 ρ 下的结果：按 (m, 族, 规则) 平均、逐实例胜负、低于 ES 饱和值的实例清单。"""
    T4b = {r["inst"]: r for r in map(json.loads, open(os.path.join(ROOT, "results", "T4b_gap.jsonl"))) if r["reset"] == 0}
    lb = {}
    for r in map(json.loads, open(os.path.join(ROOT, "results", "T6_bounds.jsonl"))): lb[r["inst"]] = r["LP_fr"]
    recs = []
    insts = sorted({k[0] for k in D if k[2] == rho}, key=lambda n: (int(n[1:3]), n))
    for n in insts:
        for rule in RULES:
            ds = [D.get((n, rule, rho, s)) for s in range(5)]
            ds = [d for d in ds if d]
            if not ds: continue
            base = fB(ds[0]); vals = [fR(d) for d in ds]
            best = min([v for v in vals if v is not None], default=None)
            ub = T4b[n]["UB_new"] if n in T4b else None
            recs.append(dict(inst=n, m=int(n[1:3]), family=n.split("-")[1], rule=rule, rho=rho, seeds=len(ds), base=base,
                             base_status=ds[0]["base_status"], rdd=vals, rdd_best=best, rdd_mean=mean(vals),
                             es_sat=ub, LP_fr=lb.get(n), ndev=ds[0]["ndev"], K=ds[0]["K"],
                             n_active=[d.get("n_active") for d in ds], proven=[d.get("proven") for d in ds],
                             sec=[d.get("wall") for d in ds]))
    L = ["| m | 族 | 规则 | 实例 | 默认规则 | RDD（种子平均） | RDD（5 种子最好） | 比默认规则少（平均） | 少 / 同 / 多（种子平均对默认） | ES 饱和值 | RDD 最好 − ES 饱和值 | 离 ⌈LP_fr⌉ %（RDD 最好） | 平均用时 min |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    by = defaultdict(list)
    for r in recs: by[(r["m"], r["family"], r["rule"])].append(r)
    import math
    for (m, f, rule), rs in sorted(by.items()):
        ok = [r for r in rs if r["base"] is not None and r["rdd_mean"] is not None]
        d = [r["base"] - r["rdd_mean"] for r in ok]
        es = [r["rdd_best"] - r["es_sat"] for r in rs if r["es_sat"] is not None and r["rdd_best"] is not None]
        gap = [100 * (r["rdd_best"] - math.ceil(r["LP_fr"] - 1e-6)) / math.ceil(r["LP_fr"] - 1e-6) for r in rs if r["rdd_best"] is not None and r["LP_fr"]]
        L.append(f"| {m} | {f} | {rule} | {len(rs)} | {fm(mean([r['base'] for r in ok]))} | {fm(mean([r['rdd_mean'] for r in ok]))} | {fm(mean([r['rdd_best'] for r in ok]))} | "
                 f"{fm(mean(d))} | {sum(x > 1e-9 for x in d)} / {sum(abs(x) <= 1e-9 for x in d)} / {sum(x < -1e-9 for x in d)} | "
                 f"{fm(mean([r['es_sat'] for r in rs if r['es_sat'] is not None]))} | {fm(mean(es))} | {fm(mean(gap), 1)} | {fm(mean([mean(r['sec']) for r in rs]) / 60 if rs else None, 0)} |")
    below = [r for r in recs if r["es_sat"] is not None and r["rdd_best"] is not None and r["rdd_best"] < r["es_sat"]]
    L += ["", f"低于 ES 饱和值的 (实例, 规则)：{len(below)} 个"]
    for r in below:
        L.append(f"- {r['inst']} {r['rule']}：RDD 最好 {r['rdd_best']}（种子 {[v for v in r['rdd']]}，严格证明 {r['proven']}），ES 饱和值 {r['es_sat']}，默认规则 {r['base']}")
    return recs, L


def main():
    D = load()
    if len(sys.argv) > 1 and sys.argv[1] == "full":
        rho = float(sys.argv[2]); recs, L = full(D, rho)
        with open(os.path.join(ROOT, "results", "T8_rdd.jsonl.tmp"), "w") as f:
            for r in recs: f.write(json.dumps(r, ensure_ascii=False) + "\n")
        os.replace(os.path.join(ROOT, "results", "T8_rdd.jsonl.tmp"), os.path.join(ROOT, "results", "T8_rdd.jsonl"))
        print("\n".join(L))
    else:
        print("\n".join(sweep(D)))


if __name__ == "__main__":
    main()
