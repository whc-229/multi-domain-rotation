"""T12 汇总：大实例（m ∈ {8, 10, 15, 20}）上"松弛规则 + RDD-FD"，场景 = 全新启动、只在启动窗口内提前大修。
输入：results/T12/jobs/（base = fresh_rlp_w 严格窗口；rdd = RDD-FD 种子 0–2）、results/T6/rlp/ 与 results/T10/variants/（原 fresh_rlp_w，m20 的 3 个用 T10 后备）、
     results/T6/es/（ES lpgplan，m = 8、10，种子 0–2）、results/T6_bounds.jsonl（⌈LP_fr^Rw(4S)⌉）。
可行 = 严格证明或 60 万格复核 0 缺岗；只在可行解上比资源数。
输出 results/T12_large_rdd.jsonl（逐作业）与 markdown 表：
  1. 按 (m, 族)：平均资源、离下界 %、比 fresh_rlp_w（严格 / 原）少多少、比 ES 少多少（同实例配对）、窗口内提前大修次数、可行率与证明率、用时；
  2. 随时间的最好值（从每个作业的 traj 读 15 分钟、30 分钟、1 小时、2 小时、4 小时时的最好可行资源数）；
  3. 核对：最终解重新判定、抽样缺岗、窗口内提前大修。
用法：python3 code/t12_summary.py"""
import os, sys, json, glob, math
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names  # noqa: E402


def feas(d):
    """可行时返回资源数，否则 None。"""
    if not d or d.get("R") is None or d["R"] < 0: return None
    return d["R"] if ((d.get("verify_cert") == 1 and d.get("exact") == 1) or d.get("verify_gaps") == 0) else None


def jl(p):
    """读 json；不存在返回 None。"""
    return json.load(open(p)) if os.path.exists(p) else None


def mean(xs):
    """忽略 None 的平均。"""
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def fm(x, nd=1):
    """格式化。"""
    return "—" if x is None else f"{x:.{nd}f}"


def at_time(traj, sec):
    """traj = [[R, 代, 评价次数, 秒, 状态], …]：sec 秒时的最好可行资源数（第一项是默认规则或后备起点）。"""
    best = None
    for R, g, ev, t, st in traj:
        if t <= sec and R >= 0: best = R
    return best


def main():
    LB = {r["inst"]: math.ceil(r["LP_frRw"] - 1e-6) for r in map(json.loads, open(os.path.join(ROOT, "results", "T6_bounds.jsonl")))}
    rows = []
    for n in sorted(train_names([8, 10, 15, 20]), key=lambda n: (int(n[1:3]), n)):
        b = jl(os.path.join(ROOT, "results", "T12", "jobs", f"{n}__base.json"))
        o = jl(os.path.join(ROOT, "results", "T6", "rlp", f"{n}__fresh_rlp_w.json"))
        if o is not None: o["proven"] = int(o.get("verify_cert") == 1 and o.get("exact") == 1)
        if feas(o) is None and os.path.exists(os.path.join(ROOT, "results", "T10", "variants", f"{n}__hs1520.json")):
            o = jl(os.path.join(ROOT, "results", "T10", "variants", f"{n}__hs1520.json"))
        es = [feas(jl(os.path.join(ROOT, "results", "T6", "es", f"{n}__es_lpgplan__s{s}.json"))) for s in range(3)]
        for s in range(3):
            d = jl(os.path.join(ROOT, "results", "T12", "jobs", f"{n}__rdd__s{s}.json"))
            if d is None: continue
            rows.append(dict(inst=n, m=d["m"], family=d["family"], seed=s, LB=LB[n], R=feas(d), status=d["status"], proven=int(d.get("proven", 0)),
                             base_strict=feas(b), base_orig=feas(o), es=es, resets=d.get("resets"), last_reset=d.get("last_reset"), W=d["W"],
                             window_ok=d.get("window_ok"), f0_src=d.get("f0_src"), base_hsfb=d.get("base_hsfb"), wall=d["wall"], gens=d.get("gens"),
                             judges=d.get("judges"), judge_fail=d.get("judge_fail"), check_judge=d.get("check_judge"),
                             check_dump=(d.get("check_dump_same"), d.get("check_dump_n")), traj=d["traj"]))
    tmp = os.path.join(ROOT, "results", "T12_large_rdd.jsonl.tmp")
    with open(tmp, "w") as f:
        for r in rows: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(tmp, os.path.join(ROOT, "results", "T12_large_rdd.jsonl"))
    by = defaultdict(list)
    for r in rows: by[(r["m"], r["family"])].append(r)
    L = ["### 按 (m, 族)\n",
         "| m | 族 | 作业 | 可行 / 证明 | ⌈LP⌉ | RDD-FD 平均资源 | 离下界 % | fresh_rlp_w（严格）| 比它少（个，配对） | 少 / 同 / 多 | fresh_rlp_w（原）| 比它少 | ES lpgplan（种子平均）| 比 ES 少（配对） | 少 / 同 / 多 | 窗口内提前大修（次） | 平均用时 h |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for (m, fam), xs in sorted(by.items()):
        ok = [x for x in xs if x["R"] is not None]
        ds = [x["R"] - x["base_strict"] for x in ok if x["base_strict"] is not None]
        do = [x["R"] - x["base_orig"] for x in ok if x["base_orig"] is not None]
        de = [x["R"] - mean(x["es"]) for x in ok if mean(x["es"]) is not None]
        L.append(f"| {m} | {fam} | {len(xs)} | {len(ok)} / {sum(x['proven'] for x in ok)} | {fm(mean([x['LB'] for x in xs]))} | {fm(mean([x['R'] for x in ok]))} | "
                 f"{fm(mean([100 * (x['R'] - x['LB']) / x['LB'] for x in ok]))} | {fm(mean([x['base_strict'] for x in ok]))} | {fm(-mean(ds) if ds else None, 2)} | "
                 f"{sum(d < 0 for d in ds)} / {sum(d == 0 for d in ds)} / {sum(d > 0 for d in ds)} | {fm(mean([x['base_orig'] for x in ok]))} | {fm(-mean(do) if do else None, 2)} | "
                 f"{fm(mean([mean(x['es']) for x in ok]))} | {fm(-mean(de) if de else None, 2)} | "
                 f"{sum(d < 0 for d in de)} / {sum(d == 0 for d in de)} / {sum(d > 0 for d in de)} | {fm(mean([x['resets'] for x in ok]))} | {fm(mean([x['wall'] for x in xs]) / 3600, 2)} |")
    marks = [900, 1800, 3600, 7200, 14400]
    L += ["", "### 随时间的最好可行资源数（按 (m, 族) 平均；只算每个时刻都有可行解的作业）\n",
          "| m | 族 | 默认规则 | 15 分钟 | 30 分钟 | 1 小时 | 2 小时 | 4 小时 |", "|---|---|---|---|---|---|---|---|"]
    for (m, fam), xs in sorted(by.items()):
        vals = [[at_time(x["traj"], 0.0 + x["traj"][0][3])] + [at_time(x["traj"], t) for t in marks] for x in xs]
        cols = []
        for c in range(len(marks) + 1):
            if c > 0 and marks[c - 1] > (14400 if m >= 15 else 7200): cols.append("—"); continue
            v = [row[c] for row in vals]
            cols.append(fm(mean(v)) if all(y is not None for y in v) else fm(mean(v)) + "*")
        L.append(f"| {m} | {fam} | " + " | ".join(cols) + " |")
    cj = [r["check_judge"] for r in rows if r["check_judge"] is not None]
    L += ["", f"核对：最终解重新判定逐字段相同 {sum(cj)} / {len(cj)}；抽样缺岗与 eval 相同 {sum(r['check_dump'][0] or 0 for r in rows)} / {sum(r['check_dump'][1] or 0 for r in rows)}；"
          f"提前大修都在窗口内 {sum(r['window_ok'] == 1 for r in rows)} / {len(rows)}；默认规则用了后备补足时域的作业 {sum(r['base_hsfb'] == 1 for r in rows)}，用了后备起点（默认规则不可行）的 {sum((r['f0_src'] or 0) > 0 for r in rows)}。"]
    print("\n".join(L))


if __name__ == "__main__":
    main()
