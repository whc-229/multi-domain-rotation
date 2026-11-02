"""T15 汇总：下界松在哪里（m ∈ {2, 3} 的 24 个训练实例，不提前大修场景）。
来源：results/T15/jobs/（t15_job.py）；完整实例的下界与上界取 T9 的新下界与 T12c 之后的上界（results/T12c_ub.jsonl，reset = 0 的 T9_LB、T9_UB_new）。
(a) 只算机器人：下界 = ⌈LP_fr(4S)⌉（tools=False）；上界 = fresh_noreset_u 规则、RDD-FD（种子 0–2）、ES curdef（种子 0–2）里可行的最少资源数；
    与完整实例（机器人 + 工具）的差距对比，绝对差距与相对差距（差距 / 下界）都给。
(b) 只算机器人的完整整数模型（H = 2·S_max，HiGHS 单线程，2 小时）：记状态、整数解、对偶界；
    ⌈对偶界⌉ 与同一 H 的 ⌈LP(2S)⌉ 比，看整数性能否抬高下界；解到最优时与 ⌈LP(2S)⌉ 比。
可行 = 严格证明或 60 万格复核 0 缺岗。
输出 results/T15_bound_diag.jsonl（逐实例）与 markdown 表（打印）。
用法：python3 code/t15_summary.py"""
import os, sys, json
from collections import defaultdict, Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names  # noqa: E402

UBK = ["ub_rlp"] + [f"ub_fd__s{s}" for s in range(3)] + [f"ub_es__s{s}" for s in range(3)]


def jl(n, k):
    """读 results/T15/jobs/<实例>__<种类>.json；不存在返回 None。"""
    p = os.path.join(ROOT, "results", "T15", "jobs", f"{n}__{k}.json")
    return json.load(open(p)) if os.path.exists(p) else None


def feas(d):
    """可行时返回资源数，否则 None。"""
    if not d or d.get("R") is None or d["R"] < 0: return None
    return d["R"] if ((d.get("verify_cert") == 1 and d.get("exact") == 1) or d.get("verify_gaps") == 0) else None


def best_of(n, prefix):
    """某一类上界（ub_rlp / ub_fd / ub_es）在各种子里的最少资源数。"""
    v = [feas(jl(n, k)) for k in UBK if k.startswith(prefix)]
    v = [x for x in v if x is not None]
    return min(v) if v else None


def main():
    full = {r["inst"]: r for r in map(json.loads, open(os.path.join(ROOT, "results", "T12c_ub.jsonl"))) if r["reset"] == 0}
    rows = []
    for n in sorted(train_names([2, 3]), key=lambda n: (int(n[1:3]), n)):
        lp, mip = jl(n, "lp"), jl(n, "mip")
        ubs = {k: feas(jl(n, k)) for k in UBK}
        ok = {k: v for k, v in ubs.items() if v is not None}
        UB = min(ok.values()) if ok else None
        src = sorted({k.split("__")[0] for k, v in ok.items() if v == UB})
        f = full[n]
        b = mip.get("bound") if mip else None
        rows.append(dict(inst=n, m=int(n[1:3]), family=n.split("-")[1],
                         r_LP4=lp["LP4"], r_LP2=lp["LP2"], r_LB=lp["LB"], r_LB2=lp["LB2"], r_N4=lp["N4"],
                         r_UB=UB, r_UB_src=src, r_UB_rlp=best_of(n, "ub_rlp"), r_UB_fd=best_of(n, "ub_fd"), r_UB_es=best_of(n, "ub_es"),
                         r_tools_in_sol=max((jl(n, k) or {}).get("tools_in_sol", 0) for k in UBK),
                         r_gap=None if UB is None else UB - lp["LB"],
                         f_LB=f["T9_LB"], f_UB=f["T9_UB_new"], f_gap=f["T9_UB_new"] - f["T9_LB"],
                         mip_status=mip and mip["status"], mip_obj=mip and mip["obj"], mip_bound=b, mip_LB=mip and mip["LB_int"], mip_sec=mip and mip["sec"],
                         mip_raises=int(mip is not None and mip["LB_int"] is not None and mip["LB_int"] > lp["LB2"]),
                         mip_opt_above=int(mip is not None and mip["status"] == "Optimal" and mip["obj"] is not None and round(mip["obj"]) > lp["LB2"])))
    tmp = os.path.join(ROOT, "results", "T15_bound_diag.jsonl.tmp")
    with open(tmp, "w") as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(tmp, tmp[:-4])

    by = defaultdict(list)
    for r in rows: by[(r["m"], r["family"])].append(r)
    L = ["### (a) 只算机器人与完整实例的差距（按 (m, 族) 平均；相对差距 = 差距 / 下界）\n",
         "| m | 族 | 实例 | 机器人：下界 | 机器人：上界 | 机器人：差距 | 机器人：相对差距 | 完整：下界 | 完整：上界 | 完整：差距 | 完整：相对差距 | 机器人差距 < 完整差距 |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for (m, fam), xs in sorted(by.items()):
        L.append(f"| {m} | {fam} | {len(xs)} | {np.mean([x['r_LB'] for x in xs]):.2f} | {np.mean([x['r_UB'] for x in xs]):.2f} | {np.mean([x['r_gap'] for x in xs]):.2f} | "
                 f"{100 * np.mean([x['r_gap'] / x['r_LB'] for x in xs]):.1f}% | {np.mean([x['f_LB'] for x in xs]):.2f} | {np.mean([x['f_UB'] for x in xs]):.2f} | "
                 f"{np.mean([x['f_gap'] for x in xs]):.2f} | {100 * np.mean([x['f_gap'] / x['f_LB'] for x in xs]):.1f}% | {sum(x['r_gap'] < x['f_gap'] for x in xs)} / {len(xs)} |")
    L.append(f"\n合计：机器人差距平均 {np.mean([x['r_gap'] for x in rows]):.2f}（相对 {100 * np.mean([x['r_gap'] / x['r_LB'] for x in rows]):.1f}%），"
             f"完整差距平均 {np.mean([x['f_gap'] for x in rows]):.2f}（相对 {100 * np.mean([x['f_gap'] / x['f_LB'] for x in rows]):.1f}%）；"
             f"机器人差距为 0 的实例 {sum(x['r_gap'] == 0 for x in rows)} / {len(rows)}；机器人差距的分布 {dict(sorted(Counter(x['r_gap'] for x in rows).items()))}。")
    L += ["", "### (a) 只算机器人：各上界方法（各种子最好；按 (m, 族) 平均）与最好值来源\n",
          "| m | 族 | 规则 fresh_noreset_u | RDD-FD | ES curdef | 最好 | 最好值来自（实例数） |", "|---|---|---|---|---|---|---|"]
    for (m, fam), xs in sorted(by.items()):
        c = Counter(s for x in xs for s in x["r_UB_src"])
        L.append(f"| {m} | {fam} | {np.mean([x['r_UB_rlp'] for x in xs]):.2f} | {np.mean([x['r_UB_fd'] for x in xs]):.2f} | {np.mean([x['r_UB_es'] for x in xs]):.2f} | "
                 f"{np.mean([x['r_UB'] for x in xs]):.2f} | " + "，".join(f"{k} {v}" for k, v in sorted(c.items())) + " |")
    L += ["", "### 逐实例\n", "| 实例 | 机器人 LP(4S) | 机器人下界 | 机器人上界 | 机器人差距 | 完整下界 | 完整上界 | 完整差距 | MIP 状态 | MIP 整数解 | ⌈MIP 对偶界⌉ | ⌈LP(2S)⌉ |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for x in rows:
        L.append(f"| {x['inst']} | {x['r_LP4']:.2f} | {x['r_LB']} | {x['r_UB']} | {x['r_gap']} | {x['f_LB']} | {x['f_UB']} | {x['f_gap']} | "
                 f"{'最优' if x['mip_status'] == 'Optimal' else '限时'} | {'—' if x['mip_obj'] is None else int(x['mip_obj'])} | {x['mip_LB']} | {x['r_LB2']} |")
    st = Counter(x["mip_status"] for x in rows)
    L += ["", f"### (b) 只算机器人的完整整数模型（H = 2·S_max）\n",
          f"状态：{dict(st)}；⌈对偶界⌉ > ⌈LP(2S)⌉ 的实例 {sum(x['mip_raises'] for x in rows)} / {len(rows)}；"
          f"解到最优且最优值 > ⌈LP(2S)⌉ 的实例 {sum(x['mip_opt_above'] for x in rows)}；"
          f"⌈LP(2S)⌉ = ⌈LP(4S)⌉ 的实例 {sum(x['r_LB2'] == x['r_LB'] for x in rows)} / {len(rows)}；"
          f"限时内的整数解比机器人上界差的实例 {sum(x['mip_obj'] is not None and x['mip_obj'] > x['r_UB'] for x in rows)}，没有整数解 {sum(x['mip_obj'] is None for x in rows)}。"]
    print("\n".join(L))


if __name__ == "__main__":
    main()
