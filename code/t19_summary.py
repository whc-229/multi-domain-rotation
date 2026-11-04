"""T19 汇总：接入构造（(a) 稳态规则 σ 作为默认规则；(c) 最好 RDD 解的循环作 σ 的诊断）。
来源：results/T19/steady/、results/T19/best/（code/t19_connector.py）。
对照：
- fresh_rlp_w（严格窗口）的规则本身：m = 4、6 取 T12b 的 base，m = 8、10 取 T12 的 base（与 T18 的 W = 0 相同）；
- 已知最好 RDD 解（T8 / T11 / T16 β = 0 / T18，与 (c) 选 σ 的口径相同）；
- 下界 ⌈LP_fr^Rw⌉（T6，允许启动期提前大修）。
每个 (实例, 来源) 的"接入构造"取 8 个时间原点 × 2 种入场口径里可行（严格证明或 60 万格 0 缺岗）且资源最少的；另给 first 口径各 θ 的中位数。
输出 results/T19_connector.jsonl（逐实例、逐来源）与 markdown 表（打印）。用法：python3 code/t19_summary.py"""
import os, sys, json, glob, math
from collections import defaultdict, Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402
from t18_summary import job as t18_job, feas  # noqa: E402
from t19_connector import best_record  # noqa: E402


def variants(d):
    """一条记录里所有 (θ, 变体)。"""
    return [(t, v) for t in d["thetas"] for v in t.get("variants", [])]


def summarize(d, lb):
    """一条 T19 记录 → 汇总行。"""
    n = d["inst"]; S = d["S"]
    vs = variants(d)
    ok = [(t, v) for t, v in vs if v["feasible"]]
    first_ok = [v["R"] for t, v in ok if v["mode"] == "first"]
    best = {}
    for mode in ("first", "last"):
        xs = [v for t, v in ok if v["mode"] == mode]
        best[mode] = min(v["R"] for v in xs) if xs else None
    b = d["best"]
    row = dict(inst=n, m=d["m"], family=d["family"], kind=d["kind"], S=S, LB=lb.get(n), src=d["src"], R_sigma=d["R_sigma"], N_sigma=d["N_sigma"], M_sigma=d["M_sigma"],
               t1=d["t1"], t1_from_cert=d["t1_from_cert"], repro_same=d["src_repro"]["same"],
               n_theta_ok=sum(1 for t in d["thetas"] if t["y_gaps"] == 0), n_var=len(vs), n_feas=len(ok),
               n_cont_fail=sum(1 for t, v in vs if v["plan_gaps"] == 0 and not v["feasible"]), n_plan_gap=sum(1 for t, v in vs if v["plan_gaps"] > 0),
               R_first=best["first"], R_last=best["last"], R_first_med=float(np.median(first_ok)) if first_ok else None,
               C=[t.get("C") for t in d["thetas"] if t.get("C") is not None],
               conn=None)
    if b is not None:
        bv = [v for t, v in vs if t["theta"] == b["theta"] and v["mode"] == b["mode"]][0]
        row["conn"] = dict(R=b["R"], dN=b["dN"], dN_tot=sum(b["dN"]), theta=b["theta"], mode=b["mode"], proven=b["proven"],
                           T0=b["T0"], Tplan=b["Tplan"], T0_over_S=round(b["T0"] / S, 3), Tplan_over_S=round(b["Tplan"] / S, 3),
                           planW=bv["planW"], widened=int(bv["planW"] > 2 * S), n_gap=bv["n_gap"], n_pre=bv["n_pre"], n_ext=bv["n_ext"], n_reset=bv["n_reset"])
    return row


def main():
    lb = {r["inst"]: math.ceil(r["LP_frRw"] - 1e-6) for r in map(json.loads, open(os.path.join(ROOT, "results", "T6_bounds.jsonl")))}
    rows = []
    for kind in ("steady", "best"):
        for f in sorted(glob.glob(os.path.join(ROOT, "results", "T19", kind, "*.json")), key=lambda f: (int(os.path.basename(f)[1:3]), f)):
            d = json.load(open(f)); r = summarize(d, lb)
            if kind == "steady":
                base = t18_job(d["inst"], 0.0, None)
                r["fresh_rlp_w"] = dict(R=base.get("R"), feasible=int(feas(base) is not None), proven=int(base.get("verify_cert") == 1 and base.get("exact") == 1),
                                        src="T12b" if d["m"] <= 6 else "T12")
                bb = best_record(d["inst"]); r["known_best"] = dict(R=bb[1]["R"], job=os.path.relpath(bb[0], ROOT)) if bb else None
            rows.append(r)
    tmp = os.path.join(ROOT, "results", "T19_connector.jsonl.tmp")
    with open(tmp, "w") as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(tmp, tmp[:-4])

    A = [r for r in rows if r["kind"] == "steady"]; Cc = [r for r in rows if r["kind"] == "best"]
    L = ["### (a) 接入构造（σ = 稳态松弛规则）作为默认规则 vs fresh_rlp_w（严格窗口）：平均资源数（按 (m, 族)）\n",
         "| m | 族 | 实例 | 下界 | σ 本身（稳态初值） | 接入构造 | ΔN | fresh_rlp_w | 构造少 / 同 / 多（fresh_rlp_w 可行的实例） | 已知最好 RDD 解 | 接入完成 T_plan / S_max |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    by = defaultdict(list)
    for r in A: by[(r["m"], r["family"])].append(r)
    for (m, fam), xs in sorted(by.items()):
        ok = [x for x in xs if x["conn"]]
        fb = [x for x in ok if x["fresh_rlp_w"]["feasible"]]
        w = sum(x["conn"]["R"] < x["fresh_rlp_w"]["R"] for x in fb); t = sum(x["conn"]["R"] == x["fresh_rlp_w"]["R"] for x in fb); l = len(fb) - w - t
        L.append(f"| {m} | {fam} | {len(ok)}/{len(xs)} | {np.mean([x['LB'] for x in ok]):.1f} | {np.mean([x['R_sigma'] for x in ok]):.1f} | {np.mean([x['conn']['R'] for x in ok]):.1f} | "
                 f"{np.mean([x['conn']['dN_tot'] for x in ok]):.1f} | {np.mean([x['fresh_rlp_w']['R'] for x in fb]):.1f}（{len(fb)} 个可行） | {w} / {t} / {l} | "
                 f"{np.mean([x['known_best']['R'] for x in ok]):.1f} | {np.median([x['conn']['Tplan_over_S'] for x in ok]):.2f}（{min(x['conn']['Tplan_over_S'] for x in ok):.2f}–{max(x['conn']['Tplan_over_S'] for x in ok):.2f}） |")
    okA = [x for x in A if x["conn"]]; fbA = [x for x in okA if x["fresh_rlp_w"]["feasible"]]
    dif = [x["conn"]["R"] - x["fresh_rlp_w"]["R"] for x in fbA]
    L.append(f"\n合计：接入构造可行 {len(okA)} / {len(A)}；与 fresh_rlp_w 比（fresh_rlp_w 可行的 {len(fbA)} 个）：少 {sum(d < 0 for d in dif)}、同 {sum(d == 0 for d in dif)}、多 {sum(d > 0 for d in dif)}，"
             f"平均差 {np.mean(dif):+.2f}，差的分布 {dict(sorted(Counter(dif).items()))}；fresh_rlp_w 不可行的实例 {[x['inst'] for x in okA if not x['fresh_rlp_w']['feasible']]}（接入构造在这些实例上可行）。")
    L.append(f"ΔN（额外机器人）分布 {dict(sorted(Counter(x['conn']['dN_tot'] for x in okA).items()))}；ΔM = 0（构造保证每格每站的 (型, 时刻) 与参照排班相同）。"
             f"严格证明 {sum(x['conn']['proven'] for x in okA)} / {len(okA)}，其余 60 万格 0 缺岗。")
    L.append(f"比已知最好 RDD 解多：平均 {np.mean([x['conn']['R'] - x['known_best']['R'] for x in okA]):+.2f}，分布 {dict(sorted(Counter(x['conn']['R'] - x['known_best']['R'] for x in okA).items()))}；"
             f"σ 本身（稳态初值）比已知最好多 {np.mean([x['R_sigma'] - x['known_best']['R'] for x in okA]):+.2f}。")
    L += ["", "### 入场口径与 θ（(a)、(c) 合并按来源列）\n",
          "| 来源 | 实例 | 可用 θ（Y_θ 无缺岗，均值 / 8） | 变体可行 / 总数 | 计划期 0 缺岗但规则接手后缺岗 | first 最好（均值） | last 最好（均值） | first 各 θ 中位（均值） | T₀ = c₀ + C_σ（/ S_max，中位） | T_plan > 2·S_max（提前大修窗口放宽） |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for kind, xs in (("steady", A), ("best", Cc)):
        both = [x for x in xs if x["R_first"] is not None and x["R_last"] is not None]
        T0s = [x["conn"]["T0_over_S"] for x in xs if x["conn"]]
        L.append(f"| {kind} | {len(xs)} | {np.mean([x['n_theta_ok'] for x in xs]):.1f} | {sum(x['n_feas'] for x in xs)} / {sum(x['n_var'] for x in xs)} | {sum(x['n_cont_fail'] for x in xs)} | "
                 f"{np.mean([x['R_first'] for x in both]):.1f} | {np.mean([x['R_last'] for x in both]):.1f}（{len(both)} 个两者都可行） | "
                 f"{np.mean([x['R_first_med'] for x in xs if x['R_first_med'] is not None]):.1f} | {np.median(T0s):.2f} | "
                 f"{', '.join(x['inst'] for x in xs if x['conn'] and x['conn']['widened']) or '无'} |")
    L += ["", "### (c) 诊断：σ = 已知最好 RDD 解的循环，接入上界 N_σ + ΔN（工具 M_σ 不变）vs 该解本身\n",
          "| m | 族 | 实例 | 有接入上界 | 最好解本身（平均） | 接入上界（平均） | 相等 | ΔN 分布 | 接入完成 T_plan / S_max（中位） |", "|---|---|---|---|---|---|---|---|---|"]
    byc = defaultdict(list)
    for r in Cc: byc[(r["m"], r["family"])].append(r)
    for (m, fam), xs in sorted(byc.items()):
        ok = [x for x in xs if x["conn"]]
        L.append(f"| {m} | {fam} | {len(xs)} | {len(ok)} | {np.mean([x['R_sigma'] for x in ok]):.1f} | {np.mean([x['conn']['R'] for x in ok]):.1f} | {sum(x['conn']['dN_tot'] == 0 for x in ok)} / {len(ok)} | "
                 f"{dict(sorted(Counter(x['conn']['dN_tot'] for x in ok).items()))} | {np.median([x['conn']['Tplan_over_S'] for x in ok]):.2f} |")
    okC = [x for x in Cc if x["conn"]]
    L.append(f"\n合计：有接入上界 {len(okC)} / {len(Cc)}（其余实例 8 个 θ 的 Y_θ 都有缺岗：{[x['inst'] for x in Cc if not x['conn']]}）；"
             f"相等 {sum(x['conn']['dN_tot'] == 0 for x in okC)} / {len(okC)}，ΔN 分布 {dict(sorted(Counter(x['conn']['dN_tot'] for x in okC).items()))}，平均 {np.mean([x['conn']['dN_tot'] for x in okC]):.2f}；"
             f"σ 来源 {dict(Counter(x['src']['job'].split('/')[1] for x in Cc))}。")
    L += ["", "### 逐实例（(a)）\n", "| 实例 | 下界 | σ | 接入构造（ΔN；θ、口径） | fresh_rlp_w | 已知最好 | T_plan / S | 证明 |", "|---|---|---|---|---|---|---|---|"]
    for x in A:
        c = x["conn"]; fr = x["fresh_rlp_w"]
        L.append(f"| {x['inst']} | {x['LB']} | {x['R_sigma']} | " + (f"{c['R']}（{c['dN_tot']}；θ{c['theta']}、{c['mode']}）" if c else "—") + " | "
                 f"{fr['R']}{'' if fr['feasible'] else '（不可行）'} | {x['known_best']['R']} | " + (f"{c['Tplan_over_S']:.2f} | {'是' if c['proven'] else '复核'}" if c else "— | —") + " |")
    L += ["", "### 逐实例（(c)）\n", "| 实例 | 最好解（来源） | 接入上界（ΔN；θ、口径） | 可用 θ | T_plan / S |", "|---|---|---|---|---|"]
    for x in Cc:
        c = x["conn"]
        L.append(f"| {x['inst']} | {x['R_sigma']}（{x['src']['job'].split('/')[1]} {os.path.basename(x['src']['job'])[:-5].split('__', 1)[1]}） | "
                 + (f"{c['R']}（{c['dN_tot']}；θ{c['theta']}、{c['mode']}）" if c else "—") + f" | {x['n_theta_ok']} | " + (f"{c['Tplan_over_S']:.2f}" if c else "—") + " |")
    L.append(f"\n核对：σ 来源的 Python 复现与 C++ 缺岗、提前大修次数相同 {sum(x['repro_same'] for x in rows)} / {len(rows)}；"
             f"计划期有缺岗的变体 {sum(x['n_plan_gap'] for x in rows)} 个；σ 进入循环的时刻由 C++ 判环给出 {sum(x['t1_from_cert'] for x in rows)} / {len(rows)}（其余取 8·S_max）。")
    print("\n".join(L))


if __name__ == "__main__":
    main()
