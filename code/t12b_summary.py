"""T12b 汇总：严格窗口（提前大修只在 [0, 2·S_max) 内执行）下重算三组数字，旧值 / 新值对照（写进 T12 报告）。
1. fresh_rlp_w：旧 = T6 / T6b 的结果（results/T6/rlp/，m = 20 定不出的 3 个用 T10 的后备补足时域结果）；
   新 = 严格窗口（m ≤ 6：results/T12b/jobs/；m ≥ 8：results/T12/jobs/<实例>__base.json）。
2. 规则与 ES（m ≤ 10，种子 0–2 的平均，同实例配对）：
   旧 = 原 fresh_rlp_w 对 T6 的 ES lpgplan；新 = 严格窗口 fresh_rlp_w 对严格窗口的 ES lpgplan（results/T12b/es/）。
   旧 ES 解是否在窗口外提前大修：用 rost5 eval 重放 3000 格（同一组选项）看最后一次提前大修时刻。
3. T9 的"启动期提前大修"场景：旧上界 = results/T9_intfleet.jsonl 的 UB（不提前大修的最好上界与原 fresh_rlp_w 的较小者）；
   新上界 = 不提前大修的最好上界（code/t9_summary.ub_noreset）与严格窗口 fresh_rlp_w、严格窗口 ES 的最小者；下界不变（T9 的新下界）。
可行 = 严格证明或 60 万格复核 0 缺岗。输出 results/T12b_strict.jsonl（逐实例）与 markdown 表（打印）。
用法：python3 code/t12b_summary.py"""
import os, sys, json, subprocess
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import t9_summary  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")


def jl(p):
    """读 json；不存在返回 None。"""
    return json.load(open(p)) if os.path.exists(p) else None


def feas(d):
    """可行时返回资源数，否则 None。"""
    if not d or d.get("R") is None or d["R"] < 0: return None
    return d["R"] if ((d.get("verify_cert") == 1 and d.get("exact") == 1) or d.get("verify_gaps") == 0) else None


def mean(xs):
    """忽略 None 的平均；全是 None 返回 None。"""
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def fm(x, nd=1):
    """格式化。"""
    return "—" if x is None else f"{x:.{nd}f}"


def old_es_last_reset(d, txt):
    """重放旧 ES 解 3000 格，返回最后一次提前大修时刻。"""
    opts = [o for o in d["cmd"].split() if o.startswith(("lpg=", "planf="))]
    a = [R5, "eval", txt, "3000", "0"] + [str(x) for x in d["Nk"]] + ["|"] + [str(x) for x in d["Mj"]] + ["|"] + [str(x) for x in d["w"]] + ["|"] + opts + ["curdef"]
    return json.loads(subprocess.run(a, capture_output=True, text=True, check=True).stdout)["last_reset"]


def main():
    T4b = {(r["inst"], r["reset"]): r for r in map(json.loads, open(os.path.join(ROOT, "results", "T4b_gap.jsonl")))}
    T9 = {r["inst"]: r for r in map(json.loads, open(os.path.join(ROOT, "results", "T9_intfleet.jsonl"))) if r["reset"] == 1}
    rows = []
    for n in sorted(train_names([2, 3, 4, 6, 8, 10, 15, 20]), key=lambda n: (int(n[1:3]), n)):
        I = load_train(n); S = S_max(I); W = 2 * S; m = I.meta["m"]; _, txt = inst_paths(n, "train")
        o = jl(os.path.join(ROOT, "results", "T6", "rlp", f"{n}__fresh_rlp_w.json"))
        if feas(o) is None and os.path.exists(os.path.join(ROOT, "results", "T10", "variants", f"{n}__hs1520.json")):
            o = jl(os.path.join(ROOT, "results", "T10", "variants", f"{n}__hs1520.json"))
        s = jl(os.path.join(ROOT, "results", "T12b" if m <= 6 else "T12", "jobs", f"{n}__base.json"))
        row = dict(inst=n, m=m, family=I.meta["family"], W=W, rule_old=feas(o), rule_new=feas(s),
                   rule_new_proven=int(bool(s) and s.get("verify_cert") == 1 and s.get("exact") == 1),
                   rule_new_window_ok=int(bool(s) and s.get("last_reset", -1) < W), rule_new_hsfb=(s or {}).get("hsfb_used", (s or {}).get("base_hsfb")))
        if m <= 10:
            eo, en, lr, wn = [], [], [], []
            for sd in range(3):
                d = jl(os.path.join(ROOT, "results", "T6", "es", f"{n}__es_lpgplan__s{sd}.json"))
                eo.append(feas(d)); lr.append(old_es_last_reset(d, txt) if feas(d) is not None else None)
                e = jl(os.path.join(ROOT, "results", "T12b", "es", f"{n}__es_lpgplan_w__s{sd}.json"))
                en.append(feas(e) if e and e.get("window_ok") == 1 else None); wn.append(e.get("window_ok") if e else None)
            row.update(es_old=eo, es_old_late=sum(1 for x in lr if x is not None and x >= W), es_new=en, es_new_window_ok=wn)
        if n in T9:
            t = T9[n]; ub = t9_summary.ub_noreset(n, T4b)
            c = [x for x in [ub, row["rule_new"]] + (row.get("es_new") or []) if x is not None]
            row.update(T9_LB=t["new_LB"], T9_UB_old=t["UB"], T9_UB_new=min(c) if c else None)
        rows.append(row)
    tmp = os.path.join(ROOT, "results", "T12b_strict.jsonl.tmp")
    with open(tmp, "w") as f:
        for r in rows: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(tmp, tmp[:-4])

    by = defaultdict(list)
    for r in rows: by[(r["m"], r["family"])].append(r)
    L = ["### 1. fresh_rlp_w：原结果 / 严格窗口\n",
         "| m | 族 | 实例 | 可行（旧 / 新） | 平均资源（旧） | 平均资源（新，配对） | 新 − 旧 | 变多 / 不变 / 变少 | 新：严格证明 | 新：提前大修都在窗口内 |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for (m, fam), xs in sorted(by.items()):
        pr = [x for x in xs if x["rule_old"] is not None and x["rule_new"] is not None]
        dd = [x["rule_new"] - x["rule_old"] for x in pr]
        L.append(f"| {m} | {fam} | {len(xs)} | {sum(x['rule_old'] is not None for x in xs)} / {sum(x['rule_new'] is not None for x in xs)} | "
                 f"{fm(mean([x['rule_old'] for x in pr]))} | {fm(mean([x['rule_new'] for x in pr]))} | {fm(mean(dd), 2)} | "
                 f"{sum(d > 0 for d in dd)} / {sum(d == 0 for d in dd)} / {sum(d < 0 for d in dd)} | {sum(x['rule_new_proven'] for x in xs)} | {sum(x['rule_new_window_ok'] for x in xs)} |")
    L += ["", "### 2. 规则与 ES（ES 取 3 个种子的平均，同实例配对；规则 − ES < 0 表示规则更少）\n",
          "| m | 族 | 旧 ES 解在窗口外提前大修（个 / 可行解） | 旧：规则 − ES | 旧：规则少 / 同 / 多 | 新 ES 可行且在窗口内 | 新 ES 平均 | 新：规则 − ES | 新：规则少 / 同 / 多 |",
          "|---|---|---|---|---|---|---|---|---|"]
    for (m, fam), xs in sorted(by.items()):
        if m > 10: continue
        do = [x["rule_old"] - mean(x["es_old"]) for x in xs if x["rule_old"] is not None and mean(x["es_old"]) is not None]
        dn = [x["rule_new"] - mean(x["es_new"]) for x in xs if x["rule_new"] is not None and mean(x["es_new"]) is not None]
        L.append(f"| {m} | {fam} | {sum(x['es_old_late'] for x in xs)} / {sum(sum(e is not None for e in x['es_old']) for x in xs)} | {fm(mean(do), 2)} | "
                 f"{sum(d < 0 for d in do)} / {sum(d == 0 for d in do)} / {sum(d > 0 for d in do)} | {sum(sum(e is not None for e in x['es_new']) for x in xs)} / {3 * len(xs)} | "
                 f"{fm(mean([mean(x['es_new']) for x in xs]))} | {fm(mean(dn), 2)} | {sum(d < 0 for d in dn)} / {sum(d == 0 for d in dn)} / {sum(d > 0 for d in dn)} |")
    L += ["", "### 3. T9\"启动期提前大修\"场景的上界与差距（下界 = T9 的新下界，不变）\n",
          "| m | 族 | 实例 | 下界 | 上界（旧） | 上界（新） | 差距（旧） | 差距（新） | 上界变化的实例 |", "|---|---|---|---|---|---|---|---|---|"]
    for (m, fam), xs in sorted(by.items()):
        t = [x for x in xs if x.get("T9_LB") is not None and x.get("T9_UB_old") is not None and x.get("T9_UB_new") is not None]
        if not t: continue
        L.append(f"| {m} | {fam} | {len(t)} | {fm(mean([x['T9_LB'] for x in t]), 2)} | {fm(mean([x['T9_UB_old'] for x in t]), 2)} | {fm(mean([x['T9_UB_new'] for x in t]), 2)} | "
                 f"{fm(mean([x['T9_UB_old'] - x['T9_LB'] for x in t]), 2)} | {fm(mean([x['T9_UB_new'] - x['T9_LB'] for x in t]), 2)} | {sum(x['T9_UB_new'] != x['T9_UB_old'] for x in t)} |")
    print("\n".join(L))


if __name__ == "__main__":
    main()
