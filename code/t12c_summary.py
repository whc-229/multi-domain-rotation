"""T12c 汇总：作为上界的 ES 解能否用全精度权重复现；更新 T4、T4b、T9 的上界与差距（旧值 / 新值对照），重新确认 T8 的 ES 饱和值。
输入：results/T12c/jobs/（逐实例，t12c_job.py）、results/T12c/es/（逐个补跑的 ES 解）、results/T4_exact.jsonl、results/T4b_gap.jsonl、
     results/T9_intfleet.jsonl、results/T12b_strict.jsonl（T9 启动期提前大修场景的严格窗口上界）、T8 / T6 的结果（t9_summary.ub_noreset 的其他来源）。
- T4：上界口径 = T3 + T4 的解；T4b：再加 es100。下界不变。
- T9（m = 2、3 才受影响）：不提前大修上界 = min(T4b 新上界, T8 的 RDD 与默认规则, T6 fresh_noreset_u)；
  启动期提前大修上界 = min(上一行, 严格窗口 fresh_rlp_w, 严格窗口 ES)（T12b 的口径）。下界 = T9 的新下界。
- T8 的 ES 饱和值 = T4b 不提前大修的上界。
输出 results/T12c_ub.jsonl（逐实例、逐场景）和 markdown 表（打印）。不改任何旧结果文件。
用法：python3 code/t12c_summary.py"""
import os, sys, json, glob
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402
import t9_summary  # noqa: E402

T8_BELOW = [("m02-K1-01", "es", 28), ("m02-K12-00", "es", 47), ("m02-K12-05", "rlp", 49), ("m03-K1-01", "es", 42),
            ("m03-K1-01", "hand", 42), ("m03-K1-04", "rlp", 43), ("m03-K1-04", "es", 43)]   # T8 报告里"比 ES 饱和值少 1 个"的 (实例, 规则, RDD 最好)


def mean(xs):
    """忽略 None 的平均。"""
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def fm(x, nd=2):
    """格式化。"""
    return "—" if x is None else f"{x:.{nd}f}"


def main():
    T4 = {(r["inst"], r["reset"]): r for r in map(json.loads, open(os.path.join(ROOT, "results", "T4_exact.jsonl")))}
    T4b = {(r["inst"], r["reset"]): r for r in map(json.loads, open(os.path.join(ROOT, "results", "T4b_gap.jsonl")))}
    T9 = {(r["inst"], r["reset"]): r for r in map(json.loads, open(os.path.join(ROOT, "results", "T9_intfleet.jsonl")))}
    S12 = {r["inst"]: r for r in map(json.loads, open(os.path.join(ROOT, "results", "T12b_strict.jsonl")))}
    J = {json.load(open(f))["inst"]: json.load(open(f)) for f in glob.glob(os.path.join(ROOT, "results", "T12c", "jobs", "*.json"))}
    rows, bad = [], []
    T4b_new = {}
    for n in sorted(J, key=lambda n: (int(n[1:3]), n)):
        j = J[n]; m = int(n[1:3]); fam = n.split("-")[1]
        for rs in (0, 1):
            a, b = j[f"r{rs}_T4"], j[f"r{rs}_T4b"]
            # 核对：走候选解时的旧上界与 T4 / T4b 汇总的上界相同
            if a["UB_old"] != T4[(n, rs)]["UB"]: bad.append((n, rs, "T4", a["UB_old"], T4[(n, rs)]["UB"]))
            if b["UB_old"] != T4b[(n, rs)]["UB_new"]: bad.append((n, rs, "T4b", b["UB_old"], T4b[(n, rs)]["UB_new"]))
            T4b_new[(n, rs)] = b["UB_new"]
            rows.append(dict(inst=n, m=m, family=fam, reset=rs, LB_T4=T4[(n, rs)]["LB"], LB_T4b=T4b[(n, rs)]["LB"],
                             T4_UB_old=a["UB_old"], T4_UB_new=a["UB_new"], T4_src=a["src"], T4b_UB_old=b["UB_old"], T4b_UB_new=b["UB_new"], T4b_src=b["src"],
                             tried=len({t["file"] for t in a["tried"] + b["tried"]})))
    # T9：m = 2、3 用新的 T4b 上界重算
    T4b_patched = {k: dict(v, UB_new=T4b_new.get(k, v["UB_new"])) for k, v in T4b.items()}
    for r in rows:
        n, rs = r["inst"], r["reset"]
        ub0_old = t9_summary.ub_noreset(n, T4b); ub0_new = t9_summary.ub_noreset(n, T4b_patched)
        if rs == 0:
            old, new = T9[(n, 0)]["UB"], ub0_new
            if old != ub0_old: bad.append((n, 0, "T9", old, ub0_old))
        else:
            s = S12[n]; extra = [x for x in [s["rule_new"]] + (s.get("es_new") or []) if x is not None]
            old = s["T9_UB_new"]; new = min([ub0_new] + extra)
            if old != min([ub0_old] + extra): bad.append((n, 1, "T9(T12b)", old, min([ub0_old] + extra)))
        r.update(T9_LB=T9[(n, rs)]["new_LB"], T9_UB_old=old, T9_UB_new=new)
    tmp = os.path.join(ROOT, "results", "T12c_ub.jsonl.tmp")
    with open(tmp, "w") as f:
        for r in rows: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(tmp, tmp[:-4])

    es = [json.load(open(f)) for f in glob.glob(os.path.join(ROOT, "results", "T12c", "es", "*.json"))]
    L = [f"### 补跑的 ES 解\n",
         f"补跑 {len(es)} 个；与原作业逐字段相同 {sum(e['same_as_orig'] for e in es)}；用 w17 重新判定与原判定相同 {sum(e['rejudge_w17'] == 1 for e in es)} / {sum(e['rejudge_w17'] is not None for e in es)}；"
         f"用 4 位小数的 w 重新判定相同 {sum(e['rejudge_w4'] for e in es)} / {len(es)}；能复现 {sum(e['reproducible'] for e in es)}。"]
    nr = [e for e in es if not e["reproducible"] or e["rejudge_w4"] != 1 or not e["same_as_orig"]]
    if nr:
        L += ["", "| 原作业 | R | 与原作业相同 | w17 复现 | 4 位小数复现 | 原（缺岗, 判环, 提前大修） | w17 | 4 位小数 |", "|---|---|---|---|---|---|---|---|"]
        for e in sorted(nr, key=lambda e: e["file"]):
            o, a, b = e["orig"], e["v17"] or {}, e["v4"]
            L.append(f"| {os.path.basename(e['file'])[:-5]} | {e['R']} | {e['same_as_orig']} | {e['rejudge_w17']} | {e['rejudge_w4']} | "
                     f"({o['verify_gaps']}, {o['verify_cert']}, {o['resets']}) | ({a.get('gaps')}, {a.get('cert')}, {a.get('resets')}) | ({b['gaps']}, {b['cert']}, {b['resets']}) |")
    by = defaultdict(list)
    for r in rows: by[(r["m"], r["family"], r["reset"])].append(r)
    for tag, lb, uo, un in (("T4（T3 + T4 的解）", "LB_T4", "T4_UB_old", "T4_UB_new"), ("T4b（再加 es100）", "LB_T4b", "T4b_UB_old", "T4b_UB_new"),
                            ("T9（下界 = T9 的新下界；启动期场景为 T12b 的口径）", "T9_LB", "T9_UB_old", "T9_UB_new")):
        L += ["", f"### {tag}：上界与差距，旧 / 新\n", "| m | 族 | 场景 | 实例 | 下界 | 上界（旧） | 上界（新） | 差距（旧） | 差距（新） | 上界变化的实例 |", "|---|---|---|---|---|---|---|---|---|---|"]
        for (m, fam, rs), xs in sorted(by.items()):
            ok = [x for x in xs if None not in (x[lb], x[uo], x[un])]
            L.append(f"| {m} | {fam} | {('启动期' if lb == 'T9_LB' else '允许') if rs else '不允许'} | {len(ok)} | {fm(mean([x[lb] for x in ok]))} | {fm(mean([x[uo] for x in ok]))} | {fm(mean([x[un] for x in ok]))} | "
                     f"{fm(mean([x[uo] - x[lb] for x in ok]))} | {fm(mean([x[un] - x[lb] for x in ok]))} | {sum(x[un] != x[uo] for x in ok)} |")
        ch = [x for x in rows if x[un] != x[uo]]
        if ch: L.append("\n上界变化：" + "；".join(f"{x['inst']} {'允许' if x['reset'] else '不允许'} {x[uo]} → {x[un]}" for x in ch))
    L += ["", "### T8：比 ES 饱和值少 1 个的 (实例, 规则)\n", "| 实例 | 规则 | RDD 最好 | ES 饱和值（旧） | ES 饱和值（新） | 饱和值来源（新） | 仍然更少 |", "|---|---|---|---|---|---|---|"]
    for n, rule, rb in T8_BELOW:
        r = next((x for x in rows if x["inst"] == n and x["reset"] == 0), None)
        if r is None: continue    # 该实例还没跑完
        src = r["T4b_src"] or {}
        L.append(f"| {n} | {rule} | {rb} | {r['T4b_UB_old']} | {r['T4b_UB_new']} | {src.get('kind')} {src.get('src')} s{src.get('seed', '')} | {'是' if rb < r['T4b_UB_new'] else '否'} |")
    L += ["", f"核对：旧上界与 T4 / T4b / T9 汇总相同，不一致 {len(bad)} 处" + (f"：{bad}" if bad else "") + f"；实例 {len(J)} / 24。"]
    print("\n".join(L))


if __name__ == "__main__":
    main()
