"""T11 汇总：RDD-FD 与 T8（完整定台数的 RDD）在 m ∈ {4, 6} 上配对比较（同实例、同默认规则、同种子，不提前大修场景）。
输入：results/T11/jobs/*.json，results/T8/jobs/<实例>__<规则>__r0.1__s<种子>.json。
输出：results/T11_rddfd.jsonl（逐作业，含配对的 T8 结果）与 markdown 表：
  1. 判据：按 (m, 规则) 的资源数平均差（T11 − T8）、胜负、用时比（T8 墙钟 / T11 墙钟）；
  2. 判定：T_eval 内 0 缺岗的候选里完整判定不可行的比例（判定失败率）；
  3. 核对：最终解重新判定逐字段相同的个数；抽样缺岗与 eval 相同的个数；默认规则的台数与 T8 相同的个数；
  4. 敏感性：te4 / noswap / selbase 对 main 的配对差。
可行 = 严格证明或 60 万格复核 0 缺岗。资源数只在两边都可行时比较。
用法：python3 code/t11_summary.py"""
import os, sys, json, glob
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402


def feas(d):
    """可行时返回资源数，否则 None。"""
    if not d or d.get("R") is None or d["R"] < 0: return None
    ok = (d.get("verify_cert") == 1 and d.get("exact") == 1) or d.get("verify_gaps") == 0
    return d["R"] if ok else None


def mean(xs):
    """忽略 None 的平均。"""
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def fm(x, nd=2):
    """格式化。"""
    return "—" if x is None else f"{x:.{nd}f}"


def main():
    J = {}
    for p in glob.glob(os.path.join(ROOT, "results", "T11", "jobs", "*.json")):
        d = json.load(open(p)); J[(d["inst"], d["rule"], d["seed"], d["setting"])] = d
    rows = []
    for (n, r, s, st), d in sorted(J.items()):
        t8 = json.load(open(os.path.join(ROOT, "results", "T8", "jobs", f"{n}__{r}__r0.1__s{s}.json")))
        rows.append(dict(inst=n, m=d["m"], family=d["family"], rule=r, seed=s, setting=st, R=feas(d), status=d["status"], proven=d["proven"],
                         base_R=d["base_R"] if d["base_R"] >= 0 else None, base_same_t8=int(d["base_R"] == t8["base_R"] and d["base_Nk"] == t8["base_Nk"] and d["base_Mj"] == t8["base_Mj"]),
                         t8_R=feas(t8), t8_proven=int(t8.get("proven") == 1), wall=d["wall"], t8_wall=t8["wall"], gens=d["gens"], evals=d["evals"],
                         judges=d["judges"], judge_fail=d["judge_fail"], judge_quick=d.get("judge_quick"), dec_steps=d["dec_steps"], backtracks=d["backtracks"],
                         sec_eval=d["sec_eval"], sec_judge=d["sec_judge"], sec_dec=d["sec_dec"], base_sec=d["base_sec"],
                         check_judge=d.get("check_judge"), check_dump=(d["check_dump_same"], d["check_dump_n"]), traj=d["traj"]))
    tmp = os.path.join(ROOT, "results", "T11_rddfd.jsonl.tmp")
    with open(tmp, "w") as f:
        for r in rows: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(tmp, os.path.join(ROOT, "results", "T11_rddfd.jsonl"))
    by = defaultdict(list)
    for r in rows: by[(r["setting"], r["m"], r["rule"])].append(r)
    L = ["### 判据：RDD-FD（main）对 T8，同实例、同规则、同种子配对\n",
         "| m | 规则 | 作业 | 两边都可行 | 默认规则 | T8 | RDD-FD | 平均差（FD − T8） | FD 少 / 同 / 多 | T8 平均用时 s | FD 平均用时 s | 用时比（逐作业 T8/FD 的中位数） | 最小用时比 | FD 严格证明 / T8 严格证明 |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    allp = []
    for (st, m, rule), xs in sorted(by.items()):
        if st != "main": continue
        pr = [x for x in xs if x["R"] is not None and x["t8_R"] is not None]; allp += pr
        ratio = sorted(x["t8_wall"] / max(x["wall"], 0.1) for x in xs)
        L.append(f"| {m} | {rule} | {len(xs)} | {len(pr)} | {fm(mean([x['base_R'] for x in pr]))} | {fm(mean([x['t8_R'] for x in pr]))} | {fm(mean([x['R'] for x in pr]))} | "
                 f"{fm(mean([x['R'] - x['t8_R'] for x in pr]))} | {sum(x['R'] < x['t8_R'] for x in pr)} / {sum(x['R'] == x['t8_R'] for x in pr)} / {sum(x['R'] > x['t8_R'] for x in pr)} | "
                 f"{fm(mean([x['t8_wall'] for x in xs]), 0)} | {fm(mean([x['wall'] for x in xs]), 1)} | {ratio[len(ratio) // 2]:.0f} | {ratio[0]:.1f} | "
                 f"{sum(x['proven'] for x in pr)} / {sum(x['t8_proven'] for x in pr)} |")
    tw8 = sum(x["t8_wall"] for x in rows if x["setting"] == "main"); tw = sum(x["wall"] for x in rows if x["setting"] == "main")
    L.append(f"\n全部 main：{len(allp)} 对，平均差 {fm(mean([x['R'] - x['t8_R'] for x in allp]))}；"
             f"FD 少 / 同 / 多 = {sum(x['R'] < x['t8_R'] for x in allp)} / {sum(x['R'] == x['t8_R'] for x in allp)} / {sum(x['R'] > x['t8_R'] for x in allp)}；"
             f"总用时 T8 {tw8 / 3600:.1f} h，FD {tw / 3600:.2f} h（{tw8 / tw:.0f} 倍）。")
    L += ["", "### 判定与用时构成（main）\n",
          "| m | 规则 | 平均代数 | 平均评价次数 | 完整判定（含快速否决）次数 | 其中不可行 | 判定失败率 | 换资源（含回退） | 默认规则定台数 s | T_eval 评价 s | 判定 s | 选资源 s |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for (st, m, rule), xs in sorted(by.items()):
        if st != "main": continue
        jn = sum(x["judges"] for x in xs); jf = sum(x["judge_fail"] for x in xs)
        L.append(f"| {m} | {rule} | {fm(mean([x['gens'] for x in xs]), 0)} | {fm(mean([x['evals'] for x in xs]), 0)} | {jn} | {jf} | {fm(100 * jf / jn if jn else None, 1)}% | "
                 f"{fm(mean([x['dec_steps'] for x in xs]), 1)}（{fm(mean([x['backtracks'] for x in xs]), 1)}） | {fm(mean([x['base_sec'] for x in xs]), 1)} | "
                 f"{fm(mean([x['sec_eval'] for x in xs]), 1)} | {fm(mean([x['sec_judge'] for x in xs]), 1)} | {fm(mean([x['sec_dec'] for x in xs]), 2)} |")
    L += ["", "### 核对\n"]
    for st in ("main", "te4", "noswap", "selbase"):
        xs = [x for x in rows if x["setting"] == st]
        if not xs: continue
        cj = [x["check_judge"] for x in xs if x["check_judge"] is not None]
        L.append(f"- {st}：{len(xs)} 个作业；最终解重新判定逐字段相同 {sum(cj)} / {len(cj)}；抽样缺岗与 eval 相同 {sum(x['check_dump'][0] for x in xs)} / {sum(x['check_dump'][1] for x in xs)}；"
                 f"默认规则台数与 T8 相同 {sum(x['base_same_t8'] for x in xs)} / {len(xs)}。")
    L += ["", "### 敏感性（对 main 配对，同实例、同规则、同种子）\n",
          "| 设置 | 改了什么 | m | 规则 | 对数 | 平均差（设置 − main） | 少 / 同 / 多 | 判定失败率 | 平均用时 s |", "|---|---|---|---|---|---|---|---|---|"]
    desc = {"te4": "T_eval = 4·S_max", "noswap": "不做台数互换（只减不换）", "selbase": "选资源用默认规则"}
    M = {(x["inst"], x["rule"], x["seed"]): x for x in rows if x["setting"] == "main"}
    for (st, m, rule), xs in sorted(by.items()):
        if st == "main": continue
        pr = [(x, M[(x["inst"], x["rule"], x["seed"])]) for x in xs if (x["inst"], x["rule"], x["seed"]) in M]
        pr = [(a, b) for a, b in pr if a["R"] is not None and b["R"] is not None]
        jn = sum(x["judges"] for x in xs); jf = sum(x["judge_fail"] for x in xs)
        L.append(f"| {st} | {desc[st]} | {m} | {rule} | {len(pr)} | {fm(mean([a['R'] - b['R'] for a, b in pr]))} | "
                 f"{sum(a['R'] < b['R'] for a, b in pr)} / {sum(a['R'] == b['R'] for a, b in pr)} / {sum(a['R'] > b['R'] for a, b in pr)} | "
                 f"{fm(100 * jf / jn if jn else None, 1)}% | {fm(mean([x['wall'] for x in xs]), 1)} |")
    print("\n".join(L))


if __name__ == "__main__":
    main()
