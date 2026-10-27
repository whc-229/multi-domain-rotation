"""LP 计划 + 进化派工 对 不提前大修的进化派工（curdef，同一随机数流）的配对比较；并对照全新启动 LP 下界。
用法：python3 analyze_plan.py <plan 结果> <curdef 结果> [freshlb_4S.jsonl]"""
import json, sys, math
from collections import defaultdict
from scipy.stats import wilcoxon
valid = lambda d: d.get("R") is not None and d["R"] >= 0 and (d["verify_cert"] == 1 or (d["verify_cert"] == -1 and d["verify_gaps"] == 0))
P = defaultdict(dict)
for fn in sys.argv[1:3]:
    for l in open(fn):
        d = json.loads(l)
        if "error" in d or d.get("cfg") not in ("plan", "curdef"): continue
        P[(d["spec"], d["seed"])][d["cfg"]] = d
LB = {}
if len(sys.argv) > 3:
    for l in open(sys.argv[3]): d = json.loads(l); LB[d["spec"]] = (d["lp"], d["lpR"])
dR = []; rows = defaultdict(list)
for (sp, sd), D in sorted(P.items()):
    if "plan" in D and "curdef" in D:
        a, b = D["plan"], D["curdef"]
        if valid(a) and valid(b): dR.append(a["R"] - b["R"])
        rows[sp].append((sd, b["R"], a["R"], a["resets"], a["cyc_resets"], valid(a), a["verify_cert"]))
for sp, L in rows.items():
    lb = LB.get(sp)
    lbs = f"  LB {math.ceil(lb[0]-1e-6)} / 允许大修 {math.ceil(lb[1]-1e-6)}" if lb else ""
    print(sp, " ".join(f"s{sd}: {c}→{p}（提前大修 {r} 次，环内 {cr}{'' if v else '，失效'}）" for sd, c, p, r, cr, v, vc in L) + lbs)
if dR:
    nz = [x for x in dR if x != 0]; p = wilcoxon(nz).pvalue if len(nz) >= 5 else float("nan")
    print(f"plan − curdef：{len(dR)} 对，均值 {sum(dR)/len(dR):+.2f}，少 {sum(x<0 for x in dR)} 同 {sum(x==0 for x in dR)} 多 {sum(x>0 for x in dR)}，p={p:.3g}")
