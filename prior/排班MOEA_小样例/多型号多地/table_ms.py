"""多站带工具实例：稳态下界、全新启动下界（不许/允许提前大修）、各方法最好的长期有效解。"""
import json, math, sys
from collections import defaultdict
from export_inst import build
from tools_model import ToolInst, tool_lb
valid = lambda d: d.get("R") is not None and d["R"] >= 0 and (d["verify_cert"] == 1 or (d["verify_cert"] == -1 and d["verify_gaps"] == 0))
LB = {}
for l in open("freshlb_all.jsonl"):
    d = json.loads(l); LB[d["spec"]] = d
best = defaultdict(dict); runs = defaultdict(lambda: defaultdict(list))
files = ["stats_small_certs.jsonl", "stats_m6_certs.jsonl", "stats_m8.jsonl", "stats_m10.jsonl", "stats4_small.jsonl", "stats4_small_pad.jsonl", "stats4_lpg.jsonl"]
for fn in files:
    for l in open(fn):
        d = json.loads(l)
        if "error" in d or d.get("gp") or d.get("mode") != "certs" or not valid(d): continue
        c = d.get("cfg", "base")
        if c in ("base", "curdef", "pad"): c = "无提前大修"
        best[d["spec"]][c] = min(best[d["spec"]].get(c, 10**9), d["R"]); runs[d["spec"]][c].append(d["R"])
print("| 实例 | 稳态下界 | 全新启动下界 | 允许提前大修的全新启动下界 | 进化派工最好（不提前大修） | LP 计划 + 进化派工最好 | 后者对允许大修下界的差距 |")
print("|---|---|---|---|---|---|---|")
order = ["L1", "L2", "L3"] + [f"{m}s{s}" for m in [2, 3, 4, 6, 8, 10] for s in [1, 2, 3, 4]]
for sp in order:
    if sp not in LB: continue
    d = LB[sp]; I = build(sp); ss = tool_lb(ToolInst(I), True)[0]
    lp = math.ceil(d["lp"] - 1e-6); lpR = math.ceil(d["lpR"] - 1e-6) if d.get("lpR") else None
    b0 = best[sp].get("无提前大修"); b1 = best[sp].get("lpgplan")
    gap = f"{100*(b1-lpR)/lpR:.1f}%" if b1 and lpR else "—"
    print(f"| {sp} | {ss:.0f} | {lp} | {lpR if lpR else '—'} | {b0 if b0 else '—'} | {b1 if b1 else '—'} | {gap} |")
