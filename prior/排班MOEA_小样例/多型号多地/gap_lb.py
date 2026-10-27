"""各实例：稳态流体下界、全新启动 LP 下界（不许/允许提前大修）与进化解（长期有效）的差距。
用法：python3 gap_lb.py freshlb_4S.jsonl 结果文件..."""
import json, sys, math
from collections import defaultdict
from export_inst import build
from tools_model import ToolInst, tool_lb
valid = lambda d: d.get("R") is not None and d["R"] >= 0 and (d["verify_cert"] == 1 or (d["verify_cert"] == -1 and d["verify_gaps"] == 0))
LB = {json.loads(l)["spec"]: json.loads(l) for l in open(sys.argv[1])}
best = defaultdict(dict)
for fn in sys.argv[2:]:
    for l in open(fn):
        d = json.loads(l)
        if "error" in d or not valid(d): continue
        if d.get("mode") not in ("certs",) or d.get("gp"): continue
        c = d.get("cfg", "base")
        sp = d["spec"]
        best[sp][c] = min(best[sp].get(c, 10**9), d["R"])
print(f"{'实例':6s} {'稳态LB':>6s} {'全新LB':>6s} {'全新LB(大修)':>10s} | 各配置最好 R | 对稳态LB差距 | 对全新LB差距")
for sp, d in LB.items():
    I = build(sp); ss = tool_lb(ToolInst(I), True)[0]
    lp, lpR = math.ceil(d["lp"] - 1e-6), math.ceil(d["lpR"] - 1e-6)
    b = best.get(sp, {})
    if not b: print(f"{sp:6s} {ss:6.0f} {d['lp']:6.2f} {d['lpR']:10.2f} | —"); continue
    bb = min(b.values()); bnr = min(v for k, v in b.items() if k != "reset") if any(k != "reset" for k in b) else None
    s = " ".join(f"{k}={v}" for k, v in sorted(b.items()))
    print(f"{sp:6s} {ss:6.0f} {d['lp']:6.2f} {d['lpR']:10.2f} | {s} | {100*(bb-ss)/ss:5.1f}% | 不许大修 {100*(bnr-lp)/lp if bnr else float('nan'):5.1f}%，允许 {100*(bb-lpR)/lpR:5.1f}%")
