"""由全新启动 LP（允许提前大修）生成"按比例"的提前大修计划：首行 #N 为 LP 各型台数，其余每行 k i t0 流量（只保留前 H/2 格）。
仿真里按实际台数换算：k 型累计流量 × N_k(实际)/N_k(LP) 四舍五入。用法：python3 make_planf.py <spec> <输出> [H倍数=4]"""
import sys, json
from export_inst import build
from freshlb import fresh_lb, S_max
from mtms import Inst, real_types
spec, out = sys.argv[1], sys.argv[2]
mult = float(sys.argv[3]) if len(sys.argv) > 3 else 4
tools = True
if "-x" in spec:
    tn, K = spec.split("-x"); T = real_types()
    I = Inst({tn: T[tn]}, [dict(E=[tn], delta=0, K=1) for _ in range(int(K))], name=spec); tools = False
else:
    I = build(spec)
S = S_max(I); H = int(mult * S)
r = fresh_lb(I, H, tools=tools, reset=True, want_plan=True, time_limit=3600)
if r.get("N") is None: sys.exit(f"LP 未解出：{r}")
with open(out, "w") as f:
    f.write("#N " + " ".join(str(x) for x in r["N"]) + "\n")
    for k, i, D, x in r["reset_flow"]:
        if D < H // 2: f.write(f"{k} {i} {D} {x}\n")
import math
tot = {}
for key, x in r["duty_flow"].items():
    k, s = map(int, key.split(",")); tot[s] = tot.get(s, 0) + x
with open(out.replace(".f.txt", ".lpg.txt"), "w") as f:
    f.write("#N0 " + " ".join(str(int(math.floor(x + 1e-6))) for x in r["N"]) + "\n")
    f.write("#M0 " + " ".join(str(int(math.floor(x + 1e-6))) for x in r["M"]) + "\n")
    for key, x in r["duty_flow"].items():
        k, s = map(int, key.split(",")); f.write(f"z {k} {s} {x / tot[s]:.4f}\n")
print(json.dumps(dict(spec=spec, H=H, lpR=r["obj"], N=r["N"], M=r["M"], flow=round(sum(x for *_, x in r["reset_flow"]), 3), sec=r["sec"])), flush=True)
