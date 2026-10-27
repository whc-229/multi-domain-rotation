"""由全新启动 LP（带工具、允许提前大修）生成提前大修计划文件：每行 k i t0 c。只保留前 cut 格（默认 H/2，去掉时域末端的伪影）。
用法：python3 make_plan.py <spec 或 单型实例文件名前缀> <输出> [H倍数=4] [tools=1]"""
import sys, json
from export_inst import build
from freshlb import fresh_lb, round_plan, S_max
from mtms import Inst, real_types
spec, out = sys.argv[1], sys.argv[2]
mult = float(sys.argv[3]) if len(sys.argv) > 3 else 4
tools = (sys.argv[4] != "0") if len(sys.argv) > 4 else True
if "-x" in spec:   # 单型、多岗位、只算机器人
    tn, K = spec.split("-x"); T = real_types()
    I = Inst({tn: T[tn]}, [dict(E=[tn], delta=0, K=1) for _ in range(int(K))], name=spec); tools = False
else:
    I = build(spec)
S = S_max(I); H = int(mult * S)
r = fresh_lb(I, H, tools=tools, reset=True, want_plan=True)
plan = round_plan(r["reset_flow"], H_cut=H // 2)
with open(out, "w") as f:
    for k, i, D, c in plan: f.write(f"{k} {i} {D} {c}\n")
print(json.dumps(dict(spec=spec, H=H, lpR=r["obj"], n_resets=sum(c for *_, c in plan), entries=len(plan), sec=r["sec"])))
