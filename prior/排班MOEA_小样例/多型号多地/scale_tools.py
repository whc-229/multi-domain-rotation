"""多地 + 共享工具池：机器人 + 工具最少总数。比较 下界、FIFO、PB、ES（全局共享）与"各站独立"（每站自带机器人和工具，各自 ES）。
稳态评估：H=3000 格、预热 1000 格；最终解在 H=20000 上复核缺岗数。"""
import json, sys, time
from gen import make, SIX
from mtms import Inst
from tools_model import ToolInst, tool_lb, size_resources, simulate_tools
from evo_tools import evolve_tools

ms = [int(x) for x in sys.argv[1].split(",")]; seeds = [int(x) for x in sys.argv[2].split(",")]
gens = int(sys.argv[3]) if len(sys.argv) > 3 else 15
H, W, HV = 3000, 1000, 20000
for m in ms:
    for sd in seeds:
        I = make(m, SIX, sd); TI = ToolInst(I); rec = dict(inst=I.name)
        lb = tool_lb(TI); rec["LB"] = round(lb[0], 2)
        for pol in ["fifo", "pb"]:
            R, Nk, Mj = size_resources(TI, pol, H, W)
            rec[pol] = R
            rec[pol + "_verify"] = simulate_tools(TI, Nk, Mj, pol, HV, W)[0] if R else None
        t0 = time.time(); R, Nk, Mj, w, hist = evolve_tools(TI, H, W, gens=gens, fresh=True)
        rec["es"] = R; rec["es_robots"] = Nk; rec["es_tools"] = Mj; rec["es_sec"] = round(time.time() - t0)
        rec["es_verify"] = simulate_tools(TI, Nk, Mj, "score", HV, W, params=w)[0] if R else None
        # 各站独立：每站自己的机器人和工具
        t0 = time.time(); ded = 0
        for s in range(m):
            sub = {t: I.T[I.tn.index(t)] for t in I.sites[s]["E"]}
            J = Inst(sub, [I.sites[s]]); TJ = ToolInst(J)
            r, _, _, _, _ = evolve_tools(TJ, H, W, gens=max(4, gens // 3))
            ded += r if r else 999
        rec["dedicated_es"] = ded; rec["ded_sec"] = round(time.time() - t0)
        print(json.dumps(rec, ensure_ascii=False), flush=True)
