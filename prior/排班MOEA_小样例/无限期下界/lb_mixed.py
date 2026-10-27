import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import math, itertools
from roster_sim import INSTANCES
from lb import robot_occ_seq, tool_occ, max_shifts_cont

def mixed_lb(name, H=None):
    inst = INSTANCES[name](); H = H or inst.horizon
    R = inst.robots
    info = []
    for r, p in enumerate(R):
        seq, mu = robot_occ_seq(inst, r)
        kc = max_shifts_cont(seq, mu, H)
        o = min(tool_occ(inst, r).values())
        info.append(dict(d=p["duty"], m=p["tool_num"], kc=kc, o=o, cap1=p["duty"]*kc, lim=p["max_num"]))
    o_min = min(i["o"] for i in info)
    # 单把工具在 H 内最多被使用次数(连续松弛, 每次占用 >= o_min)
    u_max = math.ceil(H / o_min - 1e-9)
    best = None
    for N in itertools.product(*[range(0, i["lim"]+1) for i in info]):
        cap = sum(n*i["cap1"] for n, i in zip(N, info))
        if cap < H - 1e-9: continue
        # 用工具-次数/小时最小的类型优先覆盖
        need, uses = H, 0.0
        for n, i in sorted(zip(N, info), key=lambda z: z[1]["m"]/z[1]["d"]):
            x = min(need, n*i["cap1"]); uses += x / i["d"] * i["m"]; need -= x
            if need <= 1e-9: break
        M = math.ceil(uses / u_max - 1e-9)
        tot = sum(N) + M
        if best is None or tot < best[0]: best = (tot, N, M)
    return best

for n in ["L1-720", "L2-720", "L3-720", "L3-720紧", "L1-1440", "L3-1440"]:
    print(n, "混合类型连续下界 (f1, 各型机器人数, 工具数):", mixed_lb(n))
