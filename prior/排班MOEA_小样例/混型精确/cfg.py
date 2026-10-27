import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
"""对每个实例：从混型速率下界开始，按总资源数递增枚举所有配置（各型机器人数 + 各型工具数），
用精确搜索判定无限期 / 有限期 H 可行，得到真正的最少资源数（含混型）。"""
import itertools, time, sys, json
from roster_sim import INSTANCES
from lb import robot_occ_seq
from mixed_exact import Model, infinite_feasible, finite_feasible

def rho(inst, r):
    seq, mu = robot_occ_seq(inst, r); cyc = seq[mu:]
    return inst.robots[r]["duty"] * len(cyc) / sum(cyc)

def configs(inst, total, rate_prune=True):
    nR, nTl = len(inst.robots), len(inst.tools)
    rh = [rho(inst, r) for r in range(nR)]
    for nr_tot in range(2, total):
        for nr in itertools.product(range(nr_tot + 1), repeat=nR):
            if sum(nr) != nr_tot or any(nr[r] > inst.robots[r]["max_num"] for r in range(nR)): continue
            if rate_prune and sum(n * x for n, x in zip(nr, rh)) < 1 - 1e-9: continue          # 速率下界剪枝（只对无限期有效）
            used = [r for r in range(nR) if nr[r]]
            ok_t = sorted({t for r in used for t in inst.robots[r]["tools"]})
            nt_tot = total - nr_tot
            if nt_tot < min(inst.robots[r]["tool_num"] for r in used) * 2: continue   # 相邻两班工具必不同(已证)
            for split in itertools.product(range(nt_tot + 1), repeat=len(ok_t)):
                if sum(split) != nt_tot: continue
                nt = [0] * nTl
                for t, c in zip(ok_t, split): nt[t] = c
                if any(nt[t] > inst.tools[t]["max_num"] for t in range(nTl)): continue
                # 每个用到的机器人型号都要有足够的适配工具
                if any(sum(nt[t] for t in inst.robots[r]["tools"]) < inst.robots[r]["tool_num"] for r in used): continue
                yield nr, tuple(nt)

