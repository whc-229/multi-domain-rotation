import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
"""有限时域 H 内：同型 N 台能否覆盖全部 S 个班次槽（精确 DFS，带失败状态记忆）"""
import math, sys
from functools import lru_cache
from roster_sim import INSTANCES
from lb import robot_occ_seq
sys.setrecursionlimit(100000)
def feasible(name, r, N):
    inst = INSTANCES[name](); d = inst.robots[r]["duty"]
    seq, mu = robot_occ_seq(inst, r); g = [math.ceil(o/d-1e-9) for o in seq]; L = len(g)
    nph = lambda ph: ph+1 if ph+1 < L else mu
    S = math.ceil(inst.horizon / d)
    @lru_cache(maxsize=None)
    def ok(s, t):
        if t == S: return True
        for (w, ph) in set(s):
            if w == 0:
                lst = list(s); lst.remove((w, ph))
                new = tuple(sorted([(max(w2-1,0),p2) for (w2,p2) in lst] + [(g[ph]-1, nph(ph))]))
                if ok(new, t+1): return True
        return False
    return ok(tuple([(0,0)]*N), 0)
for name, N in [("L1-720",6),("L1-1440",6),("L3-720",5),("L3-1440",5)]:
    print(name, "0型", N, "台能否覆盖整个时域:", feasible(name, 0, N))
