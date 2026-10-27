import time, json, sys
import numpy as np
from slotmodel import homogeneous, hetero_pool, rate
import kpost

def lbK(pool, K):
    rs = sorted((rate(r) for r in pool), reverse=True); s = 0
    for k, v in enumerate(rs, 1):
        s += v
        if s >= K - 1e-12: return k

def run_cyc(R, K, rng, budget=20000, pnodes=100000, tries=5):
    """循环编码找周期 → 束搜索/深搜找前缀；若该周期从全新状态进不去，就再找另一个周期（总预算不变）"""
    t = time.time(); used_tot = 0; nodes_tot = 0; why = "未找到周期"
    for _ in range(tries):
        if used_tot >= budget: break
        res, used = kpost.cyclic(R, K, budget - used_tot, rng); used_tot += used
        if res is None: break
        seq, phis, P = res; tg = kpost.cycle_targets(R, K, seq, phis, P)
        assert kpost.verify_cycle(R, K, seq, phis, P, tg[0], 0)
        used_r = sorted(set(int(x) for x in seq))
        ok, nodes = kpost.prefix(R, K, tg, pnodes // tries, used_r); nodes_tot += nodes
        if ok: return dict(ok=True, evals=used_tot, nodes=nodes_tot, P=P, secs=round(time.time() - t, 1), why="")
        why = "前缀未找到"
    return dict(ok=False, evals=used_tot, nodes=nodes_tot, secs=round(time.time() - t, 1), why=why)

