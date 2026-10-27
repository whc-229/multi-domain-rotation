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

out = {}
mode = sys.argv[1]
if mode == "A":
    cases = [("L1-720", 0, 7), ("L2-720", 0, 5), ("L3-720", 0, 6), ("L1-720", 2, 8), ("L1-720", 2, 7)]
    for name, r, N in cases:
        R = homogeneous(name, r, N); key = f"{name}/{r}型/N={N}"
        t = time.time(); ex = kpost.exact(R, 1, 3_000_000); ex = dict(ok=ex[0], nodes=ex[1], secs=time.time() - t)
        row = dict(exact=ex, full={}, cyc=[])
        for H in [144, 720, 2880]:
            succ = []; ev = []
            for s in range(5):
                ok, used = kpost.fullkey(R, 1, H, 6000, np.random.default_rng(s)); succ.append(bool(ok)); ev.append(int(used))
            row["full"][H] = dict(succ=sum(succ), evals=ev)
        for s in range(5):
            row["cyc"].append(run_cyc(R, 1, np.random.default_rng(s)))
        out[key] = row; print(key, json.dumps(row, ensure_ascii=False, default=str), flush=True)
else:
    for K in [1, 2, 3]:
        for ps in range(3):
            pool = sorted(hetero_pool(8 * K, 100 + 10 * K + ps), key=lambda r: -rate(r))
            lb = lbK(pool, K); res = dict(LB=lb)
            for meth in ["exact", "cyc", "full"]:
                best = None; log = []
                for N in range(lb, min(len(pool), lb + 4) + 1):
                    R = pool[:N]; t = time.time()
                    if meth == "exact":
                        ok, nodes = kpost.exact(R, K, 300_000); info = dict(ok=ok, nodes=nodes)
                    elif meth == "cyc":
                        info = run_cyc(R, K, np.random.default_rng(ps)); ok = info["ok"]
                    else:
                        ok, used = kpost.fullkey(R, K, 144, 3000, np.random.default_rng(ps)); info = dict(ok=bool(ok), evals=int(used))
                    info["N"] = N; info["secs"] = round(time.time() - t, 1); log.append(info)
                    if ok: best = N; break
                res[meth] = dict(best=best, log=log)
            out[f"K={K}/池{ps}"] = res
            print(f"K={K} 池{ps} LB={lb} 精确={res['exact']['best']} 循环编码+前缀={res['cyc']['best']} 全时域GA(有限144槽)={res['full']['best']}",
                  json.dumps(res, ensure_ascii=False, default=str), flush=True)
json.dump(out, open(f"exp{mode}.json", "w"), ensure_ascii=False, indent=1, default=str)
