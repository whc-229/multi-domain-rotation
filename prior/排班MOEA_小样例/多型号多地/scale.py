"""规模放大：m 个站点（6 种真实机型，每站限定 2–3 型）。比较
- 整数下界 intLB（全局共享）；
- 各站独立配机的精确解（每站单独周期 ILP，周期取若干值的最好结果，求和；不跨站共享）；
- 派工规则 FIFO、相位平衡 PB（全局共享、仿真定台数）；
- 进化策略 ES（全局共享）。
所有仿真解在 H=20000 格（约 4.6 年）上复核无缺岗。"""
import json, sys, time
from gen import make, SIX
from mtms import fluid_lb, periodic_ilp, Inst
from sim import size_fleet, simulate
from evo_policy import evolve

def single_site(I, s, Ps, tl):
    sub = {t: I.T[I.tn.index(t)] for t in I.sites[s]["E"]}
    J = Inst(sub, [I.sites[s]]); best = None
    for P in Ps:
        r = periodic_ilp(J, P, time_limit=tl)
        if r["obj"] is not None and (best is None or r["obj"] < best): best = round(r["obj"])
    return best

ms = [int(x) for x in sys.argv[1].split(",")]; seeds = [int(x) for x in sys.argv[2].split(",")]
gens = int(sys.argv[3]) if len(sys.argv) > 3 else 20
for m in ms:
    for sd in seeds:
        I = make(m, SIX, sd); rec = dict(inst=I.name)
        t0 = time.time(); rec["intLB"] = round(fluid_lb(I, True)[0]); rec["fluid"] = round(fluid_lb(I)[0], 2)
        t0 = time.time(); ded = [single_site(I, s, [30, 60, 90], 40) for s in range(m)]
        rec["dedicated_exact"] = sum(ded) if None not in ded else None; rec["ded_sec"] = round(time.time() - t0)
        for pol in ["fifo", "pb"]:
            N, Nk = size_fleet(I, pol, 3000, 1000)
            g, _ = simulate(I, Nk, pol, 20000, 1000) if N else (None, None)
            rec[pol] = N; rec[pol + "_verify_gaps"] = g
        t0 = time.time(); N, Nk, w, hist = evolve(I, gens=gens)
        g, _ = simulate(I, Nk, "score", 20000, 1000, params=w)
        rec["es"] = N; rec["es_Nk"] = Nk; rec["es_verify_gaps"] = g; rec["es_sec"] = round(time.time() - t0); rec["es_w"] = [round(x, 3) for x in w]
        print(json.dumps(rec, ensure_ascii=False), flush=True)
