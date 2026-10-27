"""规模扫描：站点数 m、随机种子；比较
- 流体下界、整数台数下界；
- 各站点独立配机（不跨站共享）：每站单独求周期 ILP，周期 P 取若干值的最好结果，求和；
- 跨站共享：整体周期 ILP（P ∈ {30, 60}），记录目标、对偶界、用时。"""
import json, sys, time
from gen import make, SIX, NINE
from mtms import fluid_lb, periodic_ilp, Inst

def single_site(I, s, Ps, tl):
    sub_types = {t: I.T[I.tn.index(t)] for t in I.sites[s]["E"]}
    J = Inst(sub_types, [I.sites[s]], cost={t: I.cost[t] for t in sub_types})
    best = None; t0 = time.time()
    for P in Ps:
        r = periodic_ilp(J, P, time_limit=tl)
        if r["obj"] is not None and (best is None or r["obj"] < best): best = round(r["obj"])
    return best, round(time.time() - t0, 1)

ms = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [2, 3, 4, 6]
seeds = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [1, 2, 3]; TL = int(sys.argv[2]) if len(sys.argv) > 2 else 300
for m in ms:
    for sd in seeds:
        I = make(m, SIX, sd)
        rec = dict(inst=I.name, sites=[(s["E"], s["delta"]) for s in I.sites])
        rec["fluid"] = round(fluid_lb(I)[0], 3); rec["intLB"] = round(fluid_lb(I, True)[0])
        ded = [single_site(I, s, [30, 60, 90, 150], 60) for s in range(m)]
        rec["dedicated"] = sum(d[0] for d in ded) if all(d[0] is not None for d in ded) else None
        rec["dedicated_sec"] = sum(d[1] for d in ded)
        for P in [30, 60]:
            t0 = time.time(); r = periodic_ilp(I, P, time_limit=TL)
            rec[f"shared_P{P}"] = dict(obj=None if r["obj"] is None else round(r["obj"], 3),
                                       bound=None if r["bound"] is None else round(r["bound"], 3), sec=round(time.time() - t0, 1), nvar=r.get("nvar"))
        print(json.dumps(rec, ensure_ascii=False), flush=True)
