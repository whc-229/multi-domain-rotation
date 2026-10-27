"""结构编码 + 组内精确：先用构造公式做分组动态规划，再对选中的每个组（单一机型、组内站点）跑周期 ILP（P 取全表、不提前停），
组台数取 min(公式, ILP)。ILP 解本身就是周期排班，所以总数是可行上界、天然无限期可行。"""
import json, sys, time
from gen import make, SIX
from mtms import Inst, periodic_ilp, fluid_lb
from struct_enc import solve
for spec in sys.argv[1].split(","):
    t0 = time.time()
    m, sd = spec.split("s"); I = make(int(m), SIX, int(sd)); v, groups = solve(I); t_struct = time.time() - t0
    tot = 0; detail = []
    for G, tn, ub, lb in groups:
        k = I.tn.index(tn); J = Inst({tn: I.T[k]}, [I.sites[s] for s in G]); best = ub
        for P in [15, 30, 45, 60, 90, 120, 150]:
            r = periodic_ilp(J, P, time_limit=60)
            if r["obj"] is not None: best = min(best, round(r["obj"]))
        tot += best; detail.append((G, tn, ub, best))
    print(json.dumps(dict(inst=I.name, intLB=round(fluid_lb(I, True)[0]), struct=v, struct_sec=round(t_struct, 1), refined=tot,
                          refine_sec=round(time.time() - t0), groups=detail), ensure_ascii=False), flush=True)
