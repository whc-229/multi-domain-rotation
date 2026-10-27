"""复核结构编码的每个分组：用周期网络 ILP 在"该机型 + 这组站点"子实例上找实际可行排班，确认构造给出的台数可以达到（ILP 找到 ≤ 该值的解即确认）。"""
import json, sys
from gen import make, SIX
from mtms import Inst, periodic_ilp
from struct_enc import solve
for spec in sys.argv[1].split(","):
    m, sd = spec.split("s"); I = make(int(m), SIX, int(sd)); v, groups = solve(I)
    for G, tn, ub, lb in groups:
        k = I.tn.index(tn)
        J = Inst({tn: I.T[k]}, [I.sites[s] for s in G])
        best = None
        for P in [15, 30, 45, 60, 90, 150]:
            r = periodic_ilp(J, P, time_limit=60)
            if r["obj"] is not None and (best is None or r["obj"] < best): best = round(r["obj"])
            if best is not None and best <= ub: break
        print(json.dumps(dict(inst=I.name, group=G, type=tn, struct=ub, lb=lb, ilp_best=best, ok=(best is not None and best <= ub))), flush=True)
