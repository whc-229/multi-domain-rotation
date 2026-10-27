"""L3 系列机型（超周期 130–190 格）的单站精确试点：周期取 ≥ S 的 15 的倍数。K = 1..3。"""
import sys, json, time
from mtms import Inst, real_types, fluid_lb
from startup import steady, fresh, G_of
from run_startup import S_of, best
c = lambda x: -(-x // 15) * 15
T = real_types(); fn = sys.argv[1]
for tn in ["L3-0", "L3-1", "L3-2"]:
    for K in [1, 2, 3]:
        I = Inst({tn: T[tn]}, [dict(E=[tn], delta=0, K=K)], name=f"{tn}-K{K}")
        S = S_of(I); PS = [c(S), c(2 * S)]
        out = dict(inst=I.name, types=I.tn, K=[K], S=S, fluid=round(fluid_lb(I, True)[0], 3))
        for reset in [False, True]:
            tag = "R" if reset else ""; t0 = time.time()
            out["steady" + tag] = best(steady(I, P, reset, 180) for P in PS)
            lb = fresh(I, 4 * S, tail=False, reset=reset, time_limit=180)
            out["freshLB" + tag] = lb["obj"] if lb["obj"] is not None and lb["status"] == 0 else lb["bound"]
            out["freshLB_exact" + tag] = lb["status"] == 0
            G = G_of(I, reset)
            out["freshUB" + tag] = best((fresh(I, G + S, P, reset=reset, time_limit=180) or {}).get("obj") for P in PS)
            out["sec" + tag] = round(time.time() - t0)
        with open(fn, "a") as f: f.write(json.dumps(out, ensure_ascii=False) + "\n")
