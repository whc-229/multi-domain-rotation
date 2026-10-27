"""启动代价试点驱动。对每个实例：流体下界、稳态（P∈{60,120,180} 取最小）、全新启动下界（H=4S 有限时域）与上界（可延拓尾部），
各自分"不许 / 允许提前大修"两种。用法：python3 run_startup.py single|multi 输出.jsonl"""
import sys, json, time, os
from mtms import Inst, real_types, fluid_lb
from gen import make, SIX
from startup import steady, fresh, G_of

def S_of(I):   # 各型超周期格数的最大值（按零路程站点计）
    return max(sum(-(-int(round(o * 10)) // 20) for o in t["occ_h"]) for t in I.T)

def best(vals):
    v = [x for x in vals if x is not None]
    return min(v) if v else None

def periods(S):
    c = lambda x: -(-x // 15) * 15     # 15 = 值岗 5 格与 3 格的最小公倍数
    return sorted({60, 120, 180, c(S), c(2 * S)})

def run(I, tl=120):
    S = S_of(I); out = dict(inst=I.name, types=I.tn, K=[s["K"] for s in I.sites], S=S, fluid=round(fluid_lb(I, True)[0], 3))
    PS = periods(S)
    for reset in [False, True]:
        tag = "R" if reset else ""
        t0 = time.time()
        out["steady" + tag] = best(steady(I, P, reset, tl) for P in PS)
        lb = fresh(I, 4 * S, tail=False, reset=reset, time_limit=tl)
        out["freshLB" + tag] = lb["obj"] if lb["obj"] is not None and lb["status"] == 0 else lb["bound"]
        out["freshLB_exact" + tag] = lb["status"] == 0
        G = G_of(I, reset); ubs = []
        for T in [G + S, G + 2 * S]:
            for P in PS:
                r = fresh(I, T, P, reset=reset, time_limit=tl)
                if r["obj"] is not None: ubs.append(r["obj"])
        out["freshUB" + tag] = best(ubs)
        out["sec" + tag] = round(time.time() - t0)
    return out

if __name__ == "__main__":
    mode, fn = sys.argv[1], sys.argv[2]
    done = set()
    if os.path.exists(fn):
        for l in open(fn): d = json.loads(l); done.add((d["inst"], tuple(d["K"])))
    T = real_types(); jobs = []
    if mode == "single":
        for tn in SIX:
            for K in [1, 2, 3, 4]:
                jobs.append(Inst({tn: T[tn]}, [dict(E=[tn], delta=0, K=K)], name=f"{tn}-K{K}"))
    else:
        for m in [2, 3]:
            for sd in [1, 2, 3, 4]: jobs.append(make(m, SIX, sd))
    for I in jobs:
        if (I.name, tuple(s["K"] for s in I.sites)) in done: continue
        with open(fn, "a") as f: f.write(json.dumps(run(I), ensure_ascii=False) + "\n")
