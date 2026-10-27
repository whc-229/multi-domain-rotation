import os, sys; _h=os.path.dirname(os.path.abspath(__file__)); sys.path[:0]=[os.path.join(_h,'..','启发式验证'), os.path.join(_h,'..')]
"""D2：K 个值班点共享同型机器人。比较 共享最少台数 N*(K)、专点分组 K·N*(1)、速率下界 ⌈K·ā⌉。全新启动。"""
import math, sys, json, time, random
from slotmodel import Robot, homogeneous
import kpost

def nstar(g, mu, K, lim):
    a = sum(g[mu:]) / len(g[mu:]); lb = math.ceil(K * a - 1e-9)
    for N in range(lb, lb + 6):
        R = [Robot(g, mu, 0) for _ in range(N)]
        ok, nodes = kpost.exact(R, K, lim)
        if ok is None: return ("?", N)
        if ok: return N
    return None

real = {}
from slotmodel import gaps_from
from roster_sim import INSTANCES
for name in ["L1-720", "L2-720", "L3-720"]:
    for r in [0, 2]:
        g, mu = gaps_from(INSTANCES[name](), r); real[f"{name}/{r}型"] = (g, mu)
random.seed(3)
rnd = {}
for i in range(14):
    L = random.randint(2, 4)
    g = [random.randint(2, 5) for _ in range(L - 1)] + [random.randint(5, 12)]
    rnd[f"随机{i}"] = (g, 0)
rows = []
for nm, (g, mu) in list(real.items()) + list(rnd.items()):
    a = sum(g[mu:]) / len(g[mu:])
    res = {}
    for K in [1, 2, 3, 4]:
        t = time.time(); res[K] = nstar(g, mu, K, 400_000)
    n1 = res[1]
    line = dict(name=nm, g=g, abar=round(a, 3), **{f"N{K}": res[K] for K in res},
                **{f"LB{K}": math.ceil(K * a - 1e-9) for K in res})
    rows.append(line)
    print(json.dumps(line, ensure_ascii=False), flush=True)
json.dump(rows, open("d2_results.json", "w"), ensure_ascii=False, indent=1)
