import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import sys, time
from roster_sim import INSTANCES
from mixed_exact import Model, finite_feasible
from cfg import configs
name = sys.argv[1]; inst = INSTANCES[name](); M = Model(inst)
unk = []
for nr, nt in configs(inst, 12, rate_prune=False):
    f, n = finite_feasible(M, nr, nt, inst.horizon, 100000)
    if f is None: unk.append((nr, nt))
print(name, "未决配置", unk, flush=True)
for nr, nt in unk:
    t = time.time(); print(nr, nt, finite_feasible(M, nr, nt, inst.horizon, 20_000_000), f"{time.time()-t:.0f}s", flush=True)
