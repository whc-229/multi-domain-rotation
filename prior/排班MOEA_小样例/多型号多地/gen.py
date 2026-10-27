"""MT-MS 实例生成：机型取自真实 9 型，站点随机限制机型（每站 2–3 种）、随机额外路程。"""
import random, json
from mtms import real_types, Inst

def make(m, types, seed, kmin=2, kmax=3, deltas=(0, 1, 2, 4), K=1, cost=None):
    rng = random.Random(seed); T = real_types()
    sel = {t: T[t] for t in types}
    sites = []
    for s in range(m):
        E = rng.sample(types, rng.randint(kmin, min(kmax, len(types))))
        sites.append(dict(E=sorted(E), delta=rng.choice(deltas), K=K))
    used = sorted({t for s in sites for t in s["E"]})
    return Inst({t: sel[t] for t in used}, sites, cost={t: (cost or {}).get(t, 1.0) for t in used}, name=f"m{m}-n{len(types)}-s{seed}")

SIX = ["L1-0", "L1-1", "L1-2", "L3-0", "L3-1", "L3-2"]
NINE = ["L1-0", "L1-1", "L1-2", "L2-0", "L2-1", "L2-2", "L3-0", "L3-1", "L3-2"]
