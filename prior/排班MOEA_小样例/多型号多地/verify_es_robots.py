"""复核 scale.py 里只算机器人的派工 ES 解：原来只在 2 万格上查过无缺岗，这里拉长到 30 万格（约 68 年）。"""
import json, sys, time
from gen import make, SIX
from sim import simulate
for fn in ["scale.log", "scale_m10.log"]:
    for l in open(fn):
        d = json.loads(l)
        if "es_w" not in d: continue
        m, sd = d["inst"].split("-")[0][1:], d["inst"].split("-s")[1]
        I = make(int(m), SIX, int(sd)); t0 = time.time()
        g, gs = simulate(I, d["es_Nk"], "score", 300000, 1000, params=d["es_w"])
        print(json.dumps(dict(inst=d["inst"], es=d["es"], gaps_300k=g, gap_sites=gs, sec=round(time.time() - t0))), flush=True)
