import os, sys; _h=os.path.dirname(os.path.abspath(__file__)); sys.path[:0]=[os.path.join(_h,'..','启发式验证'), os.path.join(_h,'..')]
import math, json
from slotmodel import Robot
import kpost
for g in [[2,9],[3,5,2,6],[4,4,5],[2,2,5]]:
    a = sum(g)/len(g); row = dict(g=g, abar=round(a,3))
    for K in range(1, 6):
        lb = math.ceil(K*a-1e-9); val = None
        for N in range(lb, lb+4):
            ok, _ = kpost.exact([Robot(g,0,0) for _ in range(N)], K, 600_000)
            if ok is None: val = ("?", N); break
            if ok: val = N; break
        row[f"K{K}"] = (val, lb)
        print(json.dumps(dict(g=g, K=K, N=val, LB=lb, 专点=None), ensure_ascii=False), flush=True)
