"""把图谱和定向搜索中"正余量、K=1、c=1"的数值结论用状态图精确判定（任意周期）严格化；状态数过大的跳过。"""
import json, math
from fractions import Fraction as F
from statecycle import feasible
cases = set()
for it in json.load(open("atlas_random_80.json")):
    g = it["g"]; S = sum(g); L = len(g)
    for r in it["rows"]:
        lb = math.ceil(F(r["K"] * S, L))
        if r["K"] == 1 and r["N_wait"] == lb + 1 and (lb * L - S) > 0: cases.add(tuple(g))
for l in open("c1_search.log"):
    if l.startswith("{"):
        d = json.loads(l)
        if d["K"] == 1 and d["c"] == 1: cases.add(tuple(d["g"]))
for key, rs in json.load(open("net_wait_K12.json")).items():
    g = json.loads(key); S = sum(g); L = len(g); r = rs[0]
    if r["gap"] == 1 and (math.ceil(F(S, L)) * L - S) > 0: cases.add(tuple(g))
done = 0; conf = 0; skip = 0
for g in sorted(cases, key=lambda x: sum(x)):
    g = list(g); S = sum(g); L = len(g); N = math.ceil(F(S, L))
    ns = math.comb(S + N - 1, N)
    if ns > 1_500_000:
        skip += 1; print(json.dumps(dict(g=g, skipped=True, states=ns)), flush=True); continue
    ok, _ = feasible(g, 1, N); done += 1; conf += (not ok)
    print(json.dumps(dict(g=g, N_lb=N, lb_feasible=ok, W=N * L - S, states=ns)), flush=True)
print(f"精确判定 {done} 例，其中下界不可行（c=1 得到严格确认）{conf} 例；状态过多跳过 {skip} 例")
