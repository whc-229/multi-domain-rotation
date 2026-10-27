"""K≥2 且正余量出现差距的情形：(a) 放大周期上限重解网络 ILP；(b) 穷尽 DFS（任意相位、等待为 0 的初态）确认下界不可行。"""
import itertools, math, json
from net import solve
from slotmodel import Robot
from conj2 import exact_from
cases = [([3,7,6],2),([3,7,6],5),([2,6,3],4),([3,6,8],4),([2,8,5,4],5),([2,2,3,3],2),([2,6,8,6],2)]
for g, K in cases:
    S, L = sum(g), len(g); lb = math.ceil(K*S/L - 1e-9); q = L // math.gcd(K, L)
    best = None
    for P in range(q, 4*S + 1, q):
        v, x = solve(g, K, P, time_limit=30)
        if v is not None:
            N = int(round(v)); best = N if best is None else min(best, N)
            if best == lb: break
    res = "跳过"
    if lb <= 13:
        R = [Robot(g, 0, 0) for _ in range(lb)]; res = False
        for phs in itertools.combinations_with_replacement(range(L), lb):
            r = exact_from(R, tuple((0, p) for p in phs), K, 200000)
            if r: res = True; break
            if r is None: res = "未决"
    print(json.dumps(dict(g=g, K=K, LB=lb, ILP_P到4S=best, DFS_LB可行=res), ensure_ascii=False), flush=True)
