import json, math, sys, time
from net import solve, build, decompose
seqs = [[3,10],[2,9],[4,4,5],[3,5,2,6],[2,2,5],[3,3,15],[4,4,5,4,9],[3,6,3,6,9],[5,5,5,21],[3,3,15,4],[4,6,5,4]]
Kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
mode = sys.argv[2] if len(sys.argv) > 2 else "wait"     # wait | dummy | reset
res = {}
for g in seqs:
    S, L = sum(g), len(g); a = S / L; row = []
    Pmax = max(2 * S, 40)
    for K in range(1, Kmax + 1):
        lb = math.ceil(K * a - 1e-9); best = None; bestP = None; lpmin = None
        q = L // math.gcd(K, L) if mode == "wait" else 1   # 精确覆盖时必须 L | K·P（每个圈的工作数是 L 的倍数）
        for P in range(q, Pmax + 1, q):
            v, x = solve(g, K, P, exact=(mode != "dummy"), reset=(mode == "reset"), time_limit=20)
            if v is None: continue
            N = int(round(v))
            if best is None or N < best:
                best, bestP = N, P
                arcs, _, _ = build(g, K, P, exact=(mode != "dummy"), reset=(mode == "reset"))
                rob, ok = decompose(g, P, x, arcs); assert ok and rob == N, (g, K, P, rob, N)
            if best == lb: break
        lv, _ = solve(g, K, bestP, exact=(mode != "dummy"), reset=(mode == "reset"), integer=False)
        row.append(dict(K=K, N=best, LB=lb, gap=best - lb, P=bestP, LP_at_P=round(lv, 4)))
    res[str(g)] = row
    print(json.dumps(dict(g=g, abar=round(a,3), gaps=[r["gap"] for r in row], N=[r["N"] for r in row], P=[r["P"] for r in row]), ensure_ascii=False), flush=True)
json.dump(res, open(f"net_{mode}_K{Kmax}.json", "w"), ensure_ascii=False, indent=1)
