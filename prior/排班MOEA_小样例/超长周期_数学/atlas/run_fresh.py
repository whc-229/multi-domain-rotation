"""T3：全新启动的有限时域最少台数 N_H（时间展开 ILP），H = 1…8S 的若干点；与稳态（网络 ILP）比较启动代价；
reset=True 时允许提前大修。"""
import json, math, sys
from fresh import fresh_min
seqs = [[3,10],[2,9],[4,4,5],[3,5,2,6],[2,2,5],[3,3,15],[4,4,5,4,9],[3,6,3,6,9],[4,6,5,4],[2,6,3]]
reset = len(sys.argv) > 1 and sys.argv[1] == "reset"
Kmax = int(sys.argv[2]) if len(sys.argv) > 2 else 10
steady = json.load(open("net_wait_K12.json"))
for g in seqs:
    S = sum(g); a = S / len(g)
    Hs = sorted(set([S // 2, S, 2 * S, 3 * S, 4 * S, 6 * S, 8 * S]))
    rows = []
    for K in range(1, Kmax + 1):
        NH = [fresh_min(g, K, H, reset=reset, time_limit=60) for H in Hs]
        st = None
        key = str(g)
        if key in steady and K <= len(steady[key]): st = steady[key][K - 1]["N"]
        # turnpike：N_H 首次达到 H=8S 时的值的 H
        lim = NH[-1]; H0 = next(H for H, n in zip(Hs, NH) if n == lim)
        rows.append(dict(K=K, LB=math.ceil(K * a - 1e-9), N_steady=st, N_fresh_8S=lim, startup=(lim - st) if st is not None and lim is not None else None, H0_upper=H0, NH=NH))
    print(json.dumps(dict(g=g, reset=reset, Hs=Hs, startup=[r["startup"] for r in rows], N_fresh=[r["N_fresh_8S"] for r in rows],
                          N_steady=[r["N_steady"] for r in rows], H0=[r["H0_upper"] for r in rows]), ensure_ascii=False), flush=True)
