"""随机间隔循环的图谱：稳态 + 允许等待（网络 ILP），固定模式（循环覆盖 ILP），零余量时的 K 重铺砌存在性；
K=1 时用穷尽 DFS（任意相位、等待为 0 的初态）交叉验证网络 ILP。"""
import math, random, json, sys, itertools, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import milp, LinearConstraint, Bounds
from net import solve

def cover_ilp(g, K, exact):
    S = sum(g); A = [0]
    for x in g[:-1]: A.append(A[-1] + x)
    M = np.zeros((S, S))
    for s in range(S):
        for a in A: M[s, (s - a) % S] += 1
    res = milp(c=np.ones(S), constraints=LinearConstraint(M, lb=K, ub=K if exact else np.inf),
               integrality=np.ones(S), bounds=Bounds(0, np.inf), options={"time_limit": 30})
    return None if res.x is None else int(round(res.fun))

def n_wait(g, K, Pmax):
    S, L = sum(g), len(g); lb = math.ceil(K * S / L - 1e-9); q = L // math.gcd(K, L); best = None; bp = None
    for P in range(q, Pmax + 1, q):
        v, x = solve(g, K, P, time_limit=20)
        if v is None: continue
        N = int(round(v))
        if best is None or N < best: best, bp = N, P
        if best == lb: break
    return best, bp

def dfs_k1(g, N):
    from slotmodel import Robot
    from conj2 import exact_from
    R = [Robot(g, 0, 0) for _ in range(N)]; und = False
    for phs in itertools.combinations_with_replacement(range(len(g)), N):
        r = exact_from(R, tuple((0, p) for p in phs), 1, 20000)
        if r: return True
        if r is None: und = True
    return None if und else False

def work(g):
    S, L = sum(g), len(g); a = S / L; Pmax = max(2 * S, 24); rows = []
    for K in range(1, 2 * L + 1):
        lb = math.ceil(K * a - 1e-9); nw, P = n_wait(g, K, Pmax)
        nf = cover_ilp(g, K, False)
        zero = abs(K * a - round(K * a)) < 1e-9
        tile = (cover_ilp(g, K, True) is not None) if zero else None
        row = dict(K=K, LB=lb, N_wait=nw, P=P, N_fix=nf, zero_slack=zero, tiling=tile, slack=round(lb - K * a, 4))
        if K == 1 and nw is not None and nw > lb and lb <= 7:
            row["dfs_at_LB"] = dfs_k1(g, lb)       # 期望 False（证明 LB 不可行）
        rows.append(row)
    return dict(g=g, S=S, L=L, abar=a, rows=rows)

if __name__ == "__main__":
    random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 11)
    n = int(sys.argv[1]); seqs = set()
    while len(seqs) < n:
        L = random.randint(2, 4)
        g = tuple(random.randint(2, 8) for _ in range(L))
        # 只保留循环意义下的代表（旋转等价）
        rep = min(tuple(g[i:] + g[:i]) for i in range(L))
        if sum(rep) <= 30: seqs.add(rep)
    seqs = [list(s) for s in sorted(seqs)]
    with Pool(2) as p:
        out = []
        for r in p.imap_unordered(work, seqs):
            out.append(r)
            print(json.dumps(dict(g=r["g"], gaps=[x["N_wait"] - x["LB"] if x["N_wait"] else None for x in r["rows"]],
                                  fixgaps=[x["N_fix"] - x["LB"] for x in r["rows"]],
                                  tiling=[x["tiling"] for x in r["rows"]], dfs=[x.get("dfs_at_LB") for x in r["rows"]][0]),
                             ensure_ascii=False), flush=True)
    json.dump(out, open(f"atlas_random_{n}.json", "w"), ensure_ascii=False, indent=1)
