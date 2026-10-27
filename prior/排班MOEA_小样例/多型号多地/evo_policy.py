"""进化派工策略（试点）：打分策略 score(r) = w·特征，(μ+λ)-ES 进化 w。
适应度 = size_fleet 得到的最少总台数（主），以及下界配比下的缺岗格数（次，作平滑）。策略与时域无关，评估时域 H 可任意加长。"""
import numpy as np, json, sys, time
from sim import size_fleet, simulate
from mtms import fluid_lb

def fitness(I, w, H, warm, start):
    N, Nk = size_fleet(I, "score", H, warm, params=w, start=start)
    g0, _ = simulate(I, start, "score", H, warm, params=w)
    return (N if N is not None else 999) + g0 / (10.0 * H), N, Nk

def evolve(I, H=3000, warm=1000, gens=30, mu=4, lam=12, seed=0, log=None):
    rng = np.random.default_rng(seed)
    _, x = fluid_lb(I, integer=True); start = [int(round(v)) for v in x]
    # 初始种群：含与 FIFO、相位平衡近似的权重
    base = [np.array([1, 0, 0, 0, 0, 0.]), np.array([0.1, 0, 10, 1, 0, 0]), np.array([1, 0.5, 1, 0.5, 0, 0])]
    pop = [b for b in base] + [rng.normal(0, 1, 6) for _ in range(mu - len(base) + lam)]
    scored = [(fitness(I, w, H, warm, start), w) for w in pop]
    scored.sort(key=lambda z: z[0][0]); par = scored[:mu]; sigma = 0.5; hist = []
    for g in range(gens):
        kids = []
        for _ in range(lam):
            a = par[rng.integers(len(par))][1]
            kids.append(a + sigma * rng.normal(0, 1, 6) * (np.abs(a) + 0.2))
        sk = [(fitness(I, w, H, warm, start), w) for w in kids]
        allp = sorted(par + sk, key=lambda z: z[0][0]); par = allp[:mu]
        hist.append(par[0][0][1])
        if log: print(json.dumps(dict(gen=g, best=par[0][0][1], Nk=par[0][0][2])), file=log, flush=True)
    return par[0][0][1], par[0][0][2], par[0][1].tolist(), hist
