"""带工具池的进化派工策略：8 维权重（前 6 维选机器人，后 2 维选工具），(μ+λ)-ES。适应度 = 机器人 + 工具最少总数。"""
import numpy as np, json
from tools_model import size_resources, simulate_tools, tool_lb

def evolve_tools(TI, H, warm, gens=20, mu=4, lam=12, seed=0, fresh=True):
    rng = np.random.default_rng(seed)
    base = [np.array([1, 0, 0, 0, 0, 0, 1, 0.]), np.array([0.1, 0, 10, 1, 0, 0, 1, 0.1]), np.array([1, 0.5, 1, 0.5, 0, 0, 1, -0.1])]
    def fit(w):
        R, Nk, Mj = size_resources(TI, "score", H, warm, params=w, fresh=fresh)
        return (R if R is not None else 999), Nk, Mj
    pop = base + [rng.normal(0, 1, 8) for _ in range(mu + lam - len(base))]
    sc = sorted([(fit(w), w) for w in pop], key=lambda z: z[0][0]); par = sc[:mu]; hist = []
    for g in range(gens):
        kids = [par[rng.integers(len(par))][1] + 0.5 * rng.normal(0, 1, 8) * (np.abs(par[rng.integers(len(par))][1]) + 0.2) for _ in range(lam)]
        par = sorted(par + [(fit(w), w) for w in kids], key=lambda z: z[0][0])[:mu]
        hist.append(par[0][0][0])
    (R, Nk, Mj), w = par[0]
    return R, Nk, Mj, w.tolist(), hist
