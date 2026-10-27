"""周期结构编码（试点，只算机器人、稳态）：
  决策 = 把站点分组，每组指定一个对组内所有站点都允许的机型；同组机器人只用这一型，按上一轮的显式构造（剩余均衡拉伸 / 混合拉伸）轮流值守组内站点。
  每组台数由公式给出：⌈K·ā⌉ ≤ N ≤ ⌈K·ā⌉ + b，b 来自构造（只依赖间隔 mod L 与 K mod L），排班是显式周期表，天然无限期可行，不需要仿真。
  间隔按"槽"计：槽长 = 该型值岗格数 d；组内跨站点时取最保守的间隔 ⌈(nx_max − lead_min + lead_max)/d⌉。
外层：站点数 ≤ 12 时对"集合划分"做精确动态规划（O(3^m)），给出这一编码下的最优值。"""
import sys, os, math, json, functools, itertools
_here = os.path.dirname(os.path.abspath(__file__))
for _p in ("../atlas", "../超长周期_数学/atlas"):   # 服务器上 / 电脑上 atlas 目录的位置
    sys.path.insert(0, os.path.join(_here, _p))
from fractions import Fraction as F
from balanced_stretch import emin
from mix_stretch import mix_bound
from gen import make, SIX
from mtms import fluid_lb

@functools.lru_cache(None)
def c_bound(res, k, L):
    """剩余类构造给出的 c(K) 上界（K ≡ k mod L）；k = 0 时为 0。"""
    if k == 0: return 0
    W = (-k * sum(res)) % L
    b1 = (emin(res, k, L) - W) // L
    b2 = mix_bound(res, k, L) if k >= 2 else None
    return min(x for x in [b1, b2] if x is not None)

def group_cost(I, k, G):
    T = I.T[k]; d = I.d[k]; L = I.L[k]
    leads = [I.lead(s) for s in G]
    gaps = []
    for i in range(L):
        nxm = max(I.step(k, i, s)[1] for s in G)
        gaps.append(math.ceil((nxm - min(leads) + max(leads)) / d))
    K = len(G); S = sum(gaps)
    lb = math.ceil(F(K * S, L))
    res = tuple(g % L for g in gaps); res = min(tuple(res[i:] + res[:i]) for i in range(L))
    b = c_bound(res, K % L, L)
    return lb + b, lb, gaps

def solve(I):
    m = len(I.sites); full = (1 << m) - 1
    best_single = {}
    for mask in range(1, full + 1):
        G = [s for s in range(m) if mask >> s & 1]
        opts = []
        for k in range(len(I.T)):
            if all(s in I.elig[k] for s in G):
                ub, lb, gaps = group_cost(I, k, G); opts.append((ub, lb, I.tn[k], gaps))
        if opts: best_single[mask] = min(opts)
    @functools.lru_cache(None)
    def dp(mask):
        if mask == 0: return (0, ())
        low = mask & -mask; best = (10**9, ())
        sub = mask
        while sub:
            if sub & low and sub in best_single:
                r = dp(mask ^ sub)
                v = best_single[sub][0] + r[0]
                if v < best[0]: best = (v, ((sub, best_single[sub]),) + r[1])
            sub = (sub - 1) & mask
        return best
    v, groups = dp(full)
    return v, [([s for s in range(m) if g >> s & 1], c[2], c[0], c[1]) for g, c in groups]

if __name__ == "__main__":
    for spec in sys.argv[1].split(","):
        m, sd = spec.split("s"); I = make(int(m), SIX, int(sd))
        v, groups = solve(I)
        print(json.dumps(dict(inst=I.name, struct_ub=v, intLB=round(fluid_lb(I, True)[0]), groups=groups), ensure_ascii=False), flush=True)
