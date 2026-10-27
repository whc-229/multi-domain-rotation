"""猜想 C4'：全新启动流体值 α_fresh = ā ⇔ A(z) 整除 1+z+…+z^{S-1} 且商的系数非负（即 A 的平移可以不越界地恰好铺满区间 [0,S)，允许分数重数）。
对随机与穷举的小间隔序列核对（α 用 H=3S 的 LP）。"""
import itertools, random, json
import numpy as np
from fractions import Fraction as F
from fresh_lp3 import fresh_lp
def interval_tile(g):
    S = sum(g); L = len(g); A = [sum(g[:i]) for i in range(L)]
    # 三角方程 Σ_a B(x-a) = 1, x∈[0,S)，唯一解；检查非负且不越界
    B = [F(0)] * S
    for x in range(S):
        B[x] = 1 - sum(B[x - a] for a in A if 0 < a <= x)
    ok_nonneg = all(b >= 0 for b in B)
    ok_noover = all(B[b] == 0 for b in range(S) if b + A[-1] >= S)
    return ok_nonneg and ok_noover
def canon(g): return min(tuple(g[i:] + g[:i]) for i in range(len(g)))  # 注意：全新启动不是循环不变的，这里不去重
rows = []; agree = 0; n = 0; bad = []; pos = 0; posex = []
cands = []
for L in [2, 3, 4]:
    for g in itertools.product(range(1, 7), repeat=L):
        cands.append(list(g))
random.seed(3)
cands = random.sample(cands, 400) + [[c] * (L - 1) + [G] for L in [3, 4] for c in [1, 2] for G in range(c, 4 * L * c)]
for g in cands:
    S = sum(g); abar = F(S, len(g))
    a = F(fresh_lp(g, 3 * S)).limit_denominator(5000)
    lhs = (a == abar); rhs = interval_tile(g)
    n += 1; agree += (lhs == rhs); pos += lhs
    if lhs: posex.append(g)
    if lhs != rhs: bad.append((g, str(a), str(abar), rhs))
print("核对", n, "一致", agree, "其中 α=ā 的实例", pos); print("不一致例:", bad[:30]); print("正例样本:", posex[:25])
