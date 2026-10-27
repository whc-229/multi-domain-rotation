"""同型机器人、K 个值班点：精确最少台数 N* 与速率下界 ⌈K·ā⌉ 的关系（ā = 周期内平均槽间隔）"""
import math, random, itertools, sys
from slotmodel import Robot
import kpost
random.seed(1)
rows = []
for trial in range(int(sys.argv[1])):
    L = random.randint(1, 4)
    g = [random.randint(2, 6) for _ in range(L - 1)] + [random.randint(4, 16)]
    K = int(sys.argv[2])
    abar = sum(g) / len(g); lb = math.ceil(K * abar - 1e-9)
    Nstar = None
    for N in range(lb, lb + 4):
        R = [Robot(g, 0, 0) for _ in range(N)]
        ok, nodes = kpost.exact(R, K, 150_000)
        if ok is None: Nstar = ("?", N); break
        if ok: Nstar = N; break
    rows.append((g, K, abar, lb, Nstar))
from collections import Counter
c = Counter(); bad = []
for g, K, a, lb, N in rows:
    if isinstance(N, tuple) or N is None: c["未决"] += 1; continue
    c[(K, N - lb)] += 1
    if N - lb >= 1: bad.append((g, K, round(a, 2), lb, N))
print(sorted(c.items(), key=str))
print("N* > 下界 的例子（前 40）:")
for b in bad[:40]: print(b)
# 检验猜想：N*=lb 当且仅当 lb - K·ā ≥ 某阈值？打印 slack 分布
import statistics
for k in [1,2,3]:
    s0 = sorted(round(lb - K*a, 3) for g, K, a, lb, N in rows if K == k and N == lb)
    s1 = sorted(round(lb - K*a, 3) for g, K, a, lb, N in rows if K == k and isinstance(N, int) and N > lb)
    print(f"K={k}: 达到下界的余量(lb-Kā) min={s0[:1]} | 未达到的余量 max={s1[-1:] if s1 else None} 个数 {len(s0)}/{len(s1)}")
