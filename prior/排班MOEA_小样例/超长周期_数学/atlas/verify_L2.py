"""定理（L=2，稳态，允许等待）数值核对：
N_w(K) = ⌈K(a+b)/2⌉ + 1，若 a、b 都是偶数、2-adic 赋值不同、且 K 为奇数；否则 = ⌈K(a+b)/2⌉。"""
import math, json, sys
from net import solve
def v2(x):
    k = 0
    while x % 2 == 0: x //= 2; k += 1
    return k
def pred(a, b, K):
    base = math.ceil(K * (a + b) / 2)
    return base + (1 if (a % 2 == 0 and b % 2 == 0 and v2(a) != v2(b) and K % 2 == 1) else 0)
bad = []; n = 0
amax = int(sys.argv[1]) if len(sys.argv) > 1 else 10
for a in range(1, amax + 1):
    for b in range(a, amax + 1):
        g = [a, b]; S = a + b
        for K in range(1, 5):
            lb = math.ceil(K * S / 2); q = 2 // math.gcd(K, 2); best = None
            for P in range(q, max(2 * S, 24) + 1, q):
                v, x = solve(g, K, P, time_limit=20)
                if v is None: continue
                N = int(round(v)); best = N if best is None else min(best, N)
                if best == lb: break
            n += 1
            if best != pred(a, b, K): bad.append((a, b, K, best, pred(a, b, K)))
print("核对", n, "组；不一致", len(bad), bad[:20])
