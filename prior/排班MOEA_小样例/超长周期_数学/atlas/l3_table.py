"""L=3 定理 C1 的证书表：对每个循环余数类与 j=1,2，给出一组最优拉伸向量 e（长度 3j）及上界 (E−W)/3。"""
import itertools, math
def best_stretch(res, j, L):
    n = j * L; best = None
    for e in itertools.product(range(L), repeat=n):
        s = 0; sums = []
        for i in range(n):
            sums.append(s % L); s += res[i % L] + e[i]
        if s % L: continue
        if all(sums.count(r) == j for r in range(L)):
            E = sum(e)
            if best is None or E < best[0]: best = (E, e)
    return best
L = 3; seen = set()
print("余数类(g mod 3) | j | W(j) | E_min | 上界 c(j) | 拉伸向量 e")
for res in itertools.product(range(L), repeat=L):
    key = min(tuple(res[i:] + res[:i]) for i in range(L))
    if key in seen: continue
    seen.add(key)
    for j in [1, 2]:
        W = (-j * sum(key)) % L; E, e = best_stretch(key, j, L)
        print(key, "|", j, "|", W, "|", E, "|", (E - W) // L, "|", e)
