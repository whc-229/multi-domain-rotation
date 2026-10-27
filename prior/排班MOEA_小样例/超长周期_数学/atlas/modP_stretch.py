"""引理 6‴（模 P 混合）：模数 P（L | P 不必要），若干类机器人，第 i 类按拉伸模式 g^{(i)}（单超周期，总长 M_i，P | M_i）循环，
偏移 c_i + P·Z_{M_i}（M_i/P 台）。若各类部分和（加偏移）模 P 的多重集恰好把 Z_P 覆盖 K 次，则 N(K) ≤ Σ M_i/P。
这里针对 L=5 的例外类 (0,2,0,4,4)：取 P = 10、K = 1、两类，遍历它模 10 的全部 32 种提升。"""
import itertools, math
def cand(res10, P, L):
    out = {}   # 残余集合（frozenset）→ 最小拉伸
    for e in itertools.product(range(P), repeat=L):
        s = 0; R = []
        for i in range(L): R.append(s % P); s += res10[i] + e[i]
        if s % P or len(set(R)) < L: continue
        fs = frozenset(R); E = sum(e)
        if fs not in out or E < out[fs]: out[fs] = E
    return out
base = (0, 2, 0, 4, 4); L = 5; P = 10
worst = -1
for lift in itertools.product([0, 1], repeat=L):
    res10 = tuple(r + 5 * b for r, b in zip(base, lift))
    S = sum(res10); W = (-S) % L
    C = cand(res10, P, L); best = None
    items = list(C.items())
    for R1, E1 in items:
        for R2, E2 in items:
            for c in range(P):
                R2s = frozenset((x + c) % P for x in R2)
                if not (R1 & R2s) and len(R1 | R2s) == P:
                    if best is None or E1 + E2 < best: best = E1 + E2
    # c(1) ≤ (2S + E1 + E2)/10 − ⌈S/5⌉，S 用代表元 res10 计算（只依赖 S mod 10 与拉伸量）
    bound = None if best is None else (2 * S + best) // 10 - math.ceil(S / 5)
    worst = max(worst, bound if bound is not None else 99)
    print(res10, "最小 E1+E2 =", best, "→ c(1) 上界", bound, flush=True)
print("最坏上界", worst)
