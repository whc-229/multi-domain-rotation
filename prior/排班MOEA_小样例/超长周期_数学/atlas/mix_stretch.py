"""混合拉伸构造（引理 6″）：K 类机器人，共同周期 M = S + E（L | M），第 i 类按某个单超周期拉伸 e^{(i)}（Σe = E）循环，
偏移为 c_i + L·Z_M（每类 M/L 台）。覆盖数 = Σ_i hist_i(x − c_i)，hist 为部分和模 L 的直方图。若和为常数 K，则 N(K) ≤ K·M/L，
c(K) ≤ (K·E − W(K))/L。与引理 6′ 的多超周期界取小，遍历全部循环余数类，输出 c 的上界分布，并单独统计 W ≥ 2 的情形（猜想 C2）。"""
import itertools, math, sys, json
from functools import lru_cache
from balanced_stretch import emin
def comps(E, L):
    for cut in itertools.combinations(range(E + L - 1), L - 1):
        prev = -1; parts = []
        for c in cut: parts.append(c - prev - 1); prev = c
        parts.append(E + L - 2 - prev); yield parts
def hists(res, E, L):
    out = set()
    for e in comps(E, L):
        s = 0; h = [0] * L
        for i in range(L): h[s % L] += 1; s += res[i] + e[i]
        assert s % L == (sum(res) + E) % L
        for c in range(L): out.add(tuple(h[(r - c) % L] for r in range(L)))
    return out
def mix_ok(res, E, K, L):
    H = list(hists(res, E, L)); target = tuple([K] * L)
    @lru_cache(None)
    def f(k, acc):
        if k == 0: return acc == target
        return any(f(k - 1, tuple(a + b for a, b in zip(acc, h))) for h in H if all(a + b <= K for a, b in zip(acc, h)))
    return f(K, tuple([0] * L))
def mix_bound(res, K, L):
    S = sum(res); W = (-K * S) % L; E0 = (-S) % L
    for E in range(E0, E0 + 3 * L, L):
        if (K * E - W) // L > 1: break
        if mix_ok(res, E, K, L): return (K * E - W) // L
    return None
if __name__ == "__main__":
    Lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for L in range(2, Lmax + 1):
        seen = set(); stats = {}; c2bad = []; c1bad = []
        for res in itertools.product(range(L), repeat=L):
            key = min(tuple(res[i:] + res[:i]) for i in range(L))
            if key in seen: continue
            seen.add(key)
            for K in range(1, L):
                W = (-K * sum(key)) % L
                b1 = (emin(key, K, L) - W) // L
                b2 = mix_bound(key, K, L) if K >= 2 else None
                b = min(x for x in [b1, b2] if x is not None)
                stats[(K, b)] = stats.get((K, b), 0) + 1
                if b >= 2: c1bad.append((key, K, b1, b2))
                if W >= 2 and b >= 1: c2bad.append((key, K, W, b1, b2))
        print(f"L={L}: (K, 上界) 计数 {dict(sorted(stats.items()))}", flush=True)
        print(f"   C1 未覆盖（上界≥2）: {c1bad[:12]}{' …' if len(c1bad) > 12 else ''} 共 {len(c1bad)}", flush=True)
        print(f"   C2 未覆盖（W≥2 但上界≥1）: {c2bad[:12]}{' …' if len(c2bad) > 12 else ''} 共 {len(c2bad)}", flush=True)
