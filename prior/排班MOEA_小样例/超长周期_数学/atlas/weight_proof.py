"""g=[1,2k]（S=2k+1）全新启动的显式权重证书：
w(t)=⌊t/2⌋+1 (0≤t<2k)，w(2k)=k+1，w(t)=⌈(4k−t)/2⌉ (2k<t<4k)，其余 0。
断言：任何一条从相位 0 出发的单机轨迹（可等待）收集的 Σw ≤ S；且 Σ_t w(t) = 2k²+2k+1。
于是每格恰好 K 台时 N·S ≥ K·Σw，即 N_fresh(K) ≥ K(S²+1)/(2S)。"""
from functools import lru_cache
def maxpath(g, w):
    L = len(g); H = len(w)
    @lru_cache(None)
    def best(t, f):            # 在 t 时刻（或之后）可以以相位 f 上岗的最大收集量
        if t >= H: return 0
        return max(best(t + 1, f), w[t] + best(t + g[f], (f + 1) % L))
    return best(0, 0)
for k in range(1, 41):
    S = 2 * k + 1; g = [1, 2 * k]
    w = [(t // 2 + 1) if t < 2 * k else (k + 1 if t == 2 * k else -(-(4 * k - t) // 2)) for t in range(4 * k)]
    mp = maxpath(g, w); tot = sum(w)
    assert tot == 2 * k * k + 2 * k + 1, (k, tot)
    assert mp <= S, (k, mp)
print("k=1..40 全部通过：单机最大收集量 ≤ S，总权重 = 2k²+2k+1")
