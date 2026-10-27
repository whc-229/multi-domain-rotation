"""多超周期的"剩余均衡拉伸"上界（引理 6′）：
把 g 重复 j 次（jL 个间隔），每个间隔可拉长 e_i ≥ 0，要求总长 M ≡ 0 (mod L)，且 jL 个部分和模 L 每个剩余类恰好出现 j 次。
则 M/L 台机器人（偏移 0, L, 2L, …）给出 j 重恰好覆盖，所以 c(j) ≤ (E_min − W(j))/L，其中 W(j) = (−jS) mod L。
E_min 只依赖 g mod L 的循环序列。动态规划求 E_min，遍历全部余数模式，输出每个 (L, j) 的最坏上界。"""
import itertools, math, functools, json, sys
def emin(res, j, L):
    n = j * L
    @functools.lru_cache(None)
    def dp(i, r, counts):
        # 已确定 s''_0..s''_i（第 i 个部分和余数为 r），counts 为各余数出现次数
        if i == n:   # s''_n = M 必须 ≡ 0；此时 counts 应全为 j（s''_0..s''_{n-1}）
            return 0 if (r == 0 and all(c == j for c in counts)) else math.inf
        best = math.inf
        gi = res[i % L]
        for target in range(L):
            e = (target - r - gi) % L
            if i + 1 < n:
                if counts[target] >= j: continue
                nc = list(counts); nc[target] += 1; nc = tuple(nc)
            else:
                if target != 0: continue
                nc = counts
            v = dp(i + 1, target, nc)
            if e + v < best: best = e + v
        return best
    c0 = [0] * L; c0[0] = 1
    return dp(0, 0, tuple(c0))
if __name__ == "__main__":
  Lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
  for L in range(2, Lmax + 1):
      for j in range(1, L):
          worst = -1; wex = None; dist = {}
          seen = set()
          for res in itertools.product(range(L), repeat=L):
              key = min(tuple(res[i:] + res[:i]) for i in range(L))
              if key in seen: continue
              seen.add(key)
              S = sum(res); W = (-j * S) % L
              E = emin(res, j, L)
              b = (E - W) // L
              dist[b] = dist.get(b, 0) + 1
              if b > worst: worst, wex = b, res
          print(f"L={L} j={j}: 上界 c(j) ≤ (E_min−W)/L 的分布 {dict(sorted(dist.items()))}，最坏 {worst}，例 {wex}", flush=True)
