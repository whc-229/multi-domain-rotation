"""复核 [4,4,5] 大 K 时 √K 是否成立：更长时域 T=60000、3 个随机种子，报告各种子的最少台数。"""
import json, math
from stoch2 import min_fleet
for gs, G0, d in [([4, 4], 5, 2), ([3], 10, 2)]:
    for K in [16, 32, 64, 128]:
        Ns = [min_fleet(gs, G0, d, K, 1e-3, 60000, "PB", seed=s)[0] for s in range(3)]
        a = (sum(gs) + G0) / (len(gs) + 1)
        print(json.dumps(dict(g=gs + [G0], spread=d, K=K, N_seeds=Ns, excess=[round(n - K * a, 2) for n in Ns], ratio=[round((n - K * a) / math.sqrt(K), 2) for n in Ns])), flush=True)
