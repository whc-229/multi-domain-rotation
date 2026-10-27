"""班次槽模型：每台机器人 j 有自己的槽间隔序列 g_j（由其维修参数算出），值班点每槽需 1 台。
状态 = 每台 (还需等待槽数 w, 超周期相位 ph)。新机器人 = (0, 0)。"""
import math, copy, random
from roster_sim import INSTANCES
from lb import robot_occ_seq

class Robot:
    def __init__(self, g, mu, group):
        self.g, self.mu, self.L, self.group = list(g), mu, len(g), group
    def nph(self, ph):
        return ph + 1 if ph + 1 < self.L else self.mu

def gaps_from(inst, r, maint=None):
    inst = copy.deepcopy(inst)
    if maint is not None:
        inst.robots[r]["maint"] = maint
    seq, mu = robot_occ_seq(inst, r); d = inst.robots[r]["duty"]
    return [math.ceil(o / d - 1e-9) for o in seq], mu

def rate(rb):   # 稳态每槽贡献 = 周期内班数 / 周期内占用槽数
    per = rb.g[rb.mu:]
    return len(per) / sum(per)

def rate_lb(pool):
    """按速率从高到低取，最少几台才能使 sum(rate) >= 1（对任意子集都成立的下界）"""
    rs = sorted((rate(r) for r in pool), reverse=True); s = 0
    for k, v in enumerate(rs, 1):
        s += v
        if s >= 1 - 1e-12: return k
    return None

def homogeneous(name, r, N):
    g, mu = gaps_from(INSTANCES[name](), r)
    return [Robot(g, mu, 0) for _ in range(N)]

def hetero_pool(n=12, seed=7):
    """以 L3 的 0 型机器人为母本，各级维修'可用时长/维修时长'各自 ±30% 随机扰动（值岗时长不变）"""
    rng = random.Random(seed); base = INSTANCES["L3-720"]()
    pool = []
    for j in range(n):
        maint = [(round(a * rng.uniform(0.7, 1.3), 1), round(b * rng.uniform(0.7, 1.3), 1))
                 for a, b in base.robots[0]["maint"]]
        maint.sort()  # 级别越高可用时长越长
        g, mu = gaps_from(base, 0, maint)
        pool.append(Robot(g, mu, j))
    return pool

# ---------- 公共仿真
def start_state(robots):
    return [(0, 0)] * len(robots)

def step(robots, st, j):
    new = [(max(w - 1, 0), ph) for (w, ph) in st]
    rb = robots[j]; ph = st[j][1]
    new[j] = (rb.g[ph] - 1, rb.nph(ph))
    return new

def canon(robots, st):
    """同组（参数完全相同）机器人可互换 → 组内排序作为规范状态"""
    by = {}
    for rb, s in zip(robots, st):
        by.setdefault(rb.group, []).append(s)
    return tuple((k, tuple(sorted(v))) for k, v in sorted(by.items()))
