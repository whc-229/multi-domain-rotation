import math, random
import numpy as np
from slotmodel import start_state, step, canon

MAXSTEP = 20000

# ================= 1. 基线：全时域随机键（与现有 MOEA 的编码思想相同），只能证明有限时域
def baseline_fullkey(robots, budget, rng, H=144, pop=30):
    def ev(x):
        st = start_state(robots)
        for k in range(H):
            av = [j for j, (w, _) in enumerate(st) if w == 0]
            if not av: return k
            st = step(robots, st, av[min(int(x[k] * len(av)), len(av) - 1)])
        return H
    return _real_ga(ev, H, budget, rng, pop, target=H)

# ================= 2. 策略编码：每个 (机器人组, 相位) 一个优先级；确定性仿真到状态重复 => 周期证书
def policy_ga(robots, budget, rng, pop=30):
    keys = sorted({(rb.group, ph) for rb in robots for ph in range(rb.L)})
    idx = {k: i for i, k in enumerate(keys)}
    def ev(x):
        st = start_state(robots); seen = {}
        for t in range(MAXSTEP):
            c = canon(robots, st)
            if c in seen: return MAXSTEP          # 进入周期：可无限执行
            seen[c] = t
            av = [j for j, (w, _) in enumerate(st) if w == 0]
            if not av: return t
            j = max(av, key=lambda j: (x[idx[(robots[j].group, st[j][1])]], -j))
            st = step(robots, st, j)
        return t
    return _real_ga(ev, len(keys), budget, rng, pop, target=MAXSTEP)

def _real_ga(ev, n, budget, rng, pop, target):
    X = rng.random((pop, n)); F = np.array([ev(x) for x in X]); used = pop
    best = F.max()
    while used < budget and best < target:
        kids = []
        for _ in range(pop):
            a, b = [max(rng.integers(pop, size=2), key=lambda i: F[i]) for _ in range(2)]
            m = rng.random(n) < 0.5
            c = np.where(m, X[a], X[b])
            mut = rng.random(n) < 1.0 / n * 2
            c = np.clip(c + mut * rng.normal(0, 0.2, n), 0, 1 - 1e-9)
            kids.append(c)
        K = np.array(kids); FK = np.array([ev(c) for c in K]); used += pop
        allX = np.vstack([X, K]); allF = np.concatenate([F, FK])
        o = np.argsort(-allF)[:pop]; X, F = allX[o], allF[o]; best = F.max()
    return best >= target, used

# ================= 3. 状态类别 ACO：逐槽构造，信息素挂在 (组, 相位) 上；构造到状态重复即成功
def aco(robots, budget, rng, ants=20, rho=0.1, alpha=1.0):
    keys = sorted({(rb.group, ph) for rb in robots for ph in range(rb.L)})
    tau = {k: 1.0 for k in keys}; used = 0
    while used < budget:
        walks = []
        for _ in range(ants):
            st = start_state(robots); seen = set(); path = []; ok = False
            for t in range(MAXSTEP):
                c = canon(robots, st)
                if c in seen: ok = True; break
                seen.add(c)
                av = [j for j, (w, _) in enumerate(st) if w == 0]
                if not av: break
                # 同类候选合并：按 (组, 相位) 去重
                cls = {}
                for j in av: cls.setdefault((robots[j].group, st[j][1]), j)
                ks = list(cls); w = np.array([tau[k] ** alpha for k in ks]); w /= w.sum()
                k = ks[rng.choice(len(ks), p=w)]
                path.append(k); st = step(robots, st, cls[k])
            used += 1
            if ok: return True, used
            walks.append((len(path), path))
        for k in tau: tau[k] *= (1 - rho)
        walks.sort(key=lambda z: -z[0])
        for L, path in walks[:3]:
            for k in path: tau[k] += 1.0 / len(path) * L / 50
    return False, used

# ================= 4. JSP 式循环工序编码：周期 = 带重复的机器人编号序列（机器人≈工件，每一班≈工序）
def cyclic_jsp(robots, budget, rng, pop=30, ls_tries=20):
    n = len(robots)
    Lp = [rb.L - rb.mu for rb in robots]; Sp = [sum(rb.g[rb.mu:]) for rb in robots]
    # 选每台在一个周期内出现的超周期数 m_j：满足 m_j*S_j <= P = sum m_k*Lp_k，P 尽量小
    m = [1] * n
    for _ in range(50):
        P = sum(mm * l for mm, l in zip(m, Lp))
        bad = [j for j in range(n) if m[j] * Sp[j] > P]
        if not bad: break
        # 增加"太闲"的机器人（S/L 最小）的出现次数，以拉长周期
        j = min(range(n), key=lambda j: Sp[j] / Lp[j]); m[j] += 1
    else:
        return False, 0
    base = [j for j in range(n) for _ in range(m[j] * Lp[j])]
    P = len(base)

    def viol(seq):
        pos = {}
        for i, j in enumerate(seq): pos.setdefault(j, []).append(i)
        tot = 0; worst = []
        for j, ps in pos.items():
            rb = robots[j]; per = rb.g[rb.mu:]; k = len(ps)
            gaps = [(ps[(i + 1) % k] - ps[i]) % P or P for i in range(k)]
            bestv, bestlist = None, None
            for phi in range(len(per)):
                v = 0; lst = []
                for i in range(k):
                    need = per[(phi + i) % len(per)]
                    if gaps[i] < need: v += need - gaps[i]; lst.append(ps[i])
                if bestv is None or v < bestv: bestv, bestlist = v, lst
            tot += bestv; worst += bestlist
        return tot, worst

    def pox(a, b):
        keep = set(rng.choice(n, size=max(1, n // 2), replace=False).tolist())
        child = [x if x in keep else None for x in a]
        fill = iter([x for x in b if x not in keep])
        return [x if x is not None else next(fill) for x in child]

    X = [list(rng.permutation(base)) for _ in range(pop)]
    F = [viol(x) for x in X]; used = pop
    while used < budget:
        if min(f[0] for f in F) == 0: return True, used
        kids = []
        for _ in range(pop // 2):
            a, b = [min(rng.integers(pop, size=2), key=lambda i: F[i][0]) for _ in range(2)]
            c = pox(X[a], X[b])
            if rng.random() < 0.5:      # 插入变异
                i, k = rng.integers(P, size=2); c.insert(k, c.pop(i))
            fc = viol(c); used += 1
            # 冲突引导局部搜索：把一个违约位置的班次插到别处，首次改进即接受
            for _ in range(ls_tries):
                if fc[0] == 0 or used >= budget: break
                i = rng.choice(fc[1]); k = rng.integers(P)
                d = c[:]; d.insert(k, d.pop(i)); fd = viol(d); used += 1
                if fd[0] <= fc[0]: c, fc = d, fd
            kids.append((c, fc))
        allp = list(zip(X, F)) + kids
        allp.sort(key=lambda z: z[1][0]); X = [z[0] for z in allp[:pop]]; F = [z[1] for z in allp[:pop]]
    return min(f[0] for f in F) == 0, used
