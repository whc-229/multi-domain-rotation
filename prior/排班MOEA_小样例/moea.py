"""NSGA-II 骨架 + 三个可开关的场景化改进
  B1 graded : 不可行解按"仿真推进到哪里"分级惩罚（原做法：一律 30000）
  B2 sim    : 仿真感知变异——利用父代解码轨迹，定点修改关键班次的基因
              · 父代不可行：把失败班次附近的机器人/工具基因改到"复用已有资源"的区间
              · 父代可行：一半概率"减资源"（把新开机器人的班次改为复用），
                          一半概率"均衡"（把出勤最多的机器人的某个班次换给其他可用机器人）
  A  level  : 按资源数 f1 分层的精英档案：每个 f1 取值保留均衡度最好的解（即使被其他层支配），
              繁殖时以概率 p_re 从档案取父代
所有变体都记录"全部已评价可行解"的非支配集作为最终输出，保证公平比较。
"""
import numpy as np
import roster_sim as R


def nd_rank(F):
    n = len(F)
    le = (F[:, None, :] <= F[None, :, :]).all(2)
    lt = (F[:, None, :] < F[None, :, :]).any(2)
    D = le & lt
    cnt = D.sum(0)
    rank = np.full(n, -1)
    r, cur = 0, np.where(cnt == 0)[0]
    while len(cur):
        rank[cur] = r
        cnt = cnt - D[cur].sum(0)
        cnt[rank >= 0] = -1
        cur = np.where(cnt == 0)[0]
        r += 1
    return rank


def crowding(F):
    n, m = F.shape
    d = np.zeros(n)
    if n <= 2:
        return np.full(n, np.inf)
    for j in range(m):
        o = np.argsort(F[:, j], kind="stable")
        rng = F[o[-1], j] - F[o[0], j]
        d[o[0]] = d[o[-1]] = np.inf
        if rng > 0:
            d[o[1:-1]] += (F[o[2:], j] - F[o[:-2], j]) / rng
    return d


def select_env(F, N):
    rank = nd_rank(F)
    keep = []
    for r in range(rank.max() + 1):
        idx = np.where(rank == r)[0]
        if len(keep) + len(idx) <= N:
            keep.extend(idx)
        else:
            cd = crowding(F[idx])
            keep.extend(idx[np.argsort(-cd)[: N - len(keep)]])
            break
    keep = np.array(keep)
    rank_k = rank[keep]
    cd_k = np.zeros(len(keep))
    for r in np.unique(rank_k):
        m = rank_k == r
        cd_k[m] = crowding(F[keep[m]])
    return keep, rank_k, cd_k


def sbx_pm(p1, p2, rng, eta_c=15, eta_m=20, pc=0.9):
    n = len(p1)
    c1, c2 = p1.copy(), p2.copy()
    if rng.random() < pc:
        u = rng.random(n)
        beta = np.where(u <= .5, (2 * u) ** (1 / (eta_c + 1)), (1 / (2 * (1 - u))) ** (1 / (eta_c + 1)))
        swap = rng.random(n) < .5
        m = (rng.random(n) < .5) & (np.abs(p1 - p2) > 1e-12)
        a = .5 * ((1 + beta) * p1 + (1 - beta) * p2)
        b = .5 * ((1 - beta) * p1 + (1 + beta) * p2)
        c1 = np.where(m, np.where(swap, b, a), p1)
        c2 = np.where(m, np.where(swap, a, b), p2)
    for c in (c1, c2):
        pm = rng.random(n) < 1.0 / n
        u = rng.random(n)
        dl = np.where(u < .5, (2 * u) ** (1 / (eta_m + 1)) - 1, 1 - (2 * (1 - u)) ** (1 / (eta_m + 1)))
        c[pm] += dl[pm]
    return np.clip(c1, 0, 0.999999), np.clip(c2, 0, 0.999999)


_DEC = {"name": "orig", "theta": 0.85}


def reuse_gene(n_avail, n_types, rng):
    """返回落在"复用已有机器人"区间内的随机键（与当前解码方式一致）"""
    if _DEC["name"] == "orig":
        return rng.random() * n_avail / (n_avail + n_types)
    return rng.random() * _DEC["theta"]


def sim_mutate(x, info, inst, rng, window=5):
    """仿真感知变异：只改一两个关键班次的基因"""
    x = x.copy(); S = inst.S; nT = len(inst.robots)
    na, no, nc = info["n_avail"], info["n_oldtool"], info["n_cand"]
    if not info["feasible"]:
        k = info["fail_shift"]
        for s in range(max(0, k - window), min(k + 1, len(na))):
            if na[s] > 0 and rng.random() < .5:
                x[s] = reuse_gene(na[s], nT, rng)
            if s < len(no) and no[s] > 0 and rng.random() < .5:
                x[s + S] = rng.random() * no[s] / nc[s]
        return x
    if rng.random() < .5 and info["new_shift"]:            # 减资源：把一次"新开"改为复用
        cands = [s for s in info["new_shift"] if na[s] > 0]
        if cands:
            s = cands[rng.integers(len(cands))]
            x[s] = reuse_gene(na[s], nT, rng)
            if no[s] > 0:
                x[s + S] = rng.random() * no[s] / nc[s]
    else:                                                   # 均衡：出勤最多的机器人让出一个班次
        sr = np.array(info["shift_robot"])
        busiest = np.bincount(sr).argmax()
        shifts = np.where((sr == busiest) & (np.array(na) > 1))[0]
        if len(shifts):
            s = shifts[rng.integers(len(shifts))]
            x[s] = reuse_gene(na[s], nT, rng)
    return x


def fitness(info, graded):
    if info["feasible"]:
        return info["f"]
    if graded:
        return np.full(2, R.INFEAS + 1000 * (1 - info["progress"]))
    return np.full(2, R.INFEAS)


def update_nd(ndF, f):
    """维护所有已评价可行解的非支配集（只存目标值）"""
    if ((ndF <= f).all(1) & (ndF < f).any(1)).any() or (ndF == f).all(1).any():
        return ndF
    keep = ~((f <= ndF).all(1) & (f < ndF).any(1))
    return np.vstack([ndF[keep], f])


def run(inst, graded=False, sim=False, level=False, pop=100, gens=100, seed=0,
        p_sim=0.5, p_re=0.3, seed_greedy=False, seed_rule=False):
    rng = np.random.default_rng(seed)
    _DEC["name"], _DEC["theta"] = inst.decoder, inst.theta
    n = inst.n_var
    X = rng.random((pop, n))
    if seed_greedy:   # 场景知识初始化：一部分个体的机器人基因压到"复用"区间
        S = inst.S
        lim = (0.3 if inst.decoder == "orig" else inst.theta)
        X[: pop // 5, :S] *= lim * rng.random((pop // 5, 1))
    if seed_rule:     # 规则播种：全 0 基因 = 每班都复用（排序解码下即"最闲优先"规则），再加 9 个小扰动版本
        X[0] = 0.0
        X[1:10] = rng.random((9, n)) * 0.1
    infos = [R.evaluate(x, inst) for x in X]
    F = np.array([fitness(i, graded) for i in infos])
    ndF = np.zeros((0, 2))
    arch = {}                                   # f1 -> (f2, x, info)
    first_feas = None

    def absorb(x, info, g):
        nonlocal ndF, first_feas
        if info["feasible"]:
            ndF = update_nd(ndF, info["f"])
            if first_feas is None:
                first_feas = g
            k = int(info["f"][0])
            if k not in arch or info["f"][1] < arch[k][0]:
                arch[k] = (info["f"][1], x.copy(), info)

    for x, i in zip(X, infos):
        absorb(x, i, 0)
    _, rank, cd = select_env(F, pop)
    for g in range(1, gens + 1):
        def tour():
            a, b = rng.integers(pop, size=2)
            return a if (rank[a], -cd[a]) < (rank[b], -cd[b]) else b

        def parent():
            if level and arch and rng.random() < p_re:
                ks = list(arch)
                _, x, info = arch[ks[rng.integers(len(ks))]]
                return x, info
            i = tour()
            return X[i], infos[i]

        kids, kid_par = [], []
        while len(kids) < pop:
            (x1, i1), (x2, i2) = parent(), parent()
            c1, c2 = sbx_pm(x1, x2, rng)
            kids += [c1, c2]; kid_par += [i1, i2]
        kids = np.array(kids[:pop]); kid_par = kid_par[:pop]
        if sim:
            for j in range(pop):
                if rng.random() < p_sim:
                    kids[j] = sim_mutate(kids[j], kid_par[j], inst, rng)
        kinfo = [R.evaluate(x, inst) for x in kids]
        for x, i in zip(kids, kinfo):
            absorb(x, i, g)
        KF = np.array([fitness(i, graded) for i in kinfo])
        allX = np.vstack([X, kids]); allF = np.vstack([F, KF]); allI = infos + kinfo
        keep, rank, cd = select_env(allF, pop)
        X, F, infos = allX[keep], allF[keep], [allI[k] for k in keep]
    return dict(ndF=ndF, first_feas=first_feas,
                level_best={k: v[0] for k, v in arch.items()})


if __name__ == "__main__":
    import time
    inst = R.INSTANCES["L3-720"]()
    for kw in [dict(), dict(graded=True, sim=True, level=True)]:
        t = time.time(); r = run(inst, gens=20, **kw)
        print(kw, "非支配解数", len(r["ndF"]), "最小f1", r["ndF"][:, 0].min(), "用时", round(time.time() - t, 1))
