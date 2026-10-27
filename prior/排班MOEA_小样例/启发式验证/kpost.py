"""K 个值班点（每槽需 K 台不同机器人）的推广版本：精确 DFS / 全时域随机键 GA / JSP 式循环编码 + 前缀搜索"""
import itertools, math
import numpy as np
from slotmodel import canon

def start(R): return tuple([(0, 0)] * len(R))

def stepK(R, st, chosen):
    new = [(max(w - 1, 0), ph) for (w, ph) in st]
    for j in chosen:
        rb = R[j]; ph = st[j][1]
        new[j] = (rb.g[ph] - 1, rb.nph(ph))
    return tuple(new)

def succK(R, st, K):
    av = [j for j, (w, _) in enumerate(st) if w == 0]
    # 同组同相位视为相同 → 按类别多重集组合去重
    cls = {}
    for j in av: cls.setdefault((R[j].group, st[j][1]), []).append(j)
    keys = list(cls); out = []
    def rec(i, left, pick):
        if left == 0: out.append(pick[:]); return
        if i == len(keys): return
        lst = cls[keys[i]]
        for c in range(min(left, len(lst)), -1, -1):
            pick.extend(lst[:c]); rec(i + 1, left - c, pick); del pick[len(pick) - c:]
    rec(0, K, [])
    return [stepK(R, st, ch) for ch in out]

# ---------------- 精确 DFS 找环
def exact(R, K, node_limit):
    s0 = start(R); color = {canon(R, s0): 1}; nodes = 0
    stack = [(s0, iter(succK(R, s0, K)))]
    while stack:
        st, it = stack[-1]
        for t in it:
            nodes += 1
            if nodes > node_limit: return None, nodes
            c = canon(R, t); col = color.get(c, 0)
            if col == 1: return True, nodes
            if col == 0:
                color[c] = 1; stack.append((t, iter(succK(R, t, K)))); break
        else:
            color[canon(R, st)] = 2; stack.pop()
    return False, nodes

# ---------------- 基线：全时域随机键 GA（只能证明有限 H）
def fullkey(R, K, H, budget, rng, pop=30):
    n = H * K
    def ev(x):
        st = start(R)
        for k in range(H):
            av = [j for j, (w, _) in enumerate(st) if w == 0]
            if len(av) < K: return k
            order = sorted(range(len(av)), key=lambda i: x[k * K + (i % K)] + i * 1e-3 * x[(k * K + i) % n])
            # 用 K 个基因各选一台（不放回）
            ch = []; pool = av[:]
            for q in range(K):
                i = min(int(x[k * K + q] * len(pool)), len(pool) - 1); ch.append(pool.pop(i))
            st = stepK(R, st, ch)
        return H
    X = rng.random((pop, n)); F = np.array([ev(x) for x in X]); used = pop
    while used < budget and F.max() < H:
        K_ = []
        for _ in range(pop):
            a, b = [max(rng.integers(pop, size=2), key=lambda i: F[i]) for _ in range(2)]
            c = np.where(rng.random(n) < 0.5, X[a], X[b])
            m = rng.random(n) < 2.0 / n
            K_.append(np.clip(c + m * rng.normal(0, 0.2, n), 0, 1 - 1e-9))
        K_ = np.array(K_); FK = np.array([ev(c) for c in K_]); used += pop
        A = np.vstack([X, K_]); B = np.concatenate([F, FK]); o = np.argsort(-B)[:pop]; X, F = A[o], B[o]
    return F.max() >= H, used

# ---------------- JSP 式循环编码（周期） + 前缀搜索（从全新状态进入该周期）
def choose_m(R, K, Pmax=120, mode='slack'):
    """选每台在一个周期内出场的超周期数 m_j（可为 0 = 不用这台）：
    约束 sum_j m_j*Lp_j = K*P（每槽恰 K 台），m_j*S_j <= P（单台占用不超过周期）。
    在 P<=Pmax 中选"最大单台利用率 max_j m_j*S_j/P"最小的方案（留出最多余量），平手取较小 P。DP 求解。"""
    n = len(R); Lp = [r.L - r.mu for r in R]; Sp = [sum(r.g[r.mu:]) for r in R]
    best = None
    for P in range(1, Pmax + 1):
        T = K * P; reach = {0: (0.0, [])}       # 和 -> (当前最大利用率, m 列表)
        for j in range(n):
            mmax = P // Sp[j]; nxt = {}
            for tot, (u, ms) in reach.items():
                for mm in range(mmax + 1):
                    t2 = tot + mm * Lp[j]
                    if t2 > T: break
                    u2 = max(u, mm * Sp[j] / P)
                    if t2 not in nxt or u2 < nxt[t2][0]: nxt[t2] = (u2, ms + [mm])
            reach = nxt
        if T in reach and (best is None or reach[T][0] < best[0] - 1e-9):
            best = (reach[T][0], reach[T][1], P)
            if mode == 'minP': break                 # 最短周期：第一个可行 P 即停
    return (best[1], best[2]) if best else (None, None)

def make_viol(R, K, P):
    def viol(seq):
        pos = {}
        for i, j in enumerate(seq): pos.setdefault(j, []).append(i // K)
        tot = 0; bad = []; phis = {}
        for j, ps in pos.items():
            rb = R[j]; per = rb.g[rb.mu:]; k = len(ps)
            gaps = [P if k == 1 else (ps[(i + 1) % k] - ps[i]) % P for i in range(k)]
            if k > 1: gaps[-1] = gaps[-1] or P   # 环回
            best = None
            for phi in range(len(per)):
                v = 0; lst = []
                for i in range(k):
                    need = per[(phi + i) % len(per)]
                    if gaps[i] < need: v += need - gaps[i]; lst.append(i)
                if best is None or v < best[0]: best = (v, lst, phi)
            tot += best[0]; phis[j] = best[2]
            idx = [i for i, x in enumerate(seq) if x == j]
            bad += [idx[i] for i in best[1]]
        return tot, bad, phis
    return viol

def jit_seed(R, K, m, P, rng, noise=0.3):
    """准时化(JIT)序列构造：按"应出场次数 − 已出场次数"的欠账从大到小、在可用机器人中选 K 台；
    跑 3 个周期取最后一个，再把出场次数修补到恰好 c_j。用作初始种群（结构化播种）。"""
    n = len(R); c = [m[j] * (R[j].L - R[j].mu) for j in range(n)]
    st = list(start(R)); cnt = [0] * n; hist = []
    for t in range(3 * P):
        av = [j for j in range(n) if st[j][0] == 0 and c[j] > 0]
        if len(av) < K:
            av = sorted(range(n), key=lambda j: (st[j][0], -c[j]))[:max(K, len(av))]
        av.sort(key=lambda j: -(c[j] * (t + 1) / P - cnt[j] + noise * rng.random()))
        ch = av[:K]
        for j in ch: cnt[j] += 1
        hist.append(ch); st = list(stepK(R, tuple(st), ch))
    seq = [j for sl in hist[-P:] for j in sl]
    have = [seq.count(j) for j in range(n)]
    over = [j for j in range(n) for _ in range(max(0, have[j] - c[j]))]
    under = [j for j in range(n) for _ in range(max(0, c[j] - have[j]))]
    rng.shuffle(under)
    for j in over:
        i = max(i for i, x in enumerate(seq) if x == j) if rng.random() < 0.5 else min(i for i, x in enumerate(seq) if x == j)
        seq[i] = under.pop()
    return seq

def cyclic(R, K, budget, rng, pop=30, ls=20, seeded=True):
    """先试最短周期（预算一半），不行再试余量最大的周期"""
    used = 0
    for mode, b in [('minP', budget // 2), ('slack', budget - budget // 2)]:
        res, u = _cyclic(R, K, b, rng, pop, ls, seeded, mode); used += u
        if res is not None: return res, used
    return None, used

def _cyclic(R, K, budget, rng, pop, ls, seeded, mode):
    m, P = choose_m(R, K, mode=mode)
    if m is None: return None, 0
    base = [j for j in range(len(R)) for _ in range(m[j] * (R[j].L - R[j].mu))]
    viol = make_viol(R, K, P); n = len(R); Ln = len(base)
    def pox(a, b):
        keep = set(rng.choice(n, size=max(1, n // 2), replace=False).tolist())
        fill = iter([x for x in b if x not in keep])
        return [x if x in keep else next(fill) for x in a]
    X = [jit_seed(R, K, m, P, rng) if (seeded and q < pop // 2) else list(rng.permutation(base)) for q in range(pop)]
    F = [viol(x) for x in X]; used = pop
    while used < budget:
        bi = min(range(len(F)), key=lambda i: F[i][0])
        if F[bi][0] == 0: return (X[bi], F[bi][2], P), used
        kids = []
        for _ in range(pop // 2):
            a, b = [min(rng.integers(pop, size=2), key=lambda i: F[i][0]) for _ in range(2)]
            c = pox(X[a], X[b])
            if rng.random() < 0.5:
                i, k = rng.integers(Ln, size=2); c.insert(k, c.pop(i))
            fc = viol(c); used += 1
            for _ in range(ls):      # 冲突引导局部搜索：违约的那一班，试插到若干随机位置，取最好的（不变差即接受）
                if fc[0] == 0 or used >= budget: break
                i = rng.choice(fc[1]); best = None
                for k in rng.choice(Ln, size=min(Ln, 8), replace=False):
                    d = c[:]; d.insert(int(k), d.pop(i)); fd = viol(d); used += 1
                    if best is None or fd[0] < best[1][0]: best = (d, fd)
                if best[1][0] <= fc[0]: c, fc = best
            kids.append((c, fc))
        allp = list(zip(X, F)) + kids; allp.sort(key=lambda z: z[1][0])
        X = [z[0] for z in allp[:pop]]; F = [z[1] for z in allp[:pop]]
    return None, used

def cycle_targets(R, K, seq, phis, P):
    """周期每个旋转起点处的状态 (w, ph)：等待=上一班使用后的剩余等待；相位=本周期第一次出场的相位"""
    slots = [[] for _ in range(P)]
    for i, j in enumerate(seq): slots[i // K].append(j)
    app = {j: [s for s in range(P) for jj in slots[s] if jj == j] for j in range(len(R))}
    tg = []
    for s0 in range(P):
        st = []
        for j, rb in enumerate(R):
            per = rb.g[rb.mu:]; ps = app[j]; k = len(ps)
            if k == 0: st.append((0, 0)); continue
            # 按旋转后的时间顺序排列出场，相位从 phis[j] 开始对应原顺序
            order = sorted(range(k), key=lambda i: (ps[i] - s0) % P)
            first = order[0]; last = order[-1]
            ph_first = (phis[j] + first) % len(per); ph_last = (phis[j] + last) % len(per)
            t_last = (ps[last] - s0) % P - P          # 上一周期中的最后一次出场（负时刻）
            w = max(0, t_last + per[ph_last])
            st.append((w, rb.mu + ph_first))
        tg.append(st)
    return tg

def dominates(R, st, tg):
    """同组内可重新配对：相同相位，等待不大于目标 → 可完全模仿目标此后的排班（等待是免费的）"""
    from collections import defaultdict
    a = defaultdict(list); b = defaultdict(list)
    for rb, s in zip(R, st): a[(rb.group, s[1])].append(s[0])
    for rb, s in zip(R, tg): b[(rb.group, s[1])].append(s[0])
    if set(a) != set(b) or any(len(a[k]) != len(b[k]) for k in a): return False
    return all(x <= y for k in a for x, y in zip(sorted(a[k]), sorted(b[k])))

def _dist(R, st, tg, used):
    d = 0
    for j in used:
        rb = R[j]; (w, ph), (w2, ph2) = st[j], tg[j]
        Lp = rb.L - rb.mu
        dph = (rb.mu - ph) + (ph2 - rb.mu) if ph < rb.mu else (ph2 - ph) % Lp
        d += dph * 10 + max(0, w - w2)
    return d

def dominates_used(R, st, tg, used):
    return all(st[j][1] == tg[j][1] and st[j][0] <= tg[j][0] for j in used)

def prefix_beam(R, K, targets, node_limit, used=None, width=300, max_depth=800):
    """束搜索：从全新状态出发，到达"支配"某个周期旋转状态的状态（参与周期的机器人相位一致、等待不更长；等待是免费的）"""
    used = list(range(len(R))) if used is None else used
    beam = [start(R)]; nodes = 0; seen = set()
    for d in range(max_depth):
        for st in beam:
            for tg in targets:
                if dominates_used(R, st, tg, used): return True, nodes
        cand = []
        for st in beam:
            for t in succK(R, st, K):
                nodes += 1
                c = canon(R, t)
                if c in seen: continue
                seen.add(c)
                cand.append((min(_dist(R, t, tg, used) for tg in targets), t))
            if nodes > node_limit: return None, nodes
        if not cand: return False, nodes
        cand.sort(key=lambda z: z[0]); beam = [t for _, t in cand[:width]]
    return None, nodes

def verify_cycle(R, K, seq, phis, P, tg, s0=0, periods=3):
    """从旋转起点 s0 的目标状态出发执行若干个周期：每次派出的机器人必须可用，且每个周期结束回到同一状态"""
    st = tuple(tg); slots = [list(seq[i * K:(i + 1) * K]) for i in range(P)]
    for _ in range(periods):
        for q in range(P):
            sl = slots[(s0 + q) % P]
            if any(st[j][0] != 0 for j in sl): return False
            st = stepK(R, st, sl)
        if st != tuple(tg): return False
    return True


def prefix_dfs(R, K, targets, node_limit, used):
    """深度优先 + 失败状态记忆，子节点按到目标的距离排序"""
    s0 = start(R); dead = set(); nodes = 0
    def kids(st):
        ks = succK(R, st, K)
        return iter(sorted(ks, key=lambda t: min(_dist(R, t, tg, used) for tg in targets)))
    stack = [(s0, kids(s0))]; onpath = {canon(R, s0)}
    while stack:
        st, it = stack[-1]
        for t in it:
            nodes += 1
            if nodes > node_limit: return None, nodes
            if any(dominates_used(R, t, tg, used) for tg in targets): return True, nodes
            c = canon(R, t)
            if c in dead or c in onpath: continue
            onpath.add(c); stack.append((t, kids(t))); break
        else:
            c = canon(R, st); dead.add(c); onpath.discard(c); stack.pop()
    return False, nodes

def prefix(R, K, targets, node_limit, used=None):
    """先束搜索，失败再深搜"""
    used = list(range(len(R))) if used is None else used
    if any(dominates_used(R, start(R), tg, used) for tg in targets): return True, 0
    ok, n1 = prefix_beam(R, K, targets, node_limit // 2, used)
    if ok: return True, n1
    ok, n2 = prefix_dfs(R, K, targets, node_limit - n1, used)
    return ok, n1 + n2
