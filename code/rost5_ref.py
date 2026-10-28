"""rost5 的 Python 参考实现（只用于一致性核对，速度不重要）。与 code/cpp/rost5.cpp 的 run() 逐步对应：
- 岗位：站 s 拆成 K_s 个岗位；
- 策略：R_LP（rlp 文件 + 可选 planf）或 rost3 式线性打分（8 维权重，相位缺额 tot/L − D_i）；
- 初值：全新启动，或 rlp 文件里 steady=1 时的稳态初值；
- 判环：状态第一次重复。返回 gaps、cert（1 证明、0 环内有缺岗、−1 未判出）、period。
状态里的整数与 C++ 相同，浮点运算的顺序也与 C++ 相同。
T7 补充：线性打分支持 curdef（相位缺额 = 当前相位人数 − N_k/L）和 LP 引导项（w13·lpz[k][s]）；RDD 解码（rdd 参数，见 decode_geno）；
返回缺岗归因 why_r / why_t（与 C++ 相同口径）和窗口内决策数 ndec_w；size_res / size_rule / near_miss 逐步对应 C++ 的定台数循环。"""
import os, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import U, ToolInst  # noqa: E402

CL = 40


def llround(x):
    """C 的 llround（四舍五入、远离 0）；这里只用于非负数。"""
    r = math.floor(x)
    return int(r) + (1 if x - r >= 0.5 else 0)


def parse_rlp(I, path):
    """读规则文件，结构与 C++ parse_rlp 相同。"""
    nk, ns = len(I.T), len(I.sites)
    R = dict(H=0, steady=0, noscale=0, NLPq=None, N0=None, M0=None, Q=None, CP=None,
             p=[[0] * ns for _ in range(nk)], CQ=[[0] * ns for _ in range(nk)],
             F={}, FP={}, W={})
    for line in open(path):
        a = line.split()
        if not a: continue
        tag = a[0]
        if tag == "H": R["H"] = int(a[1])
        elif tag == "steady": R["steady"] = int(a[1])
        elif tag == "noscale": R["noscale"] = int(a[1])      # T3b 的 _u 变体：三种计划都不按台数缩放
        elif tag in ("NLPq", "N0", "M0", "Q", "CP"): R[tag] = [int(x) for x in a[1:]]
        elif tag == "p": k, s, p, c = map(int, a[1:5]); R["p"][k][s] = p; R["CQ"][k][s] = c
        elif tag == "W": R["W"][int(a[1])] = [int(x) for x in a[2:]]
        elif tag == "F": R["F"][(int(a[1]), int(a[2]))] = [int(x) for x in a[3:]]
        elif tag == "FP": R["FP"][(int(a[1]), int(a[2]))] = [int(x) for x in a[3:]]
    return R


def parse_planf(I, path):
    """与 C++ parse_planf 相同：#N 行给 LP 台数；其余 k i t0 流量，按 (t0, 流量) 排序。"""
    nk = len(I.T); planN = [0.0] * nk
    pf = [[[] for _ in range(I.L[k])] for k in range(nk)]
    for line in open(path):
        line = line.rstrip("\n")
        if not line: continue
        if line[0] == "#":
            vals = line[2:].split()
            for k in range(nk): planN[k] = float(vals[k])
            continue
        k, i, t0, x = line.split(); k, i, t0, x = int(k), int(i), int(t0), float(x)
        if k < nk and i < I.L[k]: pf[k][i].append((t0, x))
    for k in range(nk):
        for i in range(I.L[k]): pf[k][i].sort()
    return planN, pf


def decode_geno(g, W, ndev, qon, reset_opt):
    """与 C++ decode_geno 相同：槽 (p, a) 在 a ≥ 1 − q_on 时激活，作用于第 ⌊p·N_dev⌋ 个窗口决策，u = (a − 1 + q_on)/q_on；后一个槽覆盖前一个。"""
    dev = [-1.0] * max(ndev, 0); maxpos = -1; thr = 1.0 - qon
    for a in range(0, len(g) - 1, 2):
        p, x = g[a], g[a + 1]
        if ndev <= 0 or x < thr: continue
        pos = int(math.floor(p * ndev))
        if pos >= ndev: pos = ndev - 1
        dev[pos] = (x - 1.0 + qon) / qon; maxpos = max(maxpos, pos)
    return dict(W=W, ndev=ndev, dev=dev, maxpos=maxpos, reset_opt=reset_opt)


def run(I, Nk, Mj, Tmax, certify=True, rlp=None, planf=None, w=None, gap_log=None, trace=None, curdef=False, lpz=None, rdd=None):
    """仿真一次。trace 不为 None 时（T3c 诊断用，dict）只做记录、不改变任何决策：
    trace["plan"] = 换算后的提前大修计划 {(k, i): [[t0, 次数], …]}；trace["disp"] = 每次派工
    (t, 站, 岗位, 机器人, 型, 相位, 是否按计划提前大修, 新的就绪基准 rv, 值岗结束时刻)；
    trace["gaps"] = 每次缺岗 (t, 站, 就绪的候选数, 归因：("r", 可在该站值岗的各型) 或 ("t", 首个候选型所需工具))，与 C++ 的 why_r / why_t 口径相同；
    trace["robots"] = 各机器人的型。trace["light"] 为真时不记派工（长时域判环时省内存）。"""
    TI = ToolInst(I)
    nk, ns, nt = len(I.T), len(I.sites), TI.nt
    v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]; c0 = max(v0)
    lead = [I.lead(s) for s in range(ns)]
    st = {(k, i, s): I.step(k, i, s) for k in range(nk) for i in range(I.L[k]) for s in I.elig[k]}
    pre = {(k, s): TI.pre(k, s) for k in range(nk) for s in I.elig[k]}
    elig = [set(I.elig[k]) for k in range(nk)]
    nelig = [sum(1 for k in range(nk) if s in elig[k]) for s in range(ns)]
    flex = [len(I.elig[k]) for k in range(nk)]
    R = rlp
    # 提前大修计划（按实际台数换算）
    pl = []
    if planf is not None:
        planN, pf = planf
        for k in range(nk):
            row = []
            sc = Nk[k] / planN[k] if planN[k] > 1e-9 else 0.0
            if R is not None and R["noscale"]: sc = 1.0      # _u：LP 累计流直接四舍五入成整数次
            for i in range(I.L[k]):
                F = 0.0; prev = 0; lst = []
                for (t0, x) in pf[k][i]:
                    F += x * sc; c = llround(F) - prev
                    if c > 0: lst.append([t0, c]); prev += c
                row.append(lst)
            pl.append(row)
    if trace is not None:
        trace["plan"] = {(k, i): [list(e) for e in pl[k][i]] for k in range(len(pl)) for i in range(len(pl[k])) if pl[k][i]}
        trace["disp"] = []; trace["gaps"] = []
    # 机器人
    rk, rph, rv = [], [], []
    if R is not None and R["steady"]:
        for k in range(nk):
            L = I.L[k]; W = sum(R["W"][k])
            n = [Nk[k] * R["W"][k][i] // W for i in range(L)]
            rem = [Nk[k] * R["W"][k][i] % W for i in range(L)]
            ordr = sorted(range(L), key=lambda i: -rem[i])        # 稳定排序，与 C++ stable_sort 相同
            for e in range(Nk[k] - sum(n)): n[ordr[e]] += 1
            for i in range(L):
                for j in range(n[i]): rk.append(k); rph.append(i); rv.append(j * R["W"][k][i] // n[i])
    else:
        for k in range(nk):
            for _ in range(Nk[k]): rk.append(k); rph.append(0); rv.append(v0[k])
    nr = len(rk)
    if trace is not None: trace["robots"] = list(rk); trace["rv0"] = list(rv)
    tf = [[-1000000] * Mj[j] for j in range(nt)]
    psite = [s for s in range(ns) for _ in range(I.sites[s]["K"])]
    np_ = len(psite)
    fr = [c0 + lead[psite[q]] for q in range(np_)]
    D = [[0] * I.L[k] for k in range(nk)]; tot = [0] * nk
    qmax = 0
    if planf is not None:
        for k in range(nk):
            for i in range(I.L[k]):
                for (t0, x) in planf[1][k][i]: qmax = max(qmax, t0)
    lag = [[0] * ns for _ in range(nk)]; phl = [[0] * I.L[k] for k in range(nk)]
    H = R["H"] if R is not None else 0
    if R is not None: qmax = max(qmax, H)

    def tab(k, s, tt):
        """t < H 时 k 型在 s 站的计划累计班数（单位 1/Q_s）；默认按 Nk/N_LP 缩放，noscale 时不缩放。"""
        if (k, s) not in R["F"]: return 0
        if R["noscale"]: return R["Q"][s] * R["F"][(k, s)][tt] // 1000000
        if R["NLPq"][k] <= 0: return 0
        return R["Q"][s] * Nk[k] * R["F"][(k, s)][tt] // R["NLPq"][k]

    def planQ(k, s, t):
        if t < 0: return 0
        if t < H: return tab(k, s, t)
        B = tab(k, s, H - 1) if H > 0 else 0
        return B + R["p"][k][s] * (t - H + 1)

    def phplanQ(k, i, t):
        if t < 0 or H == 0: return 0
        tt = min(t, H - 1)
        if R["noscale"]: return I.L[k] * R["FP"][(k, i)][tt] // 1000000
        if R["NLPq"][k] <= 0: return 0
        return I.L[k] * Nk[k] * R["FP"][(k, i)][tt] // R["NLPq"][k]

    clip = lambda x, c: c if x > c else (-c if x < -c else x)
    tprev = -1
    gaps = 0; resets = 0; seen = {}
    cnt = [[0] * I.L[k] for k in range(nk)]
    for k in range(nk): cnt[k][0] = Nk[k]
    why_r = [0] * nk; why_t = [0] * nt; ndec = 0

    def tool_ok(r, s):
        """机器人 r 此刻去 s 站是否有足够的空闲适配工具（不改状态）。"""
        if nt == 0: return True
        k = rk[r]; need = t - pre[(k, s)]
        return sum(1 for j in TI.R[k]["C"] for f in tf[j] if f <= need) >= TI.R[k]["q"]

    def pick_dev(order, s, u):
        """与 C++ pick_dev 相同：−1 不偏离；≥ 0 改派的机器人；−2 默认机器人本班后提前大修。"""
        fe = [r for r in order if tool_ok(r, s)]
        if not fe: return -1
        r0 = fe[0]; ro = rdd["reset_opt"] and rph[r0] != I.L[rk[r0]] - 1
        nA = len(fe) - 1 + (1 if ro else 0)
        if nA <= 0: return -1
        a = int(math.floor(u * nA))
        if a >= nA: a = nA - 1
        return fe[a + 1] if a < len(fe) - 1 else -2

    def out(**kw):
        """返回值：附上归因与窗口决策数。"""
        kw.update(why_r=why_r, why_t=why_t, ndec_w=ndec); return kw
    t = min(fr)
    while t < Tmax:
        if R is not None:
            for k in range(nk):
                for s in I.elig[k]: lag[k][s] = clip(lag[k][s] + planQ(k, s, t) - planQ(k, s, tprev), R["CQ"][k][s])
                for i in range(I.L[k]): phl[k][i] = clip(phl[k][i] + phplanQ(k, i, t) - phplanQ(k, i, tprev), R["CP"][k])
            tprev = t
        if certify and t >= qmax and (rdd is None or t >= rdd["W"] or ndec > rdd["maxpos"]):
            cur = []
            for r in range(nr): cur += [rk[r], rph[r], max(rv[r] - t, -CL)]
            for row in pl:
                for lst in row:
                    for e in lst: cur.append(e[1])
            for q in range(np_): cur.append(fr[q] - t)
            if R is not None:
                for k in range(nk):
                    for s in I.elig[k]: cur.append(lag[k][s])
                    cur += phl[k]
            elif not curdef:
                for k in range(nk):
                    for i in range(I.L[k]): cur.append(tot[k] - I.L[k] * D[k][i])
            for j in range(nt): cur += sorted(max(f - t, -CL) for f in tf[j])
            key = tuple(cur)
            if key in seen:
                t1 = seen[key]
                return out(gaps=gaps, cert=1 if gaps == 0 else 0, period=t - t1, enter=t1, resets=resets)
            if len(seen) >= (1 << 20): return out(gaps=gaps, cert=-1, period=-1, enter=-1, resets=resets)
            seen[key] = t
        due = sorted([q for q in range(np_) if fr[q] <= t], key=lambda q: nelig[psite[q]])
        for q in due:
            s = psite[q]
            dv = -1.0; force = False
            if rdd is not None and t < rdd["W"]:
                c = ndec; ndec += 1
                if c < rdd["ndev"]: dv = rdd["dev"][c]
            if R is not None:
                cand = sorted((-lag[rk[r]][s], rk[r], -phl[rk[r]][rph[r]], max(rv[r] - t, -CL), r)
                              for r in range(nr) if s in elig[rk[r]] and rv[r] + lead[s] <= t)
                order = [c[4] for c in cand]
            else:
                def score(r):
                    k, i = rk[r], rph[r]
                    sc = w[0] * max(rv[r] - t, -CL)
                    sc = sc + w[1] * st[(k, i, s)][1] / I.d[k]
                    sc = sc + w[2] * flex[k]
                    deficit = (cnt[k][i] - Nk[k] / I.L[k]) if curdef else (tot[k] / I.L[k] - D[k][i])
                    sc = sc - w[3] * deficit
                    sc = sc + w[4] * I.d[k]
                    sc = sc + w[5] * i / I.L[k]
                    if lpz is not None: sc = sc - w[13] * lpz[k][s]
                    return (sc, r)
                cand = sorted(score(r) for r in range(nr) if s in elig[rk[r]] and rv[r] + lead[s] <= t)
                order = [c[1] for c in cand]
            if dv >= 0:
                pk = pick_dev(order, s, dv)
                if pk == -2: force = True
                elif pk >= 0: order.remove(pk); order.insert(0, pk)
            done = False
            for r in order:
                k, i = rk[r], rph[r]
                need = t - pre[(k, s)]
                opts = []
                for j in TI.R[k]["C"]:
                    for idx, f in enumerate(tf[j]):
                        if f <= need:
                            rel = max(f - t, -CL)
                            opts.append((TI.TT[j]["maintain"] if R is not None else w[6] * TI.TT[j]["maintain"] + w[7] * rel, rel, j, idx))
                q_ = TI.R[k]["q"]
                if len(opts) < q_: continue
                opts.sort(); ch = opts[:q_]
                unin = 0.0
                for c in ch: unin += TI.TT[c[2]]["uninstall"]
                for c in ch: tf[c[2]][c[3]] = t + I.d[k] + TI.post(k, s, c[2], unin)
                ld, nx = st[(k, i, s)]; nph = (i + 1) % I.L[k]
                planned = False
                if pl and i != I.L[k] - 1:
                    for e in pl[k][i]:
                        if e[0] > t: break
                        if e[1] > 0: e[1] -= 1; planned = True; break
                if planned: nx = st[(k, I.L[k] - 1, s)][1]; nph = 0; resets += 1
                if force and nph != 0: nx = st[(k, I.L[k] - 1, s)][1]; nph = 0; resets += 1
                rv[r] = t - ld + nx; rph[r] = nph; cnt[k][i] -= 1; cnt[k][nph] += 1
                if trace is not None and not trace.get("light"): trace["disp"].append((t, s, q, r, k, i, planned, rv[r], t + I.d[k]))
                D[k][i] += 1; tot[k] += 1
                if R is not None:
                    lag[k][s] = clip(lag[k][s] - R["Q"][s], R["CQ"][k][s])
                    phl[k][i] -= I.L[k]
                    if t >= H:
                        for i2 in range(I.L[k]): phl[k][i2] += 1
                    for i2 in range(I.L[k]): phl[k][i2] = clip(phl[k][i2], R["CP"][k])
                fr[q] = t + I.d[k]; done = True; break
            if not done:
                gaps += 1; fr[q] = t + 1
                if order:
                    for j in TI.R[rk[order[0]]]["C"]: why_t[j] += 1
                else:
                    for k in range(nk):
                        if s in elig[k]: why_r[k] += 1
                if trace is not None:
                    why = ("t", list(TI.R[rk[order[0]]]["C"])) if order else ("r", [k for k in range(nk) if s in elig[k]])
                    trace["gaps"].append((t, s, len(order), why))
                if gap_log is not None and len(gap_log) < 50:
                    gap_log.append(dict(t=t, s=s, cand=len(order), lag={k: lag[k][s] for k in range(nk) if s in elig[k]}))
        t = min(fr)
    return out(gaps=gaps, cert=-1, period=-1, enter=-1, resets=resets)


CERT_PATIENCE = 5


def size_res(I, pol, Hs, Tcert, N0, M0, certadd=1, max_add=60):
    """逐步对应 C++ size_res(mode = certs)：补足（试加取缺岗最少）→ 回删 → 证书阶段。
    证书阶段补台：certadd=1（默认）试加取 Tcert 格内缺岗最少，当前台数的 Tcert 格缺岗连续 5 次没有刷新最低值时判失败；
    2 = 试加、不设停止条件；0 = 直接加归因最大。
    pol 是传给 run 的策略参数（rlp / planf / w / curdef / lpz / rdd）。返回 (R, Nk, Mj, 状态, 周期)。"""
    nk = len(I.T); Nk, Mj = list(N0), list(M0); nres = nk + len(Mj)

    def bump(idx, d):
        if idx < nk: Nk[idx] += d
        else: Mj[idx - nk] += d
    r = run(I, Nk, Mj, Hs, False, **pol); added = 0
    while r["gaps"] > 0 and added < max_add:
        sig = r["why_r"] + r["why_t"]; bg = bi = bs = None
        for idx in range(nres):
            bump(idx, 1); g = run(I, Nk, Mj, Hs, False, **pol)["gaps"]; bump(idx, -1)
            if bi is None or g < bg or (g == bg and sig[idx] > bs): bg, bi, bs = g, idx, sig[idx]
        bump(bi, 1); added += 1; r = run(I, Nk, Mj, Hs, False, **pol)
    if r["gaps"] > 0: return -1, Nk, Mj, "fail", -1
    improved = True
    while improved:
        improved = False
        for idx in range(nres):
            v = Nk[idx] if idx < nk else Mj[idx - nk]
            if v == 0: continue
            bump(idx, -1)
            if run(I, Nk, Mj, Hs, False, **pol)["gaps"] == 0: improved = True; break
            bump(idx, 1)
    best_g0, stall = None, 0
    for _ in range(max_add):
        c = run(I, Nk, Mj, Tcert, True, **pol)
        if c["cert"] == -1 and c["gaps"] == 0: c = run(I, Nk, Mj, 4 * Tcert, True, **pol)
        if c["cert"] == 1: return sum(Nk) + sum(Mj), Nk, Mj, "cert", c["period"]
        if c["cert"] == -1 and c["gaps"] == 0: return sum(Nk) + sum(Mj), Nk, Mj, "uncert", -1
        sig = c["why_r"] + c["why_t"]
        if certadd == 0:
            bi = max(range(nres), key=lambda a: (sig[a], -a)); bump(bi, 1); continue
        if certadd == 1:
            g0 = run(I, Nk, Mj, Tcert, False, **pol)["gaps"]
            if best_g0 is None or g0 < best_g0: best_g0, stall = g0, 0
            else:
                stall += 1
                if stall >= CERT_PATIENCE: return -1, Nk, Mj, "fail", -1
        bg = bi = bs = None
        for idx in range(nres):
            bump(idx, 1); g = run(I, Nk, Mj, Tcert, False, **pol)["gaps"]; bump(idx, -1)
            if bi is None or g < bg or (g == bg and sig[idx] > bs): bg, bi, bs = g, idx, sig[idx]
        bump(bi, 1)
    return -1, Nk, Mj, "fail", -1


def size_rule(I, pol, Hs, Tcert, N0, M0, is_rlp, certadd=1, facloop=False):
    """对应 C++ size_rule：R_LP（或 facloop 的线性规则）从 (N0, M0) × 1.00 / 1.15 / 1.30 / 1.50 起，取第一个成功的；
    其他线性规则只从 (N0, M0) 起一次。返回 size_res 的结果加起点倍数。"""
    if not is_rlp and not facloop: return size_res(I, pol, Hs, Tcert, N0, M0, certadd) + (100,)
    last = None
    for fac in (100, 115, 130, 150):
        r = size_res(I, pol, Hs, Tcert, [(x * fac + 99) // 100 for x in N0], [(x * fac + 99) // 100 for x in M0], certadd)
        if r[0] >= 0: return r + (fac,)
        last = r
    return last + (-1,)


def near_miss(I, pol, Nk, Mj, Hs):
    """对应 C++ near_miss：任减 1 个资源后 Hs 格内缺岗数的最小值。"""
    nk = len(Nk); best = None; Nk, Mj = list(Nk), list(Mj)
    for idx in range(nk + len(Mj)):
        arr, a = (Nk, idx) if idx < nk else (Mj, idx - nk)
        if arr[a] == 0: continue
        arr[a] -= 1; g = run(I, Nk, Mj, Hs, False, **pol)["gaps"]; arr[a] += 1
        best = g if best is None else min(best, g)
    return 0 if best is None else best
