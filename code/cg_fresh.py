"""T2：全新启动时间展开 LP（freshlb.fresh_lb 的模型）的列生成求解。

主问题（与 freshlb.py 相同的约束，只是机器人弧流换成路径列）：
  min Σ_k c_k Σ_{p∈P_k} λ_p + Σ_j M_j + BIG·Σ 人工松弛
  覆盖   Σ_p a_p(s,t) λ_p + 松弛_{s,t} = K_s          (s,t)：t ∈ [c0 + lead_s, H)
  耦合   Σ_j y_{k,s,D,j} − q_k Σ_p b_p(k,s,D) λ_p = 0  每个值岗键 (k,s,D)，D < H
  工具流 工具节点 (j,v) 守恒，M_j 在 (j,0) 注入        （工具弧全部显式放在主问题里，规模只和 m·H 成正比）
列：k 型单台机器人从 (相位 0, v0_k) 出发、越过 H 进入汇点的一条路径；N_k = Σ_{p∈P_k} λ_p。
在 DAG 上路径分解与弧流等价，所以列生成收敛后的值等于 freshlb 整体 LP 的值。
定价：c_k − max_p Σ_{值岗弧} [Σ_{覆盖格} π_{s,t} − q_k μ_{k,s,D}]，最长路在 code/cpp/price.cpp。
初始只有人工松弛（BIG 足够大，收敛后检查松弛为 0）。每轮每型最多取 ncol 条列：取出一条后把它覆盖格上的正对偶清零再求下一条。
主问题用 highspy，加列后从上一轮基热启动。"""
import os, sys, json, math, time, ctypes, argparse
import numpy as np
from scipy.sparse import coo_matrix
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, U, ToolInst  # noqa: E402
import highspy  # noqa: E402

LIB = os.path.join(ROOT, "code", "cpp", "libprice.so")
BIG = 100.0


def _lib():
    """加载定价库；不存在就编译。"""
    if not os.path.exists(LIB):
        os.system(f"g++ -O2 -shared -fPIC -o {LIB} {os.path.join(ROOT, 'code', 'cpp', 'price.cpp')}")
    lib = ctypes.CDLL(LIB)
    P = np.ctypeslib.ndpointer
    lib.longest_path_w.restype = ctypes.c_int
    lib.longest_path_w.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                   P(np.int32, flags="C"), P(np.int32, flags="C"), P(np.int32, flags="C"), ctypes.c_int,
                                   P(np.float64, flags="C"), P(np.float64, flags="C"),
                                   P(np.int32, flags="C"), P(np.int32, flags="C"), P(np.int32, flags="C"), P(np.int32, flags="C"),
                                   ctypes.c_int]
    lib.longest_path.restype = ctypes.c_int
    lib.longest_path.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                 P(np.int32, flags="C"), P(np.int32, flags="C"), ctypes.c_int,
                                 P(np.float64, flags="C"), P(np.float64, flags="C"),
                                 P(np.int32, flags="C"), P(np.int32, flags="C"), P(np.int32, flags="C"), P(np.int32, flags="C"),
                                 ctypes.c_int]
    return lib


class FreshCG:
    """一个 (实例, H, 是否允许提前大修) 的列生成。W 不为 None 时（T6 的 LP_fr^Rw）提前大修弧只在值岗开始 D < W 时存在。"""

    def __init__(self, I, H, reset=False, W=None):
        self.I, self.H, self.reset, self.W = I, H, reset, W
        TI = self.TI = ToolInst(I)
        nk, ns, nt = len(I.T), len(I.sites), TI.nt
        self.nk, self.ns, self.nt = nk, ns, nt
        self.v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]
        self.c0 = max(self.v0)
        self.lead = [I.lead(s) for s in range(ns)]
        self.cost = [I.cost[t] for t in I.tn]
        # 覆盖行
        self.cov_first = [self.c0 + self.lead[s] for s in range(ns)]
        self.cov_row = {}
        for s in range(ns):
            for t in range(self.cov_first[s], H): self.cov_row[(s, t)] = len(self.cov_row)
        nc = len(self.cov_row)
        # 耦合行：所有 (k, s, D)，D ∈ [v0_k + lead_s, H)
        self.duty_key = {}
        for k in range(nk):
            for s in I.elig[k]:
                for D in range(self.v0[k] + self.lead[s], H): self.duty_key[(k, s, D)] = nc + len(self.duty_key)
        nq = len(self.duty_key)
        self.tool_row0 = nc + nq
        nrow = nc + nq + nt * H
        self.nc, self.nq, self.nrow = nc, nq, nrow
        # 静态列：工具等待弧、工具值岗弧、M_j、人工松弛
        R, C, V, cost = [], [], [], []
        col = 0
        trow = lambda j, v: self.tool_row0 + j * H + v
        for j in range(nt):
            for v in range(H):
                R.append(trow(j, v)); C.append(col); V.append(-1)
                if v + 1 < H: R.append(trow(j, v + 1)); C.append(col); V.append(1)
                cost.append(0.0); col += 1
        self.tool_arc = []   # (列, k, s, D, j)
        for (k, s, D), row in self.duty_key.items():
            Rk = TI.R[k]; pre = TI.pre(k, s)
            unin_lb = Rk["q"] * min(TI.TT[jj]["uninstall"] for jj in Rk["C"])
            for j in Rk["C"]:
                post = TI.post(k, s, j, unin_lb)
                a, b = max(0, D - pre), D + I.d[k] + post
                R.append(trow(j, a)); C.append(col); V.append(-1)
                if b < H: R.append(trow(j, b)); C.append(col); V.append(1)
                R.append(row); C.append(col); V.append(1)
                cost.append(0.0); self.tool_arc.append((col, k, s, D, j, a, b)); col += 1
        self.iM = col
        for j in range(nt):
            R.append(trow(j, 0)); C.append(col); V.append(1); cost.append(1.0); col += 1
        self.islack = col
        for (s, t), row in self.cov_row.items():      # 覆盖不足松弛（对偶箱上沿）
            R.append(row); C.append(col); V.append(1); cost.append(BIG); col += 1
        for (s, t), row in self.cov_row.items():      # 覆盖过量松弛（对偶箱下沿）
            R.append(row); C.append(col); V.append(-1); cost.append(BIG); col += 1
        self.nstatic = col
        A = coo_matrix((V, (R, C)), shape=(nrow, col)).tocsc()
        rhs = np.zeros(nrow)
        for (s, t), row in self.cov_row.items(): rhs[row] = I.sites[s]["K"]
        self.rhs = rhs
        h = self.h = highspy.Highs()
        h.setOptionValue("output_flag", False); h.setOptionValue("threads", 1)
        lp = highspy.HighsLp()
        lp.num_col_ = col; lp.num_row_ = nrow
        lp.col_cost_ = np.array(cost); lp.col_lower_ = np.zeros(col); lp.col_upper_ = np.full(col, highspy.kHighsInf)
        lp.row_lower_ = rhs.copy(); lp.row_upper_ = rhs.copy()
        lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
        lp.a_matrix_.start_ = A.indptr; lp.a_matrix_.index_ = A.indices; lp.a_matrix_.value_ = A.data
        lp.a_matrix_.num_col_ = col; lp.a_matrix_.num_row_ = nrow
        h.passModel(lp)
        h.setOptionValue("simplex_strategy", 4)   # 原始单纯形：加列后上一轮的基仍原始可行，热启动最省（对偶单纯形实测很慢）
        # 定价用的步长表
        self.lib = _lib()
        self.tab = []
        for k in range(nk):
            el = I.elig[k]; L = I.L[k]
            nx = np.array([[I.step(k, i, s)[1] for s in el] for i in range(L)], dtype=np.int32).ravel()
            nxR = np.array([I.step(k, L - 1, s)[1] for s in el], dtype=np.int32)
            self.tab.append((el, L, nx, nxR))
        self.cols = []   # 路径列：(k, [(v, s, i, r), ...])
        self.ncol_k = [0] * nk

    def column(self, k, duties):
        """路径 → (行下标, 系数)。"""
        I = self.I; rows = {}
        for (v, s, i, r) in duties:
            D = v + self.lead[s]
            for q in range(I.d[k]):
                row = self.cov_row.get((s, D + q))
                if row is not None: rows[row] = rows.get(row, 0.0) + 1.0
            row = self.duty_key.get((k, s, D))
            if row is not None: rows[row] = rows.get(row, 0.0) - self.TI.R[k]["q"]
        idx = np.array(sorted(rows), dtype=np.int32)
        return idx, np.array([rows[i] for i in idx])

    def weights(self, k, pi_cov, mu):
        """w[e][v] = 在 v 出发去站 elig[e] 值岗的收益。pi_cov: (ns, H) 覆盖对偶；mu: dict 或 (k,s)→长 H 数组。"""
        I, H = self.I, self.H
        el = self.tab[k][0]; d = I.d[k]; q = self.TI.R[k]["q"]
        W = np.zeros((len(el), H))
        for e, s in enumerate(el):
            ld = self.lead[s]
            cs = np.concatenate([[0.0], np.cumsum(np.concatenate([pi_cov[s], np.zeros(d + ld + 1)]))])
            Ds = np.arange(H) + ld                      # 值岗开始
            cov = cs[np.minimum(Ds + d, len(cs) - 1)] - cs[np.minimum(Ds, len(cs) - 1)]
            W[e] = cov - q * mu[(k, s)][np.minimum(Ds, H)]
        return np.ascontiguousarray(W.ravel())

    def price(self, k, W):
        """最长路；返回 (收益, [(v, s, i, r)])。"""
        el, L, nx, nxR = self.tab[k]; H = self.H
        mx = H + 2
        ov, oe, oi, orr = (np.zeros(mx, dtype=np.int32) for _ in range(4))
        best = np.zeros(1)
        if self.reset and self.W is not None:      # 窗口内才有提前大修弧
            ld = np.array([self.lead[s] for s in el], dtype=np.int32)
            n = self.lib.longest_path_w(L, H, self.v0[k], len(el), nx, nxR, ld, int(self.W), W, best, ov, oe, oi, orr, mx)
        else:
            n = self.lib.longest_path(L, H, self.v0[k], len(el), nx, nxR, int(self.reset), W, best, ov, oe, oi, orr, mx)
        assert n >= 0
        return best[0], [(int(ov[a]), el[oe[a]], int(oi[a]), int(orr[a])) for a in range(n)]

    def greedy_columns(self, prefs=("occ", "fifo"), n_rand=8, seed=0):
        """初始列：机队不限量的贪心排班。每站拆成 K_s 个岗位，岗位到点时从已就绪的允许机器人里选，
        没有就新开一台（全新启动，相位 0、就绪 v0_k）。每格每站恰好 K_s 台在岗，所以这些列使主问题可行、不用人工松弛。
        prefs 里每个规则各生成一套列：occ = 优先每格值岗占用少的机型；fifo = 优先就绪最早的机器人。"""
        I, H = self.I, self.H
        added = 0
        rng = np.random.default_rng(seed)
        rate0 = [np.mean([I.step(k, i, s)[1] for i in range(I.L[k]) for s in I.elig[k]]) / I.d[k] for k in range(self.nk)]
        for pref in list(prefs) + ["rand"] * n_rand:     # rand：机型偏好加随机扰动、选人随机，给主问题更多样的初始列
            rate = rate0 if pref != "rand" else [r * rng.uniform(0.7, 1.3) for r in rate0]
            robots = []            # [k, 相位, 就绪 v, 值岗列表]
            pos = [[s, self.cov_first[s]] for s in range(self.ns) for _ in range(I.sites[s]["K"])]
            while True:
                t = min(p[1] for p in pos)
                if t >= H: break
                for p in pos:
                    if p[1] != t: continue
                    s = p[0]; ld = self.lead[s]
                    cand = [r for r in robots if s in I.elig[r[0]] and r[2] + ld <= t]
                    if cand:
                        if pref == "occ": r = min(cand, key=lambda r: (rate[r[0]], r[2]))
                        elif pref == "fifo": r = min(cand, key=lambda r: (r[2], rate[r[0]]))
                        else: r = cand[rng.integers(len(cand))]
                    else:
                        k = min((k for k in range(self.nk) if s in I.elig[k]), key=lambda k: rate[k])
                        r = [k, 0, self.v0[k], []]; robots.append(r)
                    k, i = r[0], r[1]; v = t - ld
                    r[3].append((v, s, i, 0))
                    r[2] = v + I.step(k, i, s)[1]; r[1] = (i + 1) % I.L[k]
                    p[1] = t + I.d[k]
            for k, _, _, duties in robots:
                idx, coef = self.column(k, duties)
                self.h.addCol(self.cost[k], 0.0, highspy.kHighsInf, len(idx), idx, coef)
                self.cols.append((k, duties)); self.ncol_k[k] += 1; added += 1
        return added

    def run_master(self, mode):
        """mode = 'ipm'：内点法、不做 crossover，得到居中的对偶（迭代中用，减少振荡）；
        'exact'：内点法 + crossover，得到基本最优解和精确对偶（收尾复核和最终取值用）；
        'simplex'：原始单纯形热启动。"""
        h = self.h; tm = time.time()
        if mode == "simplex":
            h.setOptionValue("solver", "simplex"); h.setOptionValue("simplex_strategy", 4)
        else:
            h.setOptionValue("solver", "ipm"); h.setOptionValue("run_crossover", "off" if mode == "ipm" else "on")
        h.run(); self.t_master += time.time() - tm
        return h.modelStatusToString(h.getModelStatus())

    def fix_duals(self, y):
        """把耦合对偶 μ_{k,s,D} 换成它在对偶可行域里的最大值 min_j (φ_{j,尾} − φ_{j,头})（头 ≥ H 时 φ 记 0）。
        主问题里耦合对偶只受工具弧约束 μ ≤ φ_尾 − φ_头，没有下界；没被路径用到的键，求解器可能给出任意负值，
        定价时就成了虚假收益。换成上界后：工具弧仍对偶可行，已有路径列的约化成本只增不减，右端为 0 所以对偶目标不变，
        仍是最优对偶。"""
        y = y.copy(); H = self.H
        phi = y[self.tool_row0:].reshape(self.nt, H)
        best = {}
        for (col, k, s, D, j, a, b) in self.tool_arc:
            val = phi[j, a] - (phi[j, b] if b < H else 0.0)
            key = (k, s, D)
            if key not in best or val < best[key]: best[key] = val
        for key, row in self.duty_key.items(): y[row] = best[key]
        return y

    def split_duals(self, y):
        """对偶向量 y → (π: ns×H 覆盖对偶, μ: (k,s) → 长 H+1 的耦合对偶)。"""
        pi = np.zeros((self.ns, self.H))
        for (s, t), row in self.cov_row.items(): pi[s, t] = y[row]
        mu = {(k, s): np.zeros(self.H + 1) for k in range(self.nk) for s in self.I.elig[k]}
        for (k, s, D), row in self.duty_key.items(): mu[(k, s)][D] = y[row]
        return pi, mu

    def price_all(self, y_price, y_true, ncol, tol):
        """用 y_price 定价，每型最多 ncol 条；按真实对偶 y_true 的约化成本 < −tol 才加列。返回 (加列数, 真实对偶下最小约化成本)。"""
        tp = time.time(); h = self.h
        pi, mu = self.split_duals(y_price)
        new = 0; min_rc = 0.0
        for k in range(self.nk):
            pk = pi.copy()
            for c in range(ncol):
                W = self.weights(k, pk, mu)
                val, path = self.price(k, W)
                if not path: break
                idx, coef = self.column(k, path)
                rc = self.cost[k] - float(y_true[idx] @ coef)
                min_rc = min(min_rc, rc)
                if rc < -tol:
                    h.addCol(self.cost[k], 0.0, highspy.kHighsInf, len(idx), idx, coef)
                    self.cols.append((k, path)); self.ncol_k[k] += 1; new += 1
                # 把这条路径覆盖格上的正对偶清零，下一条走别处
                for (v, s, i, r) in path:
                    D = v + self.lead[s]
                    pk[s, D:D + self.I.d[k]] = np.minimum(pk[s, D:D + self.I.d[k]], 0.0)
        self.t_price += time.time() - tp
        return new, min_rc

    def set_delta(self, delta):
        """对偶箱：覆盖行双向松弛的成本 = δ，于是 −δ ≤ π ≤ δ。"""
        self.delta = delta
        idx = np.arange(self.islack, self.nstatic, dtype=np.int32)
        self.h.changeColsCost(len(idx), idx, np.full(len(idx), float(delta)))

    def solve(self, ncol=5, time_limit=6 * 3600, tol=1e-7, log=None, log_every=20, master="ipm", alpha=0.0, delta0=0.05, delta_max=BIG):
        """列生成主循环。
        master：迭代中主问题的解法（'ipm' 居中对偶 / 'simplex'）。alpha：Wentges 平滑系数（0 = 不平滑），
        定价用 α·中心 + (1−α)·当前对偶，中心取上一轮的定价对偶；平滑后找不到负列就退回真实对偶再定价。
        收尾：真实对偶下找不到负列时，用 'exact'（crossover）重解主问题，用精确对偶再定价一次，仍无负列才判最优。"""
        t0 = time.time(); h = self.h
        it = 0; status = "optimal"; hist = []
        self.t_master = 0.0; self.t_price = 0.0; self.n_exact = 0
        self.n_init = self.greedy_columns()
        self.set_delta(delta0); self.delta_hist = []
        st = self.run_master(master); center = None
        while True:
            it += 1
            if st != "Optimal": status = "master_" + st; break
            y = self.fix_duals(np.array(h.getSolution().row_dual)); obj = h.getInfo().objective_function_value
            yp = y if (alpha <= 0 or center is None) else alpha * center + (1 - alpha) * y
            new, min_rc = self.price_all(yp, y, ncol, tol)
            if new == 0 and yp is not y:
                yp = y; new, min_rc = self.price_all(y, y, ncol, tol)
            center = yp
            if new == 0 and master != "simplex" or (new == 0 and self.n_exact == 0):
                # 收尾复核：精确对偶
                st = self.run_master("exact"); self.n_exact += 1
                if st != "Optimal": status = "master_" + st; break
                y = self.fix_duals(np.array(h.getSolution().row_dual)); obj = h.getInfo().objective_function_value
                new, min_rc = self.price_all(y, y, ncol, tol)
                center = None
                if new == 0:
                    # 惩罚问题已最优；松弛非零说明 δ 太小，放大后继续（松弛为 0 时就是原 LP 的最优解）
                    xs = np.array(h.getSolution().col_value)[self.islack:self.nstatic]
                    self.delta_hist.append((self.delta, float(xs.sum()), obj))
                    if xs.sum() > 1e-9 and self.delta < delta_max:
                        self.set_delta(min(delta_max, self.delta * 4)); new = -1
            hist.append((it, round(obj, 6), round(min_rc, 8), new, round(time.time() - t0, 1)))
            if log and (it % log_every == 1 or log_every == 1 or new == 0):
                print(json.dumps(dict(it=it, obj=obj, min_rc=min_rc, new=new, ncols=len(self.cols), sec=round(time.time() - t0, 1),
                                      t_master=round(self.t_master, 1), t_price=round(self.t_price, 1))), file=log, flush=True)
            if new == 0: break
            if time.time() - t0 > time_limit: status = "timeout"; break
            if new < 0: new = 0
            h.setOptionValue("time_limit", max(1.0, time_limit - (time.time() - t0)))
            st = self.run_master(master)
        x = np.array(h.getSolution().col_value)
        obj = h.getInfo().objective_function_value
        slack = float(x[self.islack:self.nstatic].sum())
        if slack > 1e-9 and status == "optimal": status = "slack"
        lam = x[self.nstatic:]
        N = [0.0] * self.nk
        for (k, _), v in zip(self.cols, lam): N[k] += v
        M = x[self.iM:self.iM + self.nt].tolist()
        self.x, self.lam = x, lam
        return dict(obj=obj, status=status, iters=it, ncols=len(self.cols), slack=slack, N=N, M=M,
                    sec=round(time.time() - t0, 1), t_master=round(self.t_master, 1), t_price=round(self.t_price, 1),
                    n_init=self.n_init, n_exact=self.n_exact, delta_hist=self.delta_hist, hist=hist[-5:])

    def save_duals(self, path):
        """保存对偶（槽权重）与 T3 要用的 LP 流：覆盖对偶 π_{s,t}、耦合对偶、N、M，
        按 (k, i, s, D, 是否提前大修) 聚合的值岗流。"""
        y = self.fix_duals(np.array(self.h.getSolution().row_dual)); pi, mu = self.split_duals(y)
        flow = {}
        for (k, path), v in zip(self.cols, self.lam):
            if v <= 1e-12: continue
            for (vv, s, i, r) in path:
                key = (k, i, s, vv + self.lead[s], r)
                flow[key] = flow.get(key, 0.0) + v
        keys = sorted(flow)
        fl = np.array([list(kk) + [flow[kk]] for kk in keys]) if keys else np.zeros((0, 6))
        mu_arr = np.array([[k, s, D, y[row]] for (k, s, D), row in self.duty_key.items()])
        np.savez_compressed(path, pi=pi, cov_first=np.array(self.cov_first), mu=mu_arr, duty_flow=fl,
                            N=np.array([sum(v for (k, _), v in zip(self.cols, self.lam) if k == kk) for kk in range(self.nk)]),
                            M=self.x[self.iM:self.iM + self.nt], H=self.H, reset=int(self.reset))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inst"); ap.add_argument("--mult", type=float, default=4); ap.add_argument("--reset", type=int, default=0)
    ap.add_argument("--ncol", type=int, default=5); ap.add_argument("--time", type=float, default=6 * 3600)
    ap.add_argument("--log_every", type=int, default=20)
    ap.add_argument("--master", default="ipm"); ap.add_argument("--alpha", type=float, default=0.0)
    ap.add_argument("--delta0", type=float, default=0.05)
    a = ap.parse_args()
    from mdrr import load_train
    from export_inst import build
    from freshlb import S_max
    I = load_train(a.inst) if a.inst.startswith("m") else build(a.inst)
    H = int(a.mult * S_max(I))
    cg = FreshCG(I, H, bool(a.reset))
    r = cg.solve(ncol=a.ncol, time_limit=a.time, log=sys.stderr, log_every=a.log_every, master=a.master, alpha=a.alpha, delta0=a.delta0)
    r.update(inst=a.inst, H=H, reset=a.reset)
    print(json.dumps(r))


if __name__ == "__main__":
    main()
