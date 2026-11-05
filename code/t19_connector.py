"""T19 接入构造（11-04 审核：只做 (a) 与 (c)）。只用训练实例。
按接入定理构造启动期：稳态部分照抄一个轮转，启动期按"照抄接入 + 入场前预工作补缺口"构造，[0, T₀) 内按计划派（planx），之后按稳态规则。

σ 的来源（两种实验共用同一条流程）：
- (a) kind = steady：T3 的稳态松弛规则 steady.rlp（H = 0），台数取 T3 定出的 N_σ、M_σ（稳态初值下定台数）；
- (c) kind = best：已知最好 RDD 解（T8 / T11 / T16（β = 0）/ T18 中资源最少的可行解），规则取它的规则在 t ≥ H 时的形式（rlp 文件改成 H = 0；线性规则不变）。
流程（每个 σ 来源）：
1. C++ eval（60 万格判环）取该解进入循环的时刻 t₁；Python 参考实现复现到 t₁ + 8·Δ（Δ = ⌈S_max/8⌉），
   核对 Python 与 C++ 的缺岗、提前大修次数相同；在 τ_θ = t₁ + θ·Δ（θ = 0..7）记下机器人状态（σ 在一个维修循环里均匀取 8 个时间原点）。
2. 参照排班 Y_θ：用 σ 的规则（H = 0 形式），机器人初值取 τ_θ 时的状态（型、相位，就绪基准平移到全新启动的时间原点），
   岗位、工具、规则计数按全新启动，在 Python 里跑到 c₀ + 8·S_max。Y_θ 有缺岗的 θ 不用。
   技术选择：岗位在全新启动时同时开放（c₀ + lead_s），σ 本身的岗位在任一时刻都有剩余班长，无法直接对齐；
   所以用"σ 的机器人分布 + 全新岗位"跑出的 Y_θ 作接入目标。Y_θ 的机器人分布与 σ 相同，规则相同。
   τ_θ 时正在值岗的机器人退回到这一班派工之前的状态（rollback），让它们去值全新岗位开头的班；不退回的话开头可用的机器人不够。
3. C_σ：Y_θ 中任一机器人相邻两次回到相位 0 的最长间隔（[c₀, c₀ + 8·S_max) 内）。
4. 接入（两种入场口径，都算）：
   - last（任务书）：T₀ = c₀ + C_σ，机器人 r 的新机器人 r′ 在 b_r（r 在 T₀ 前最后一次回到相位 0 的时刻）之后、r 的下一班入场；
   - first：r′ 在 r 第一次回到相位 0 之后入场（初始就在相位 0 的，从第一班起照抄）。
   r 在入场之前的班就是缺口。按时间逐个指派（贪心）：
     先给已有的额外机器人（就绪最晚的）；再给尚未入场、此刻就绪、值完后提前大修（或正好在相位 L−1 自然大修）
     还能在入场前回到相位 0 并就绪的新机器人（剩余余量最小的）；都不行就加 1 台额外机器人。
   每个缺口都由同型机器人补，所以每格每站的 (型, 时刻) 与 Y_θ 完全相同：值岗时长、岗位时间线、工具用法都不变，ΔM = 0。
5. 计划文件 planx：[0, T_plan) 内 Y_θ 的每一班一条（T_plan = 所有机器人都入场的时刻），执行者的 (型, 相位, 是否本班后提前大修)，
   加第 6 列"记账相位" = 这一班在 Y_θ 里的相位：规则的历史计数（D、phl）按它记，这样 T_plan 时规则计数与 Y_θ 相同
   （rost5 / rost5_ref 11-04 加的可选列；不加时预工作的相位会让 phl 与 Y_θ 差很多，规则接手后马上走偏）。之后按 σ 的规则。C++ eval：台数 N_σ + ΔN、M_σ，选项 = 规则 + planx + planxlate + planW = max(2·S_max, T_plan)；
   计划期缺岗（应为 0）与 60 万格判环。可行 = 严格证明或 60 万格 0 缺岗。
结果：results/T19/<kind>/<实例>.json（原子写入，已存在就跳过）；计划文件与改写的规则文件放 results/T19/plans、results/T19/rules。
用法：python3 code/t19_connector.py <steady|best> <实例>"""
import os, sys, json, math, subprocess, tempfile, time, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths, U  # noqa: E402
from freshlb import S_max  # noqa: E402
import rost5_ref as ref  # noqa: E402
from t11_job import write_geno  # noqa: E402

R5 = os.environ.get("ROST5_BIN", os.path.join(ROOT, "code", "cpp", "rost5"))     # 环境变量只用于开发期试算
OUT = os.path.join(ROOT, "results", "T19")
NTH = 8                                   # 一个维修循环里取的时间原点个数
SEARCH_OPTS = ("stall=", "swap=", "fdsel=", "fddump=", "hsfb=")
BEST_DIRS = ["T8", "T11", "T16", "T18"]   # (c) 的候选来源（RDD 解）


def c_json(args):
    """运行 rost5，返回 JSON。"""
    return json.loads(subprocess.run([R5] + [str(a) for a in args], capture_output=True, text=True, check=True).stdout)


def h0_rule(src, tag):
    """把 rlp 规则文件改写成 H = 0、全新启动（steady 0）的形式，写到 results/T19/rules/，返回路径。
    H = 0 时规则只用末段速率 p 与计数（与原规则在 t ≥ H 时的动力学相同）；表 F、FP 不再被用到，去掉（C++ 按 H 读表长）。"""
    d = os.path.join(OUT, "rules"); os.makedirs(d, exist_ok=True)
    dst = os.path.join(d, os.path.basename(src).replace(".rlp", f"__{tag}.rlp"))
    if not os.path.exists(dst):
        lines = []
        for ln in open(src):
            a = ln.split()
            if a and a[0] == "H": ln = "H 0\n"
            elif a and a[0] == "steady": ln = "steady 0\n"
            elif a and a[0] in ("F", "FP"): continue          # C++ 按 H 读表长，H = 0 时去掉（本来也用不到）
            lines.append(ln)
        tmp = dst + f".{os.getpid()}.tmp"
        open(tmp, "w").writelines(lines); os.replace(tmp, dst)
    return dst


def feas_rec(d):
    """作业记录可行时返回资源数，否则 None。"""
    if not d or d.get("R") is None or d["R"] < 0: return None
    return d["R"] if ((d.get("verify_cert") == 1 and d.get("exact") == 1) or d.get("verify_gaps") == 0) else None


def parse_cmd(rec):
    """从作业记录的命令行取规则选项：返回 (w 或 None, 其余规则选项, 是否 resetopt)。去掉搜索参数与 facloop。"""
    toks = rec["cmd"].split()
    w = None; opts = []; reset = False
    for t in toks:
        if "=" not in t and t not in ("curdef", "resetopt"): continue
        if t.startswith("w="): w = [float(x) for x in t[2:].split(",")]
        elif t == "resetopt": reset = True
        elif t.startswith("W=") or t.startswith(SEARCH_OPTS): continue
        elif t.startswith(("rlp=", "planf=", "planW=")) or t == "curdef": opts.append(t)
        else: raise ValueError(f"不认识的选项 {t}（{rec.get('inst')}）")
    return w, opts, reset


def py_kw(I, w, opts):
    """规则选项 → rost5_ref.run 的参数（rlp / planf / plan_until / w / curdef）。"""
    kw = {}
    for o in opts:
        if o.startswith("rlp="): kw["rlp"] = ref.parse_rlp(I, o[4:])
        elif o.startswith("planf="): kw["planf"] = ref.parse_planf(I, o[6:])
        elif o.startswith("planW="): kw["plan_until"] = int(o[6:])
        elif o == "curdef": kw["curdef"] = True
    if w is not None: kw["w"] = list(w) + [0.0] * 6
    return kw


def cpp_eval(txt, T, cert, Nk, Mj, w, opts):
    """rost5 eval：线性规则的权重放第 3 段，其余选项放第 4 段。"""
    return c_json(["eval", txt, T, cert] + list(Nk) + ["|"] + list(Mj) + ["|"] + (list(w) if w is not None else []) + ["|"] + list(opts))


def best_record(name):
    """(c)：T8 / T11 / T16 / T18 里资源最少的可行 RDD 解（并列时取严格证明的，再按目录与文件名顺序）。"""
    cands = []
    for di, dname in enumerate(BEST_DIRS):
        for f in sorted(glob.glob(os.path.join(ROOT, "results", dname, "jobs", f"{name}__*.json"))):
            d = json.load(open(f)); R = feas_rec(d)
            if R is None or not d.get("genes") or d.get("cdf"): continue      # T16 β > 0 的偏离位置按对偶分布解码，不用
            proven = int(d.get("verify_cert") == 1 and d.get("exact") == 1)
            cands.append((R, -proven, di, f, d))
    if not cands: return None
    cands.sort(key=lambda c: c[:4])
    return cands[0][3], cands[0][4]


def source_spec(kind, name, I, tmpd):
    """σ 来源的设定：台数、C++ 选项（含 rdd 基因型文件）、Python 参数、σ 规则的 H = 0 形式（C++ 选项与 Python 参数）。"""
    if kind == "steady":
        rp = os.path.join(ROOT, "results", "T3", "rlp", f"{name}__steady.rlp")
        jb = json.load(open(os.path.join(ROOT, "results", "T3", "jobs", f"{name}__steady.json")))
        src = dict(file=os.path.relpath(rp, ROOT), job=f"results/T3/jobs/{name}__steady.json", R=jb["R"], status=jb["status"])
        h0 = h0_rule(rp, "h0fresh")
        return dict(src=src, Nk=jb["Nk"], Mj=jb["Mj"], w=None, opts=[f"rlp={rp}"], kw=dict(rlp=ref.parse_rlp(I, rp)),
                    y_w=None, y_opts=[f"rlp={h0}"], y_kw=dict(rlp=ref.parse_rlp(I, h0)), feasible=feas_rec(jb) is not None)
    f, rec = best_record(name)
    w, opts, reset = parse_cmd(rec)
    gp = os.path.join(tmpd, "best.geno"); write_geno(gp, rec["W"], rec["ndev"], rec["qon"], rec["genes"], int(reset))
    kw = py_kw(I, w, opts)
    kw["rdd"] = ref.decode_geno(rec["genes"], rec["W"], rec["ndev"], rec["qon"], reset)
    y_opts = [o for o in opts if not o.startswith(("planf=", "planW="))]
    y_opts = [f"rlp={h0_rule(o[4:], 'h0')}" if o.startswith("rlp=") else o for o in y_opts]
    src = dict(job=os.path.relpath(f, ROOT), R=rec["R"], status=rec["status"], proven=int(rec.get("verify_cert") == 1 and rec.get("exact") == 1),
               resets=rec.get("resets"))
    return dict(src=src, Nk=rec["Nk"], Mj=rec["Mj"], w=w, opts=opts + [f"W={rec['W']}"] + (["resetopt"] if reset else []) + [f"rdd={gp}"], kw=kw,
                y_w=w, y_opts=y_opts, y_kw=py_kw(I, w, y_opts), feasible=True)


def robot_runs(I, disp, nr):
    """Y 的派工表 → 每台机器人的班 [(t, 站, 型, 相位, 是否回到相位 0, 派工序号)]（按时间）。"""
    per = [[] for _ in range(nr)]
    for n, e in enumerate(disp):
        t, s, q, r, k, i, planned, rvn, end = e[:9]; early = e[11]
        per[r].append((t, s, k, i, int(i == I.L[k] - 1 or early == 1), n))
    return per


def c_sigma(per, c0):
    """C_σ：任一机器人相邻两次回到相位 0 的最长间隔（只看已观察到的两次之间）。"""
    C = 0
    for lst in per:
        ret = [e[0] for e in lst if e[4]]
        for a, b in zip(ret, ret[1:]): C = max(C, b - a)
    return C


def entries(per, rph0, mode, T0):
    """每台 Y 机器人的入场：返回 [(入场班的派工序号或 None, 入场时刻 e_r（无班为 None）, 入场站)]，以及缺口派工序号集合。
    first：初始在相位 0 → 从第一班起照抄；否则第一次回到相位 0 后的下一班入场。
    last：T₀ 之前最后一次回到相位 0 后的下一班入场；T₀ 之前没有回到相位 0 的，按 first。"""
    ent = []; gaps = set()
    for r, lst in enumerate(per):
        if not lst: ent.append((None, None, None)); continue
        ret_idx = [j for j, e in enumerate(lst) if e[4]]
        if mode == "last":
            before = [j for j in ret_idx if lst[j][0] < T0]
            if before: j0 = before[-1] + 1
            elif rph0[r] == 0: j0 = 0
            else: j0 = ret_idx[0] + 1 if ret_idx else len(lst)
        else:
            j0 = 0 if rph0[r] == 0 else (ret_idx[0] + 1 if ret_idx else len(lst))
        for j in range(j0): gaps.add(lst[j][5])
        if j0 < len(lst): ent.append((lst[j0][5], lst[j0][0], lst[j0][1]))
        else: ent.append((None, None, None))
    return ent, gaps


def construct(I, disp, rk, rph0, mode, C, c0, S, win=None):
    """贪心接入：返回 dict（T₀、计划条目、各型额外台数、缺口数、各类执行者计数）。计划条目 (t, 站, 型, 相位, 提前大修, 派工序号)。
    win（T20）：提前大修窗口。None = T19 的做法（窗口放宽到 max(2·S_max, T_plan)）；给定时窗口固定为 [0, win)，
    预工作需要的提前大修落在窗口外时这个新机器人不能用，缺口改由额外机器人补（自然大修不受窗口限制）。"""
    nr = len(rk); nk = len(I.T)
    per = robot_runs(I, disp, nr)
    T0 = c0 + C
    ent, gaps = entries(per, rph0, mode, T0)
    e_times = [e[1] for e in ent if e[1] is not None]
    Tplan = max(e_times) + 1 if e_times else c0
    if mode == "last": Tplan = max(Tplan, T0)
    planW = max(2 * S, Tplan) if win is None else win
    v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]
    lead = [I.lead(s) for s in range(len(I.sites))]
    # 新机器人（每台 Y 机器人一台）：入场前的预工作状态
    fresh = [dict(k=rk[r], ph=0, rv=v0[rk[r]], e=ent[r][1], se=ent[r][2], last=None) for r in range(nr)]
    extras = []                                        # 额外机器人：dict(k, ph, rv)
    plan = {}                                          # 派工序号 → [t, s, k, ph, reset]
    who = {}                                           # 派工序号 → ("copy" | "pre" | "extra", 编号)
    for r in range(nr):
        for e in per[r]:
            if e[5] not in gaps and e[0] < Tplan:
                plan[e[5]] = [e[0], e[1], e[2], e[3], 0]; who[e[5]] = ("copy", r)
    n_pre = n_ext = 0
    for n in sorted(gaps, key=lambda n: (disp[n][0], n)):
        t, s, k = disp[n][0], disp[n][1], disp[n][4]
        L = I.L[k]
        # 1. 已有的额外机器人（就绪最晚的）
        cand = [x for x in extras if x["k"] == k and x["rv"] + lead[s] <= t]
        if cand:
            x = max(cand, key=lambda x: x["rv"]); ld, nx = I.step(k, x["ph"], s)
            plan[n] = [t, s, k, x["ph"], 0]; who[n] = ("extra", extras.index(x))
            x["rv"] = t - ld + nx; x["ph"] = (x["ph"] + 1) % L; n_ext += 1
            continue
        # 2. 尚未入场的新机器人：值完后（提前大修或自然大修）还能在入场前回到相位 0 并就绪；取余量最小的
        best = None
        for fi, F in enumerate(fresh):
            if F["k"] != k or F["e"] is None or F["e"] <= t or F["rv"] + lead[s] > t: continue
            ld, nx = I.step(k, F["ph"], s)
            if F["ph"] == L - 1: rv_ret = t - ld + nx
            elif t < planW: rv_ret = t - ld + I.step(k, L - 1, s)[1]
            else: continue
            slack = F["e"] - (rv_ret + lead[F["se"]])
            if slack < 0: continue
            if best is None or slack < best[0]: best = (slack, fi)
        if best is not None:
            F = fresh[best[1]]; ld, nx = I.step(k, F["ph"], s)
            if F["last"] is not None: plan[F["last"]][4] = 0          # 前一班不再提前大修
            plan[n] = [t, s, k, F["ph"], int(F["ph"] != L - 1)]; who[n] = ("pre", best[1])
            F["last"] = n if F["ph"] != L - 1 else None
            F["rv"] = t - ld + nx; F["ph"] = (F["ph"] + 1) % L; n_pre += 1
            continue
        # 3. 加 1 台额外机器人
        x = dict(k=k, ph=0, rv=v0[k]); extras.append(x)
        ld, nx = I.step(k, 0, s)
        plan[n] = [t, s, k, 0, 0]; who[n] = ("extra", len(extras) - 1)
        x["rv"] = t - ld + nx; x["ph"] = 1 % L; n_ext += 1
    dN = [sum(1 for x in extras if x["k"] == k) for k in range(nk)]
    n_reset = sum(1 for v in plan.values() if v[4])
    rows = sorted(((v[0], v[1], v[2], v[3], v[4], disp[n][5], n) for n, v in plan.items()), key=lambda e: (e[0], e[6]))
    return dict(mode=mode, T0=T0, Tplan=Tplan, planW=planW, dN=dN, n_gap=len(gaps), n_pre=n_pre, n_ext=n_ext, n_reset=n_reset,
                n_plan=len(rows), n_noentry=sum(1 for e in ent if e[1] is None), rows=rows)


def robot_hist(I, tr):
    """σ 来源运行的逐台派工史：[(派工时刻, 值岗结束时刻, 派工前相位, 派工前就绪基准)]（按时间）。派工前的 rv 由上一班的 rv 或初值得到。"""
    rv = list(tr["rv0"]); hist = [[] for _ in rv]
    for e in tr["disp"]:
        t, s, q, r, k, i, planned, rvn, end = e[:9]
        hist[r].append((t, end, i, rv[r])); rv[r] = rvn
    return hist


def rollback(I, sn, hist):
    """快照时刻正在值岗的机器人（派工时刻 < 快照时刻 < 值岗结束）退回到这一班派工之前的状态（相位、rv）。
    全新启动时所有岗位同时开放，退回后这些机器人正好去值开头的班。返回 ([(相位, rv)], 退回台数)。"""
    out = []; nb = 0
    for r in range(len(sn["rk"])):
        cur = (sn["rph"][r], sn["rv"][r])
        for (t, end, i, rvp) in reversed(hist[r]):
            if t >= sn["t"]: continue
            if end > sn["t"]: cur = (i, rvp); nb += 1
            break
        out.append(cur)
    return out, nb


def write_planx(path, Tplan, rows):
    """计划文件：首行 end，之后每行 "s D k i r a"（a = 记账相位：这一班在 Y 里的相位），同一 (站, 时刻) 内按 Y 的派工顺序。"""
    tmp = path + f".{os.getpid()}.tmp"
    with open(tmp, "w") as f:
        f.write(f"end {Tplan}\n")
        for (t, s, k, i, r, a, n) in rows: f.write(f"{s} {t} {k} {i} {r} {a}\n")
    os.replace(tmp, path)


def main():
    """一个 (来源, 实例)：复现 σ、取 8 个时间原点、跑 Y_θ、两种入场口径构造、C++ 判定，写 json。"""
    kind, name = sys.argv[1], sys.argv[2]
    assert kind in ("steady", "best")
    od = os.path.join(OUT, kind); os.makedirs(od, exist_ok=True); pd_ = os.path.join(OUT, "plans"); os.makedirs(pd_, exist_ok=True)
    path = os.path.join(od, f"{name}.json")
    if os.path.exists(path): print("已完成，跳过", kind, name); return
    t_start = time.time()
    I = load_train(name); S = S_max(I); _, txt = inst_paths(name, "train")
    v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]; c0 = max(v0)
    lead_min = min(I.lead(s) for s in range(len(I.sites)))
    rec = dict(inst=name, m=I.meta["m"], family=I.meta["family"], kind=kind, S=S, c0=c0)
    with tempfile.TemporaryDirectory(prefix=f"t19_{kind}_{name}_") as tmpd:
        sp = source_spec(kind, name, I, tmpd)
        Nk, Mj = sp["Nk"], sp["Mj"]
        rec.update(src=sp["src"], N_sigma=Nk, M_sigma=Mj, R_sigma=sum(Nk) + sum(Mj))
        # 1. σ 进入循环的时刻（C++）
        ce = cpp_eval(txt, 600000, 1, Nk, Mj, sp["w"], sp["opts"])
        rec["src_eval"] = {k: ce.get(k) for k in ("gaps", "cert", "period", "t_end", "resets")}
        if ce["cert"] == 1: t1 = ce["t_end"] - ce["period"]
        else: t1 = 8 * S                                 # 判不出环时取 8·S_max 当作已进入稳态（写进记录）
        rec["t1"] = t1; rec["t1_from_cert"] = int(ce["cert"] == 1)
        D = math.ceil(S / NTH); taus = [t1 + th * D for th in range(NTH)]
        Tsrc = taus[-1] + 1
        tr = dict(snap_ts=taus)
        py = ref.run(I, Nk, Mj, Tsrc, False, trace=tr, **sp["kw"])
        hist = robot_hist(I, tr)
        cc = cpp_eval(txt, Tsrc, 0, Nk, Mj, sp["w"], sp["opts"])
        rec["src_repro"] = dict(T=Tsrc, py_gaps=py["gaps"], cpp_gaps=cc["gaps"], py_resets=py["resets"], cpp_resets=cc["resets"],
                                same=int(py["gaps"] == cc["gaps"] and py["resets"] == cc["resets"]))
        thetas = []
        for th, sn in enumerate(tr["snaps"]):
            sh = sn["t"] - (c0 + lead_min)
            st, n_back = rollback(I, sn, hist)
            init = [(sn["rk"][r], st[r][0], st[r][1] - sh) for r in range(len(sn["rk"]))]
            TY = c0 + 8 * S
            trY = {}
            y = ref.run(I, Nk, Mj, TY, False, init=init, trace=trY, **sp["y_kw"])
            row = dict(theta=th, tau=sn["at"], t_snap=sn["t"], n_rollback=n_back, y_gaps=y["gaps"], y_resets=y["resets"])
            if y["gaps"] > 0: thetas.append(row); continue
            per = robot_runs(I, trY["disp"], len(init))
            C = c_sigma(per, c0); row["C"] = C
            rph0 = [e[1] for e in init]
            row["phase0_init"] = sum(1 for p in rph0 if p == 0)
            row["variants"] = []
            for mode in ("last", "first"):
                cons = construct(I, trY["disp"], [e[0] for e in init], rph0, mode, C, c0, S)
                pf = os.path.join(pd_, f"{kind}__{name}__th{th}__{mode}.planx")
                write_planx(pf, cons["Tplan"], cons["rows"])
                N2 = [Nk[k] + cons["dN"][k] for k in range(len(Nk))]
                opts = sp["y_opts"] + [f"planx={pf}", "planxlate", f"planW={cons['planW']}"]
                g_plan = cpp_eval(txt, cons["Tplan"], 0, N2, Mj, sp["y_w"], opts)["gaps"]
                ev = cpp_eval(txt, 600000, 1, N2, Mj, sp["y_w"], opts)
                v = {k: cons[k] for k in cons if k != "rows"}
                v.update(plan_file=os.path.relpath(pf, ROOT), N=N2, R=sum(N2) + sum(Mj), plan_gaps=g_plan,
                         eval={k: ev.get(k) for k in ("gaps", "cert", "period", "t_end", "resets")},
                         feasible=int(ev["gaps"] == 0), proven=int(ev["cert"] == 1))
                row["variants"].append(v)
            thetas.append(row)
        rec["thetas"] = thetas
        ok = [(v["R"], -v["proven"], t["theta"], v["mode"], v) for t in thetas for v in t.get("variants", []) if v["feasible"]]
        rec["best"] = None
        if ok:
            ok.sort(key=lambda e: e[:4]); b = ok[0][4]
            rec["best"] = dict(theta=ok[0][2], mode=b["mode"], R=b["R"], dN=b["dN"], N=b["N"], T0=b["T0"], Tplan=b["Tplan"], proven=b["proven"],
                               T0_over_S=round(b["T0"] / S, 3), Tplan_over_S=round(b["Tplan"] / S, 3))
    rec["wall"] = round(time.time() - t_start, 1)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: rec[k] for k in ("inst", "kind", "R_sigma", "best", "wall")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
