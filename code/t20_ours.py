"""T20 本方法的预处理（按实例算一次，结果缓存；用时记下来，从每个作业的预算里扣）。只在第五批的测试实例上用。
步骤（与 T3 steady、T19 (a) 相同，只有第 4 步的窗口改了）：
1. 速率 LP（T2 的 rate_job，prior 的 tool_lb）：班率 y[k,s] 与 LP 台数 N、M；
2. 稳态松弛规则 σ 的文件（T3 的 write_rlp，H = 0、steady 1，起点台数 ⌈N⌉、⌈M⌉）；
3. σ 定台数：rost5 size（rlp=σ、rescue），与 T3 的 rlp 命令同一个定台数过程，加 11-05 的补救；
4. 接入构造（code/t19_connector.py 的同一套函数）：8 个 θ × 2 种入场口径，提前大修窗口固定 [0, 2·S_max)
   （11-06 用户定：预工作需要的提前大修落在窗口外时，这个缺口改由额外机器人补；T19 放宽到 T_plan 的做法不再用）；
   取 60 万格可行且资源最少的（同资源取严格证明的，再按 θ、口径）；都不可行时取资源最少的，交给定台数（rescue）补；
5. 规则文件：conn.rlp（σ 的 H = 0 形式，起点台数 = 构造台数）、noconn.rlp（同一规则，起点台数 = σ 的台数，消融"去掉接入构造"用）。
输出 results/T20/<阶段>/pre/<实例>.json 与 results/T20/<阶段>/pre/<实例>/（规则、计划文件）。并发的作业用文件锁等第一个算完。"""
import os, sys, json, math, time, fcntl, subprocess, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, U, load_inst, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import rost5_ref as ref  # noqa: E402
import t2_job  # noqa: E402
import t3_make_rlp  # noqa: E402
import t19_connector as C19  # noqa: E402

R5 = os.environ.get("ROST5_BIN", os.path.join(ROOT, "code", "cpp", "rost5"))     # 环境变量只用于安装前的核对
SPLIT = "test"


def c_json(args):
    """运行 rost5，返回 JSON。"""
    return json.loads(subprocess.run([R5] + [str(a) for a in args], capture_output=True, text=True, check=True).stdout)


def h0_with_fleet(src, dst, Nk=None, Mj=None):
    """把 rlp 规则改写成 H = 0、全新启动（steady 0）的形式（同 t19_connector.h0_rule），可选地把起点台数 N0 / M0 改成给定值。"""
    lines = []
    for ln in open(src):
        a = ln.split()
        if a and a[0] == "H": ln = "H 0\n"
        elif a and a[0] == "steady": ln = "steady 0\n"
        elif a and a[0] in ("F", "FP"): continue
        elif a and a[0] == "N0" and Nk is not None: ln = "N0 " + " ".join(map(str, Nk)) + "\n"
        elif a and a[0] == "M0" and Mj is not None: ln = "M0 " + " ".join(map(str, Mj)) + "\n"
        lines.append(ln)
    tmp = dst + f".{os.getpid()}.tmp"; open(tmp, "w").writelines(lines); os.replace(tmp, dst)
    return dst


def paths(stage, name):
    """预处理的结果文件与产物目录。"""
    d = os.path.join(ROOT, "results", "T20", stage, "pre")
    return os.path.join(d, f"{name}.json"), os.path.join(d, name)


def build(stage, name, split=SPLIT):
    """按顺序做第 1–5 步，返回记录（含用时 wall）。split = train 只用于一致性核对（code/t20_check.py，阶段 check_train）。"""
    t0 = time.time()
    I = load_train(name) if split == "train" else load_inst(name, split); S = S_max(I); _, txt = inst_paths(name, split); W = 2 * S
    _, art = paths(stage, name); os.makedirs(art, exist_ok=True)
    rec = dict(inst=name, m=I.meta["m"], family=I.meta["family"], S=S, W_reset=W)
    # 1–2
    rate = t2_job.rate_job(I)
    y = {}
    for key, v in rate["y"].items():
        k, s = map(int, key.split(",")); y.setdefault(s, {})[k] = v
    ceil = lambda v: int(math.ceil(v - 1e-6))
    rp = os.path.join(art, "steady.rlp")
    t3_make_rlp.write_rlp(rp, I, 0, True, rate["N"], [ceil(v) for v in rate["N"]], [ceil(v) for v in rate["M"]], np.zeros((0, 6)), y)
    rec.update(R_rate=rate["R_rate"], rate_sec=rate["sec"])
    # 3
    sg = c_json(["size", txt, 365, 60000, 600000, f"rlp={rp}", "rescue"])
    sig_feas = sg["R"] >= 0 and ((sg["verify_cert"] == 1 and sg["exact"] == 1) or sg["verify_gaps"] == 0)
    rec.update(sigma_R=sg["R"], sigma_status=sg["status"], sigma_Nk=sg["Nk"], sigma_Mj=sg["Mj"], sigma_feasible=int(sig_feas),
               sigma_factor=sg["start_factor"], sigma_rescue_add=sg.get("rescue_add"), sigma_sec=round(time.time() - t0, 1))
    if sg["R"] < 0:
        rec["wall"] = round(time.time() - t0, 1); return rec
    Nk, Mj = sg["Nk"], sg["Mj"]
    h0 = h0_with_fleet(rp, os.path.join(art, "steady_h0.rlp"))
    h0_noconn = h0_with_fleet(rp, os.path.join(art, "noconn.rlp"), Nk, Mj)
    # 4：与 t19_connector.main 的 steady 分支相同，窗口固定为 2·S_max
    v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]; c0 = max(v0)
    lead_min = min(I.lead(s) for s in range(len(I.sites)))
    kw = dict(rlp=ref.parse_rlp(I, rp)); y_kw = dict(rlp=ref.parse_rlp(I, h0)); y_opts = [f"rlp={h0}"]
    ce = C19.cpp_eval(txt, 600000, 1, Nk, Mj, None, [f"rlp={rp}"])
    t1 = ce["t_end"] - ce["period"] if ce["cert"] == 1 else 8 * S
    D = math.ceil(S / C19.NTH); taus = [t1 + th * D for th in range(C19.NTH)]; Tsrc = taus[-1] + 1
    tr = dict(snap_ts=taus)
    py = ref.run(I, Nk, Mj, Tsrc, False, trace=tr, **kw)
    cc = C19.cpp_eval(txt, Tsrc, 0, Nk, Mj, None, [f"rlp={rp}"])
    rec.update(t1=t1, t1_from_cert=int(ce["cert"] == 1), src_repro=int(py["gaps"] == cc["gaps"] and py["resets"] == cc["resets"]))
    hist = C19.robot_hist(I, tr); variants = []
    for th, sn in enumerate(tr["snaps"]):
        sh = sn["t"] - (c0 + lead_min)
        st, _ = C19.rollback(I, sn, hist)
        init = [(sn["rk"][r], st[r][0], st[r][1] - sh) for r in range(len(sn["rk"]))]
        trY = {}
        yv = ref.run(I, Nk, Mj, c0 + 8 * S, False, init=init, trace=trY, **y_kw)
        if yv["gaps"] > 0: continue
        per = C19.robot_runs(I, trY["disp"], len(init)); Cs = C19.c_sigma(per, c0)
        rph0 = [e[1] for e in init]
        for mode in ("last", "first"):
            cons = C19.construct(I, trY["disp"], [e[0] for e in init], rph0, mode, Cs, c0, S, win=W)
            pf = os.path.join(art, f"th{th}__{mode}.planx"); C19.write_planx(pf, cons["Tplan"], cons["rows"])
            N2 = [Nk[k] + cons["dN"][k] for k in range(len(Nk))]
            opts = y_opts + [f"planx={pf}", "planxlate", f"planW={W}"]
            ev = C19.cpp_eval(txt, 600000, 1, N2, Mj, None, opts)
            variants.append(dict(theta=th, mode=mode, plan=pf, N=N2, R=sum(N2) + sum(Mj), dN=cons["dN"], Tplan=cons["Tplan"], n_ext=cons["n_ext"],
                                 n_pre=cons["n_pre"], n_reset=cons["n_reset"], feasible=int(ev["gaps"] == 0), proven=int(ev["cert"] == 1)))
    rec["n_variants"] = len(variants); rec["n_feasible"] = sum(v["feasible"] for v in variants)
    if variants:
        ok = [v for v in variants if v["feasible"]] or variants
        b = min(ok, key=lambda v: (v["R"], -v["proven"], v["theta"], v["mode"]))
        conn = h0_with_fleet(rp, os.path.join(art, "conn.rlp"), b["N"], Mj)
        rec["conn"] = dict(b, plan=os.path.relpath(b["plan"], ROOT), rule=os.path.relpath(conn, ROOT))
    rec["thetas"] = [{k: v[k] for k in ("theta", "mode", "R", "feasible", "proven", "Tplan")} for v in variants]
    rec["wall"] = round(time.time() - t0, 1)
    return rec


def prepare(stage, name, split=SPLIT):
    """取缓存的预处理结果；没有就算（加锁，并发的作业等第一个算完）。"""
    jp, _ = paths(stage, name); os.makedirs(os.path.dirname(jp), exist_ok=True)
    if os.path.exists(jp): return json.load(open(jp))
    with open(jp + ".lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        if os.path.exists(jp): return json.load(open(jp))
        rec = build(stage, name, split)
        tmp = jp + ".tmp"
        with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        os.replace(tmp, jp)
    return rec


def spec(stage, name, I, S, method):
    """方法的 rost5 规则选项、eval 选项、Python 参数。ours / ours_w* / conn_only：接入构造；noconn：σ 规则从全新启动。"""
    jp, art = paths(stage, name); pre = json.load(open(jp)); W = 2 * S
    if method == "noconn" or "conn" not in pre:
        rule = os.path.join(art, "noconn.rlp")
        opts = [f"rlp={rule}", f"planW={W}"]
        return dict(opts=opts, eval_opts=opts, py=dict(rlp=ref.parse_rlp(I, rule), plan_until=W), fallback=int(method != "noconn"))
    rule = os.path.join(ROOT, pre["conn"]["rule"]); plan = os.path.join(ROOT, pre["conn"]["plan"])
    opts = [f"rlp={rule}", f"planx={plan}", "planxlate", f"planW={W}", "devfirst"]
    X = ref.parse_planx(plan); X["late"] = True
    return dict(opts=opts, eval_opts=opts, py=dict(rlp=ref.parse_rlp(I, rule), planx=X, plan_until=W, devfirst=True), fallback=0)


if __name__ == "__main__":
    print(json.dumps(prepare(sys.argv[1], sys.argv[2], *sys.argv[3:4]), ensure_ascii=False))
