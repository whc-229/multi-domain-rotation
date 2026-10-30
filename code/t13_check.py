"""T13 一致性核对：每个对手的解码路径，rost5（C++）与 Python 参考实现（code/rost5_ref.py）逐字段比较。
对手与随机化：
- fifo、pb：固定权重；台数随机 + 0..1；
- esdw：PB 权重加正态扰动（与 ES 变异同形），curdef；
- eslp：同上（14 维，w13 = LP 引导）+ lpz + 窗口计划（planf，严格窗口 planW）+ curdef；
- gphh：随机 GP 树（深度 ≤ 4，与 rost5 的 grow 同样的算子集），工具权重 W_TOOL；
- rk、d5：窗口内每个决策一个均匀随机基因（d5 带回溯修复），工具权重 W_TOOL；
- rhmip：随机计划（planx，到 4·S_max）+ PB 兜底 + 严格窗口。RH-MIP 的排班通过"按计划派"进入仿真，所以也能与 Python 逐字段比较。
每个 (对手, 实例) 10 组；实例 m02-K12-00（含多岗位站）与 m04-K1-00，共 20 组 / 对手。
每组比较两次运行：Tmax = 3000 不判环（缺岗、提前大修次数、窗口决策数）与 Tcert = 20000 判环（缺岗、判环结果、周期）。
台数：该对手默认（全 0 基因 / 基础权重）用 rost5 size 定出的台数 + 0..1 随机增量；定不出时用实例默认台数 + 0..2。
输出 results/T13/check.jsonl，打印不一致数。用法：python3 code/t13_check.py [rost5 路径]"""
import os, sys, json, random, subprocess, tempfile, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import rost5_ref as ref  # noqa: E402
import t13_methods as TM  # noqa: E402
from t7_check import txt_N0M0  # noqa: E402

R5 = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "code", "cpp", "rost5")
INSTS = ["m02-K12-00", "m04-K1-00"]
OPS = ["+", "-", "*", "/", "max", "min"]; FEATS = ["wait", "nextocc", "flex", "deficit", "duty", "phase"]


def c_json(args):
    """运行 rost5，解析一行 JSON。"""
    r = subprocess.run([R5] + [str(a) for a in args], capture_output=True, text=True, check=True)
    return json.loads(r.stdout)


def rand_tree(rng, d=0, maxd=4):
    """随机 GP 树（前缀式文本），常数用 repr 保证两边读到同一个数。"""
    if d >= maxd or (d > 1 and rng.random() < 0.3):
        return rng.choice(FEATS) if rng.random() < 0.8 else repr(rng.uniform(-2, 2))
    return f"({rng.choice(OPS)} {rand_tree(rng, d + 1, maxd)} {rand_tree(rng, d + 1, maxd)})"


def lpz_of(I, path):
    """读 lpg 文件的 z 矩阵（与 C++ parse_lpz 相同）。"""
    z = [[0.0] * len(I.sites) for _ in range(len(I.T))]
    for line in open(path):
        if line.startswith("z "): _, k, s2, v = line.split(); z[int(k)][int(s2)] = float(v)
    return z


def main():
    rng = random.Random(20261030); out = []; bad = 0
    tmpd = tempfile.mkdtemp(prefix="t13chk_")
    for meth in ["fifo", "pb", "esdw", "eslp", "gphh", "rk", "d5", "rhmip"]:
        for n in INSTS:
            I = load_train(n); S = S_max(I); W = 2 * S; _, txt = inst_paths(n, "train")
            tf = os.path.join(tmpd, "tree.txt")
            # 默认的台数
            if meth in ("fifo", "pb", "rhmip"): base_opts = TM.opts("pb" if meth == "rhmip" else meth, n, S)
            elif meth in ("esdw", "eslp"): base_opts = [TM.wstr(TM.W_PB + ([0] * 5 + [1.0] if meth == "eslp" else []))] + TM.opts(meth, n, S)
            elif meth == "gphh": open(tf, "w").write("(- (* 10 flex) deficit)"); base_opts = TM.opts(meth, n, S, tf)
            else: base_opts = TM.opts(meth, n, S) + [f"pd={meth}", f"W={W}"]
            b = c_json(["size", txt, 365, 20000, 20000] + base_opts)
            if b["R"] >= 0: N0, M0 = b["Nk"], b["Mj"]
            else: N0, M0 = txt_N0M0(I, txt); N0 = [x + rng.randint(0, 2) for x in N0]; M0 = [x + rng.randint(0, 2) for x in M0]
            ndev = c_json(["eval", txt, W, 0] + N0 + ["|"] + M0 + ["|", "|", "pd=rk", f"pdW={W}"])["ndec_w"] if meth in ("rk", "d5") else 0
            for c in range(10):
                Nk = [x + rng.randint(0, 1) for x in N0]; Mj = [x + rng.randint(0, 1) for x in M0]
                pol = {}; ev_w = []; ev_o = []
                if meth in ("fifo", "pb"):
                    w = TM.W_FIFO if meth == "fifo" else TM.W_PB; pol["w"] = [float(x) for x in w] + [0.0] * 6; ev_w = w
                elif meth in ("esdw", "eslp"):
                    nw = 14 if meth == "eslp" else 8
                    w = [x + 0.5 * rng.gauss(0, 1) * (abs(x) + 0.2) for x in (TM.W_PB + [0.0] * (nw - 8))]
                    if meth == "eslp": w[8] = -1.0; w[13] = 1.0 + 0.5 * rng.gauss(0, 1)
                    pol.update(w=w + [0.0] * (14 - nw), curdef=True); ev_w = w; ev_o = ["curdef"]
                    if meth == "eslp":
                        lp = TM.lpg_window(n); pp = f"{TM.RULES6}/{n}__fresh_rlp_w.planf"
                        pol.update(lpz=lpz_of(I, lp), planf=ref.parse_planf(I, pp), plan_until=W); ev_o += [f"lpg={lp}", f"planf={pp}", f"planW={W}"]
                elif meth == "gphh":
                    tr = rand_tree(rng); open(tf, "w").write(tr)
                    pol.update(w=[float(x) for x in TM.W_TOOL] + [0.0] * 6, tree=ref.parse_tree(tr)); ev_w = TM.W_TOOL; ev_o = [f"tree=@{tf}"]
                elif meth == "rhmip":     # 随机计划：每个 (站, 时刻) 以 0.3 的概率放 1–2 条 (可派机型, 随机相位, 窗口内随机提前大修)
                    pf = os.path.join(tmpd, "p.planx"); e = {}
                    with open(pf, "w") as f:
                        f.write(f"end {4 * S}\n")
                        for D in range(4 * S):
                            for s2 in range(len(I.sites)):
                                if rng.random() < 0.3:
                                    for _ in range(rng.randint(1, 2)):
                                        k = rng.choice([k for k in range(len(I.T)) if s2 in I.elig[k]]); i = rng.randrange(I.L[k]); r = int(D < W and rng.random() < 0.2)
                                        f.write(f"{s2} {D} {k} {i} {r}\n")
                    pol.update(w=[float(x) for x in TM.W_PB] + [0.0] * 6, planx=ref.parse_planx(pf), plan_until=W); ev_w = TM.W_PB; ev_o = [f"planx={pf}", f"planW={W}"]
                else:
                    g = [rng.random() for _ in range(ndev)]
                    pol.update(w=[float(x) for x in TM.W_TOOL] + [0.0] * 6, curdef=True, pd=dict(mode=1 if meth == "rk" else 2, W=W, g=g, repair=(meth == "d5")))
                    ev_w = TM.W_TOOL; ev_o = [f"pd={meth}", f"pdW={W}", "pdg=" + ",".join(repr(x) for x in g)]
                rec = dict(method=meth, inst=n, case=c); diff = []
                for Tmax, cert in ((3000, False), (20000, True)):
                    py = ref.run(I, Nk, Mj, Tmax, cert, **pol)
                    cc = c_json(["eval", txt, Tmax, int(cert)] + Nk + ["|"] + Mj + ["|"] + list(ev_w) + ["|"] + ev_o)
                    for k in ("gaps", "cert", "period", "resets", "ndec_w"):
                        if not cert and k in ("cert", "period"): continue
                        if py[k] != cc[k]: diff.append((Tmax, k, py[k], cc[k]))
                    rec[f"c{Tmax}"] = {k: cc[k] for k in ("gaps", "cert", "period", "resets", "ndec_w")}
                rec["same"] = int(not diff); rec["diff"] = diff; bad += bool(diff); out.append(rec)
                print(meth, n, c, "一致" if not diff else f"不一致 {diff}", rec["c20000"], flush=True)
    os.makedirs(os.path.join(ROOT, "results", "T13"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "T13", "check.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(json.dumps(dict(cases=len(out), mismatch=bad)))


if __name__ == "__main__":
    main()
