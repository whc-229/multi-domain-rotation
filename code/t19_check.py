"""T19 一致性核对（rost5 / rost5_ref 11-04 的改动：planx 可选第 6 列"记账相位"，规则的历史计数 D、phl 按它记）。
(a) 新路径：C++（新）与 Python 参考实现逐字段比较。计划取自 T19 的真实构造（code/t19_connector.py 的同一套函数，在临时目录里生成），
    规则为 σ 规则的 H = 0 形式（松弛规则 rlp 与线性规则 w + curdef 两种都覆盖），台数为构造台数，或随机 + 0..1；
    一半的组把第 6 列换成随机相位（检验记账在相位不一致时也相同），少数组去掉第 6 列（5 列格式）。
    每组两次运行：3000 格不判环（缺岗、提前大修次数）与 2 万格判环（缺岗、判环结果、周期）。
(b) 旧路径：不用第 6 列时新旧程序 eval 输出逐字段相同（去掉 sec）：
    5 列计划 + PB（RH-MIP 的路径）、松弛规则 rlp、T12 路径（fresh_rlp_w + planf + planW + 随机 RDD 基因型 + resetopt）、
    只算机器人（notools + 5 列计划）、稳态初值（steady.rlp）、线性规则 w + curdef；3000 格不判环与 6 万格判环。
输出 results/T19/check.jsonl，打印不一致数。用法：python3 code/t19_check.py <旧 rost5> <新 rost5>"""
import os, sys, json, math, random, tempfile, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths, U  # noqa: E402
from freshlb import S_max  # noqa: E402
import rost5_ref as ref  # noqa: E402
import t19_connector as T  # noqa: E402
import t13_methods as TM  # noqa: E402
from t11_job import write_geno  # noqa: E402
from t15b_check import rand_plan  # noqa: E402

NEW_CASES = [("steady", "m02-K12-00"), ("steady", "m04-K12-00"), ("steady", "m06-K1-00"), ("best", "m04-K1-00"), ("best", "m02-K12-00"), ("best", "m06-K1-00")]
OLD_INSTS = ["m02-K1-00", "m03-K12-00", "m04-K1-01"]


def run_c(r5, args):
    """运行 rost5，返回 JSON（去掉 sec）。"""
    d = json.loads(subprocess.run([r5] + [str(a) for a in args], capture_output=True, text=True, check=True).stdout)
    d.pop("sec", None); return d


def ev_args(txt, Tm, cert, Nk, Mj, w, opts):
    """eval 的参数：线性规则的权重放第 3 段，其余选项放第 4 段。"""
    return ["eval", txt, Tm, cert] + list(Nk) + ["|"] + list(Mj) + ["|"] + (list(w) if w is not None else []) + ["|"] + list(opts)


def build_plan(kind, name, tmpd):
    """按 t19_connector 的流程在临时目录里构造一个计划（取第一个可用的 θ、first 口径）：返回 (I, txt, spec, 构造, 计划文件, 台数)。"""
    T.OUT = tmpd
    I = load_train(name); S = S_max(I); _, txt = inst_paths(name, "train")
    sp = T.source_spec(kind, name, I, tmpd); Nk, Mj = sp["Nk"], sp["Mj"]
    v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]; c0 = max(v0); lmin = min(I.lead(s) for s in range(len(I.sites)))
    ce = T.cpp_eval(txt, 600000, 1, Nk, Mj, sp["w"], sp["opts"])
    t1 = ce["t_end"] - ce["period"] if ce["cert"] == 1 else 8 * S
    D = math.ceil(S / T.NTH); taus = [t1 + th * D for th in range(T.NTH)]
    tr = dict(snap_ts=taus); ref.run(I, Nk, Mj, taus[-1] + 1, False, trace=tr, **sp["kw"]); hist = T.robot_hist(I, tr)
    for sn in tr["snaps"]:
        sh = sn["t"] - (c0 + lmin); st, _ = T.rollback(I, sn, hist)
        init = [(sn["rk"][r], st[r][0], st[r][1] - sh) for r in range(len(st))]
        trY = {}; y = ref.run(I, Nk, Mj, c0 + 8 * S, False, init=init, trace=trY, **sp["y_kw"])
        if y["gaps"] > 0: continue
        per = T.robot_runs(I, trY["disp"], len(init)); C = T.c_sigma(per, c0)
        cons = T.construct(I, trY["disp"], [e[0] for e in init], [e[1] for e in init], "first", C, c0, S)
        N2 = [Nk[k] + cons["dN"][k] for k in range(len(Nk))]
        return I, txt, sp, cons, N2, Mj
    raise RuntimeError(f"{kind} {name} 没有可用的 θ")


def write_variant(path, cons, I, rng, mode):
    """计划文件的变体：mode = acct（第 6 列 = Y 的相位）、rand（第 6 列随机相位）、five（5 列）。"""
    with open(path, "w") as f:
        f.write(f"end {cons['Tplan']}\n")
        for (t, s, k, i, r, a, n) in cons["rows"]:
            if mode == "five": f.write(f"{s} {t} {k} {i} {r}\n")
            else: f.write(f"{s} {t} {k} {i} {r} {a if mode == 'acct' else rng.randrange(I.L[k])}\n")


def main():
    """(a) 新路径 C++ 与 Python；(b) 旧路径新旧程序。"""
    old, new = sys.argv[1], sys.argv[2]
    rng = random.Random(20261104); out = []; bad = 0
    os.environ["ROST5_BIN"] = new; T.R5 = new
    with tempfile.TemporaryDirectory(prefix="t19chk_") as tmpd:
        # (a) 新路径
        for (kind, name) in NEW_CASES:
            I, txt, sp, cons, N0, M0 = build_plan(kind, name, tmpd)
            for c in range(4):
                mode = ["acct", "rand", "acct", "five"][c] if c < 3 or kind == "steady" else "rand"
                pf = os.path.join(tmpd, f"{kind}_{name}_{c}.planx"); write_variant(pf, cons, I, rng, mode)
                Nk = list(N0) if c == 0 else [x + rng.randint(0, 1) for x in N0]; Mj = list(M0) if c == 0 else [x + rng.randint(0, 1) for x in M0]
                X = ref.parse_planx(pf); X["late"] = True
                opts = sp["y_opts"] + [f"planx={pf}", "planxlate", f"planW={cons['planW']}"]
                diff = []
                for Tm, cert in ((3000, 0), (20000, 1)):
                    cc = run_c(new, ev_args(txt, Tm, cert, Nk, Mj, sp["y_w"], opts))
                    py = ref.run(I, Nk, Mj, Tm, bool(cert), planx=X, plan_until=cons["planW"], **sp["y_kw"])
                    keys = ("gaps", "cert", "period", "resets") if cert else ("gaps", "resets")
                    diff += [(Tm, k, py[k], cc[k]) for k in keys if py[k] != cc[k]]
                bad += bool(diff)
                out.append(dict(part="new_path", kind=kind, inst=name, case=c, plan=mode, rule="w" if sp["y_w"] is not None else "rlp",
                                cpp=dict(gaps=cc["gaps"], cert=cc["cert"], period=cc["period"]), same=int(not diff), diff=diff))
                print("新路径", kind, name, c, mode, "一致" if not diff else diff, flush=True)
        # (b) 旧路径
        pf = os.path.join(tmpd, "old.planx"); gp = os.path.join(tmpd, "g.geno")
        for n in OLD_INSTS:
            I = load_train(n); S = S_max(I); W = 2 * S; _, txt = inst_paths(n, "train")
            pb = run_c(new, ["size", txt, 365, 60000, 60000] + TM.opts("pb", n, S))
            st = json.load(open(os.path.join(ROOT, "results", "T3", "jobs", f"{n}__steady.json")))
            for path in ("rhmip", "rlp", "t12", "notools", "steady", "wcurdef"):
                for c in range(2):
                    base = st if path == "steady" else pb
                    Nk = [x + rng.randint(0, 1) for x in base["Nk"]]; Mj = [x + rng.randint(0, 1) for x in base["Mj"]]
                    w = None
                    if path in ("rhmip", "notools"):
                        rand_plan(rng, I, 4 * S, W, pf); w = TM.W_PB; o = [f"planx={pf}", f"planW={W}"] + (["notools"] if path == "notools" else [])
                        if path == "notools": Mj = []
                    elif path == "rlp": o = [f"rlp={ROOT}/results/T3/rlp/{n}__fresh_noreset_u.rlp"]
                    elif path == "steady": o = [f"rlp={ROOT}/results/T3/rlp/{n}__steady.rlp"]
                    elif path == "wcurdef": w = [round(x + rng.gauss(0, 0.3), 4) for x in TM.W_PB]; o = ["curdef"]
                    else:
                        R6 = os.path.join(ROOT, "results", "T6", "rules")
                        o12 = [f"rlp={R6}/{n}__fresh_rlp_w.rlp", f"planf={R6}/{n}__fresh_rlp_w.planf", f"planW={W}"]
                        write_geno(gp, W, 0, 0.4, [], 1)
                        ndev = run_c(new, ev_args(txt, W, 0, Nk, Mj, None, o12 + [f"rdd={gp}"]))["ndec_w"]
                        write_geno(gp, W, ndev, 0.4, [rng.random() for _ in range(2 * max(1, ndev // 4))], 1)
                        o = o12 + ["resetopt", f"rdd={gp}"]
                    diffs = []
                    for Tm, cert in ((3000, 0), (60000, 1)):
                        eo = run_c(old, ev_args(txt, Tm, cert, Nk, Mj, w, o)); en = run_c(new, ev_args(txt, Tm, cert, Nk, Mj, w, o))
                        if eo != en: diffs.append((Tm, {k: (eo.get(k), en.get(k)) for k in set(eo) | set(en) if eo.get(k) != en.get(k)}))
                    bad += bool(diffs)
                    out.append(dict(part="old_path", inst=n, path=path, case=c, same=int(not diffs), diff=diffs))
                    print("旧路径", n, path, c, "相同" if not diffs else diffs, flush=True)
    os.makedirs(os.path.join(ROOT, "results", "T19"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "T19", "check.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(json.dumps(dict(cases=len(out), new_path=sum(1 for r in out if r["part"] == "new_path"), old_path=sum(1 for r in out if r["part"] == "old_path"), mismatch=bad)))


if __name__ == "__main__":
    main()
