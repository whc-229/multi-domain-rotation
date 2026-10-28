"""T7 一致性核对：rost5 的 RDD 解码 vs Python 参考实现（code/rost5_ref.py）。
对每个 (实例, 默认规则, 随机基因型)：
  (e) 仿真层：同一组台数（默认规则的台数 + 0..1 的随机增量），Tmax = 3000 不判环、Tcert = 60000 判环两次运行，
      比较缺岗数、判环结果、周期、提前大修次数、窗口内决策数；
  (s) 定台数层：rost5 size（规则 + 基因型，新补台规则）与 Python size_rule 比较资源数、各型台数、各工具数、状态。
另外核对"全部槽关闭 = 默认规则"：rost5 size 带全关基因型与不带基因型，输出除用时外逐字段相同；eval 同理。
基因型：K = ⌈ρ·N_dev / q_on⌉，ρ 从 {0.02, 0.1, 0.4} 随机取，q_on = 0.4，(p, a) 均匀随机；一半的基因型带"本班后提前大修"备选项。
N_dev 取默认规则在其台数下窗口 [0, 2·S_max) 内的决策数（rost5 eval 输出），Python 的同一数值也参与比较。
结果逐行写进 results/T7/check_rdd.jsonl（每个进程写自己的分片，最后合并）。
第二轮（证书阶段补台加了"连续 5 次未刷新最低缺岗即停"之后）写到 check_rdd_v2.jsonl：第 4 个参数给输出文件名。
用法：python3 code/t7_check.py <实例,...> <每个 (实例, 规则) 的基因型数> <进程数> [输出文件名]"""
import os, sys, json, random, subprocess, math
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
import rost5_ref as ref  # noqa: E402
import t8_rules  # noqa: E402
from freshlb import S_max  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")
OUT = os.path.join(ROOT, "results", "T7")
HS, TC = 365, 60000


def c_json(args):
    """运行 rost5，解析一行 JSON。"""
    r = subprocess.run([R5] + [str(a) for a in args], capture_output=True, text=True)
    return json.loads(r.stdout)


def txt_N0M0(I, txt):
    """实例文本文件末尾的默认台数 N0、M0（C++ Inst::load 读到的同一组数）。"""
    text = open(txt).read(); nk, ns, nt = map(int, text.split("\n")[1].split()[:3]); tok = text.split()
    tail = [int(x) for x in tok[len(tok) - nk - nt:]]
    return tail[:nk], tail[nk:]


def py_rule(I, name, rule, txt):
    """把规则选项翻译成 Python 参考实现的参数：返回 (pol, N0, M0, is_rlp, facloop)。"""
    opts = t8_rules.rule_opts(name, rule)
    pol = {}; N0 = M0 = None; is_rlp = facloop = False
    for o in opts:
        if o.startswith("rlp="):
            R = ref.parse_rlp(I, o[4:]); pol["rlp"] = R; N0, M0 = R["N0"], R["M0"]; is_rlp = True
        elif o.startswith("w="):
            w = [float(x) for x in o[2:].split(",")]; pol["w"] = w + [0.0] * (14 - len(w))
        elif o.startswith("lpg="):
            nk, ns = len(I.T), len(I.sites); z = [[0.0] * ns for _ in range(nk)]
            for line in open(o[4:]):
                if line.startswith("#N0"): N0 = [int(x) for x in line[3:].split()]
                elif line.startswith("#M0"): M0 = [int(x) for x in line[3:].split()]
                elif line.startswith("z "): _, k, s2, v = line.split(); z[int(k)][int(s2)] = float(v)
            pol["lpz"] = z
        elif o == "curdef": pol["curdef"] = True
        elif o == "facloop": facloop = True
    if N0 is None: N0, M0 = txt_N0M0(I, txt)
    return opts, pol, N0, M0, is_rlp, facloop


def write_geno(path, g, W, ndev, qon, ro):
    """写基因型文件（%.17g，保证 C++ 读回的数与 Python 相同）。"""
    with open(path, "w") as f:
        f.write(f"W {W}\nndev {ndev}\nqon {qon!r}\nresetopt {int(ro)}\n")
        for a in range(0, len(g), 2): f.write(f"slot {g[a]!r} {g[a + 1]!r}\n")


def one(case):
    """一个 (实例, 规则, 基因型序号) 的全部比较；返回记录。"""
    name, rule, gi, seed = case
    I = load_train(name); _, txt = inst_paths(name, "train")
    opts, pol, N0, M0, is_rlp, facloop = py_rule(I, name, rule, txt)
    W = 2 * S_max(I)
    base = c_json(["size", txt, HS, TC, 600000] + opts)
    if base["R"] < 0: return dict(inst=name, rule=rule, g=gi, skip="默认规则定台数失败")
    gdir = os.path.join(OUT, "geno"); os.makedirs(gdir, exist_ok=True)
    off = os.path.join(gdir, f"{name}__{rule}__off_g{gi}.txt"); write_geno(off, [0.0, 0.0], W, 0, 0.4, False)   # 每例单独一个文件（多进程共用一个文件时曾读到半个文件）
    ev = c_json(["eval", txt, W, 0] + base["Nk"] + ["|"] + base["Mj"] + ["|"] + opts_eval(opts) + [f"rdd={off}"])
    ndev = ev["ndec_w"]
    rng = random.Random(seed)
    rho = rng.choice([0.02, 0.1, 0.4]); qon = 0.4; ro = gi % 2 == 1
    K = max(1, math.ceil(rho * ndev / qon - 1e-9))
    g = [rng.random() for _ in range(2 * K)]
    gp = os.path.join(gdir, f"{name}__{rule}__g{gi}.txt"); write_geno(gp, g, W, ndev, qon, ro)
    D = ref.decode_geno(g, W, ndev, qon, ro)
    rec = dict(inst=name, rule=rule, g=gi, rho=rho, K=K, ndev=ndev, resetopt=int(ro), n_active=sum(1 for u in D["dev"] if u >= 0))
    Nk = [x + rng.randint(0, 1) for x in base["Nk"]]; Mj = [x + rng.randint(0, 1) for x in base["Mj"]]
    bad = []
    for Tmax, cert in ((3000, False), (TC, True)):
        py = ref.run(I, Nk, Mj, Tmax, cert, rdd=D, **pol)
        cc = c_json(["eval", txt, Tmax, int(cert)] + Nk + ["|"] + Mj + ["|"] + opts_eval(opts) + [f"rdd={gp}"])
        for k in ("gaps", "cert", "period", "resets", "ndec_w"):
            if not cert and k in ("cert", "period"): continue
            if py[k] != cc[k]: bad.append((Tmax, k, py[k], cc[k]))
    rec["eval_ok"] = int(not bad)
    cs = c_json(["size", txt, HS, TC, 600000] + opts + [f"rdd={gp}"])
    R, pN, pM, st, per, fac = ref.size_rule(I, dict(pol, rdd=D), HS, TC, N0, M0, is_rlp, facloop=facloop)
    rec.update(size_ok=int((R, pN, pM, st) == (cs["R"], cs["Nk"], cs["Mj"], cs["status"])), R_py=R, R_c=cs["R"], status=st)
    if not rec["size_ok"]: bad.append(("size", (R, pN, pM, st), (cs["R"], cs["Nk"], cs["Mj"], cs["status"])))
    rec["bad"] = bad
    return rec


def opts_eval(opts):
    """eval 命令的第 4 段只认 rlp= / lpg= / curdef / rdd=，权重放在第 3 段：把 w= 拆出来放前面。"""
    w = [o for o in opts if o.startswith("w=")]
    rest = [o for o in opts if not o.startswith("w=") and o != "facloop"]
    ws = w[0][2:].split(",") if w else []
    return ws + ["|"] + rest


def off_check(name, rule):
    """全关基因型 = 默认规则：size 与 eval 输出逐字段相同（除用时）。"""
    I = load_train(name); _, txt = inst_paths(name, "train")
    opts = t8_rules.rule_opts(name, rule); W = 2 * S_max(I)
    gdir = os.path.join(OUT, "geno"); os.makedirs(gdir, exist_ok=True)
    off = os.path.join(gdir, f"{name}__{rule}__alloff.txt")
    write_geno(off, [0.3, 0.1, 0.7, 0.5], W, 500, 0.4, True)       # 两个槽都未激活（a < 0.6）
    a = c_json(["size", txt, HS, TC, 600000] + opts); b = c_json(["size", txt, HS, TC, 600000] + opts + [f"rdd={off}"])
    for d in (a, b): d.pop("sec"); d.pop("ndec_w")
    ok = a == b
    if a["R"] >= 0:
        for Tmax, cert in ((3000, 0), (TC, 1)):
            x = c_json(["eval", txt, Tmax, cert] + a["Nk"] + ["|"] + a["Mj"] + ["|"] + opts_eval(opts))
            y = c_json(["eval", txt, Tmax, cert] + a["Nk"] + ["|"] + a["Mj"] + ["|"] + opts_eval(opts) + [f"rdd={off}"])
            x.pop("ndec_w"); y.pop("ndec_w"); ok = ok and x == y
    return dict(inst=name, rule=rule, alloff_ok=int(ok), R=a["R"])


def main():
    names = sys.argv[1].split(","); ng = int(sys.argv[2]); P = int(sys.argv[3])
    os.makedirs(OUT, exist_ok=True)
    if os.path.exists(os.path.join(OUT, sys.argv[4] if len(sys.argv) > 4 else "check_rdd.jsonl")): sys.exit("输出文件已存在，不覆盖")
    out = os.path.join(OUT, sys.argv[4] if len(sys.argv) > 4 else "check_rdd.jsonl")
    cases = [(n, r, gi, sum(map(ord, n + r)) * 1000 + gi) for n in names for r in ("rlp", "es", "hand") for gi in range(ng)]   # 种子只依赖实例名、规则、序号（可复现）
    with ProcessPoolExecutor(P) as ex, open(out + ".tmp", "w") as f:
        for rec in ex.map(off_check, [n for n in names for r in ("rlp", "es", "hand")], [r for n in names for r in ("rlp", "es", "hand")]):
            f.write(json.dumps(rec, ensure_ascii=False) + "\n"); f.flush(); print(rec, flush=True)
        for rec in ex.map(one, cases):
            f.write(json.dumps(rec, ensure_ascii=False) + "\n"); f.flush()
            if rec.get("bad"): print("不一致", rec, flush=True)
    os.replace(out + ".tmp", out)
    recs = [json.loads(l) for l in open(out)]
    g = [r for r in recs if "g" in r and "skip" not in r]
    print(json.dumps(dict(alloff=[sum(r["alloff_ok"] for r in recs if "alloff_ok" in r), sum(1 for r in recs if "alloff_ok" in r)],
                          eval_mismatch=sum(1 - r["eval_ok"] for r in g), size_mismatch=sum(1 - r["size_ok"] for r in g), cases=len(g),
                          skipped=sum(1 for r in recs if "skip" in r)), ensure_ascii=False))


if __name__ == "__main__":
    main()
