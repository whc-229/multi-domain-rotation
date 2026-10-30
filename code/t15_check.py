"""T15 / T13 补充核对（rost5 新增路径）：
(1) notools：rost5 带 notools 与"每类工具 1000 件"的同一组运行，缺岗、提前大修次数、值岗次数相同（工具足够多时不影响派工），
    规则 rlp / es / hand，实例 m02-K1-00、m03-K12-00，每个 (规则, 实例) 5 组随机台数，Tmax = 3000 与 60000 格，都不判环
    （工具多时状态里多了工具，判出循环的时刻可以不同）；
(2) 判环路径：T13 的逐决策编码（rk、d5）、GP 树、按计划派在 60000 格判环下与 Python 参考实现相同（T13 核对只跑到 20000 格，都没判出循环），
    实例 m02-K1-00，每种 3 组。
输出 results/T15/check.jsonl，打印不一致数。用法：python3 code/t15_check.py [rost5 路径]"""
import os, sys, json, random, subprocess, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import rost5_ref as ref  # noqa: E402
import t8_rules  # noqa: E402
import t13_methods as TM  # noqa: E402
from t13_check import rand_tree  # noqa: E402

R5 = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "code", "cpp", "rost5")


def c_json(args):
    """运行 rost5，解析一行 JSON。"""
    return json.loads(subprocess.run([R5] + [str(a) for a in args], capture_output=True, text=True, check=True).stdout)


def opts_eval(opts):
    """eval 第 3、4 段。"""
    w = [o for o in opts if o.startswith("w=")]
    return (w[0][2:].split(",") if w else []) + ["|"] + [o for o in opts if not o.startswith("w=") and o != "facloop"]


def main():
    rng = random.Random(15); out = []; bad = 0; tmpd = tempfile.mkdtemp(prefix="t15chk_")
    keys = ("gaps", "resets", "duties")
    for n in ["m02-K1-00", "m03-K12-00"]:
        I = load_train(n); _, txt = inst_paths(n, "train")
        for rule in ("rlp", "es", "hand"):
            opts = t8_rules.rule_opts(n, rule)
            b = c_json(["size", txt, 365, 60000, 60000] + opts)
            for c in range(5):
                Nk = [x + rng.randint(-1, 1) for x in b["Nk"]]; Nk = [max(0, x) for x in Nk]
                diff = []
                for Tm, cert in ((3000, 0), (60000, 0)):     # 不判环：工具多时状态里多了工具，判出循环的时刻可以不同，只比派工动态
                    a = c_json(["eval", txt, Tm, cert] + Nk + ["|"] + [1000] * len(b["Mj"]) + ["|"] + opts_eval(opts))
                    z = c_json(["eval", txt, Tm, cert] + Nk + ["|", "|"] + opts_eval(opts) + ["notools"])
                    diff += [(Tm, k, a[k], z[k]) for k in keys if a[k] != z[k]]
                bad += bool(diff); out.append(dict(part="notools", inst=n, rule=rule, case=c, same=int(not diff), diff=diff))
                print("notools", n, rule, c, "相同" if not diff else diff, flush=True)
    n = "m02-K1-00"; I = load_train(n); S = S_max(I); W = 2 * S; _, txt = inst_paths(n, "train")
    for meth in ("rk", "d5", "gphh", "rhmip"):
        b = c_json(["size", txt, 365, 60000, 60000] + TM.opts("pb", n, S))
        ndev = c_json(["eval", txt, W, 0] + b["Nk"] + ["|"] + b["Mj"] + ["|", "|", "pd=rk", f"pdW={W}"])["ndec_w"]
        for c in range(3):
            Nk = [x + rng.randint(0, 1) for x in b["Nk"]]; Mj = [x + rng.randint(0, 1) for x in b["Mj"]]
            if meth in ("rk", "d5"):
                g = [rng.random() * 0.3 for _ in range(ndev)]
                pol = dict(w=[float(x) for x in TM.W_TOOL] + [0.0] * 6, curdef=True, pd=dict(mode=1 if meth == "rk" else 2, W=W, g=g, repair=(meth == "d5")))
                ev = TM.W_TOOL + ["|", f"pd={meth}", f"pdW={W}", "pdg=" + ",".join(repr(x) for x in g)]
            elif meth == "gphh":
                tr = rand_tree(rng); tf = os.path.join(tmpd, "t.tree"); open(tf, "w").write(tr)
                pol = dict(w=[float(x) for x in TM.W_TOOL] + [0.0] * 6, tree=ref.parse_tree(tr)); ev = TM.W_TOOL + ["|", f"tree=@{tf}"]
            else:
                pf = os.path.join(tmpd, "p.planx")
                with open(pf, "w") as f:
                    f.write(f"end {4 * S}\n")
                    for D in range(4 * S):
                        for s2 in range(len(I.sites)):
                            if rng.random() < 0.3:
                                k = rng.choice([k for k in range(len(I.T)) if s2 in I.elig[k]]); f.write(f"{s2} {D} {k} {rng.randrange(I.L[k])} {int(D < W and rng.random() < 0.2)}\n")
                pol = dict(w=[float(x) for x in TM.W_PB] + [0.0] * 6, planx=ref.parse_planx(pf), plan_until=W); ev = TM.W_PB + ["|", f"planx={pf}", f"planW={W}"]
            cc = c_json(["eval", txt, 60000, 1] + Nk + ["|"] + Mj + ["|"] + ev)
            py = ref.run(I, Nk, Mj, 60000, True, **pol)
            diff = [(k, py[k], cc[k]) for k in ("gaps", "cert", "period", "resets") if py[k] != cc[k]]
            bad += bool(diff); out.append(dict(part="cert", method=meth, case=c, cpp={k: cc[k] for k in ("gaps", "cert", "period")}, same=int(not diff), diff=diff))
            print("判环", meth, c, "一致" if not diff else diff, {k: cc[k] for k in ("gaps", "cert", "period")}, flush=True)
    os.makedirs(os.path.join(ROOT, "results", "T15"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "T15", "check.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(json.dumps(dict(cases=len(out), mismatch=bad)))


if __name__ == "__main__":
    main()
