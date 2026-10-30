"""T16 的核对：
(1) 位置分布：rost5 rddfd（0 代）输出的决策权重 w 与累积分布 cdf，与 Python 参考实现按同一规则、同一台数算出的逐个相同；
(2) 解码：随机基因型带 cdf，rost5 eval 与 Python run（decode_geno 带 cdf）逐字段相同（3000 格不判环、20000 格判环），每个实例 10 组；
(3) β = 0：rddfd 带对偶文件、β = 0 与不带对偶文件，输出除用时外逐字段相同（100 代）。
实例 m04-K1-00、m06-K12-04，β = 0.5。输出 results/T16/check.jsonl。用法：python3 code/t16_check.py [rost5 路径]"""
import os, sys, json, random, subprocess, tempfile, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import rost5_ref as ref  # noqa: E402
from t11_job import write_geno  # noqa: E402
from t12_job import rule_opts  # noqa: E402
from t16_job import dual_file  # noqa: E402

R5 = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "code", "cpp", "rost5")
BETA = 0.5


def c_json(args):
    """运行 rost5，解析一行 JSON。"""
    return json.loads(subprocess.run([R5] + [str(a) for a in args], capture_output=True, text=True, check=True).stdout)


def main():
    rng = random.Random(16); out = []; bad = 0; tmpd = tempfile.mkdtemp(prefix="t16chk_")
    for n in ["m04-K1-00", "m06-K12-04"]:
        I = load_train(n); S = S_max(I); W = 2 * S; _, txt = inst_paths(n, "train")
        opts = rule_opts(n, S); df = dual_file(n, W)
        base = ["rddfd", txt, 365, 60000, 600000, 6 * S]
        a = c_json(base + [0, 3600, 4, 12, 0, 0.1, 0.4] + opts + [f"dualw={df}", f"beta={BETA}", f"W={W}", "resetopt"])
        # (1) Python：默认规则在 F0 下的窗口决策
        R = ref.parse_rlp(I, opts[0][4:]); PF = ref.parse_planf(I, opts[1][6:])
        tr = {"light": False}
        ref.run(I, a["base_Nk"], a["base_Mj"], W, False, rlp=R, planf=PF, plan_until=W, rdd=ref.decode_geno([], W, 0, 0.4, True), trace=tr)
        pi = {}
        for line in open(df): s, t, v = line.split(); pi[(int(s), int(t))] = float(v)
        w = [sum(max(pi.get((s, u), 0.0), 0.0) for u in range(t, t + d)) for (s, t, d) in tr["decs"][:a["ndev"]]]
        sw = sum(w); cdf = []; acc = 0.0
        for c in range(a["ndev"]): acc += (1.0 - BETA) / a["ndev"] + (BETA * w[c] / sw if sw > 0 else BETA / a["ndev"]); cdf.append(acc)
        same_w = int(w == a["dualw"] and cdf == a["cdf"])
        out.append(dict(inst=n, part="weights", same=same_w, ndev=a["ndev"], n_pos=sum(x > 0 for x in w), max_abs_diff=max(abs(x - y) for x, y in zip(cdf, a["cdf"]))))
        bad += (not same_w); print(n, "权重与 cdf", "一致" if same_w else "不一致", flush=True)
        # (2) 解码
        cf = os.path.join(tmpd, f"{n}.cdf"); open(cf, "w").write("\n".join(repr(x) for x in a["cdf"]) + "\n")
        for c in range(10):
            K = a["K"]; g = [rng.random() for _ in range(2 * K)]
            gp = os.path.join(tmpd, "g.geno"); write_geno(gp, W, a["ndev"], 0.4, g, 1)
            D = ref.decode_geno(g, W, a["ndev"], 0.4, True, cdf=a["cdf"])
            Nk = [x + rng.randint(0, 1) for x in a["base_Nk"]]; Mj = [x + rng.randint(0, 1) for x in a["base_Mj"]]
            diff = []
            for Tm, cert in ((3000, False), (20000, True)):
                cc = c_json(["eval", txt, Tm, int(cert)] + Nk + ["|"] + Mj + ["|", "|"] + opts[:3] + [f"cdf={cf}", f"rdd={gp}"])
                py = ref.run(I, Nk, Mj, Tm, cert, rlp=R, planf=PF, plan_until=W, rdd=D)
                for k in ("gaps", "cert", "period", "resets", "ndec_w"):
                    if not cert and k in ("cert", "period"): continue
                    if py[k] != cc[k]: diff.append((Tm, k, py[k], cc[k]))
            bad += bool(diff); out.append(dict(inst=n, part="decode", case=c, n_active=sum(u >= 0 for u in D["dev"]), same=int(not diff), diff=diff))
            print(n, "解码", c, "一致" if not diff else diff, flush=True)
        # (3) β = 0 与原 RDD 逐位相同
        x = c_json(base + [100, 3600, 4, 12, 1, 0.1, 0.4] + opts + [f"dualw={df}", "beta=0", f"W={W}", "resetopt"])
        y = c_json(base + [100, 3600, 4, 12, 1, 0.1, 0.4] + opts + [f"W={W}", "resetopt"])
        skip = {"sec", "sec_eval", "sec_judge", "sec_dec", "base_sec", "sec_f0", "traj"}
        d = [k for k in x if k not in skip and x[k] != y.get(k)]
        same0 = int(not d and [t[:3] for t in x["traj"]] == [t[:3] for t in y["traj"]])
        bad += (not same0); out.append(dict(inst=n, part="beta0", same=same0, diff=d)); print(n, "β = 0", "相同" if same0 else d, flush=True)
    os.makedirs(os.path.join(ROOT, "results", "T16"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "T16", "check.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(json.dumps(dict(cases=len(out), mismatch=bad)))


if __name__ == "__main__":
    main()
