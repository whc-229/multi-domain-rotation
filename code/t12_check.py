"""T12 的一致性核对：rost5 与 Python 参考实现（code/rost5_ref.py）在 T12 用到的组合上逐字段比较。
组合 = fresh_rlp_w 规则（rlp + planf）+ 严格窗口（planW = W）+ RDD 基因型（带"本班后提前大修"备选项 resetopt）。
T7 的核对不含 planf 与 planW，所以单独核对这一条路径。
每组：随机基因型（ρ ∈ {0.02, 0.1, 0.4}，q_on = 0.4）、台数 = 该规则在 T12 base 作业里定出的台数 + 0..1 随机增量（缺 base 时用 T6 的台数）；
  (1) Tmax = 3000 不判环：缺岗数、提前大修次数、窗口内决策数；
  (2) Tcert = 20000 判环：缺岗数、判环结果、周期。
输出 results/T12/check_planw.jsonl，打印不一致的组数。
用法：python3 code/t12_check.py <实例,...> <每个实例的组数>"""
import os, sys, json, random, subprocess, tempfile, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import rost5_ref as ref  # noqa: E402
from t11_job import write_geno  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")
RULES = os.path.join(ROOT, "results", "T6", "rules")


def c_eval(txt, T, cert, Nk, Mj, opts):
    """rost5 eval，返回 JSON。"""
    a = [R5, "eval", txt, str(T), str(cert)] + [str(x) for x in Nk] + ["|"] + [str(x) for x in Mj] + ["|", "|"] + opts
    return json.loads(subprocess.run(a, capture_output=True, text=True, check=True).stdout)


def main():
    names, n_per = sys.argv[1].split(","), int(sys.argv[2])
    rng = random.Random(20261030)
    out = []; bad = 0
    with tempfile.TemporaryDirectory(prefix="t12chk_") as tmpd:
        for n in names:
            I = load_train(n); S = S_max(I); W = 2 * S; _, txt = inst_paths(n, "train")
            rp, pp = f"{RULES}/{n}__fresh_rlp_w.rlp", f"{RULES}/{n}__fresh_rlp_w.planf"
            R = ref.parse_rlp(I, rp); PF = ref.parse_planf(I, pp)
            bp = os.path.join(ROOT, "results", "T12", "jobs", f"{n}__base.json")
            b = json.load(open(bp)) if os.path.exists(bp) else json.load(open(os.path.join(ROOT, "results", "T6", "rlp", f"{n}__fresh_rlp_w.json")))
            opts0 = [f"rlp={rp}", f"planf={pp}", f"planW={W}"]
            ndev = c_eval(txt, W, 0, b["Nk"], b["Mj"], opts0 + [f"rdd={_alloff(tmpd, W)}"])["ndec_w"]
            for c in range(n_per):
                rho = rng.choice([0.02, 0.1, 0.4]); K = max(1, math.ceil(rho * ndev / 0.4 - 1e-9))
                g = [rng.random() for _ in range(2 * K)]
                Nk = [x + rng.randint(0, 1) for x in b["Nk"]]; Mj = [x + rng.randint(0, 1) for x in b["Mj"]]
                gp = os.path.join(tmpd, f"g{c}.geno"); write_geno(gp, W, ndev, 0.4, g, 1)
                D = ref.decode_geno(g, W, ndev, 0.4, True)
                ce = c_eval(txt, 3000, 0, Nk, Mj, opts0 + [f"rdd={gp}"])
                pe = ref.run(I, Nk, Mj, 3000, False, rlp=R, planf=PF, rdd=D, plan_until=W)
                cc = c_eval(txt, 20000, 1, Nk, Mj, opts0 + [f"rdd={gp}"])
                pc = ref.run(I, Nk, Mj, 20000, True, rlp=R, planf=PF, rdd=D, plan_until=W)
                a = dict(gaps=ce["gaps"], resets=ce["resets"], ndec=ce["ndec_w"], cgaps=cc["gaps"], cert=cc["cert"], period=cc["period"])
                p = dict(gaps=pe["gaps"], resets=pe["resets"], ndec=pe["ndec_w"], cgaps=pc["gaps"], cert=pc["cert"], period=pc["period"])
                ok = a == p; bad += (not ok)
                out.append(dict(inst=n, case=c, rho=rho, K=K, n_active=sum(u >= 0 for u in D["dev"]), cpp=a, py=p, same=int(ok)))
                print(n, c, "一致" if ok else "不一致", a, p if not ok else "", flush=True)
    os.makedirs(os.path.join(ROOT, "results", "T12"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "T12", "check_planw.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(json.dumps(dict(cases=len(out), mismatch=bad)))


def _alloff(tmpd, W):
    """全关基因型文件（ndev = 0），用来数窗口内决策数。"""
    p = os.path.join(tmpd, "alloff.geno"); write_geno(p, W, 0, 0.4, [], 1); return p


if __name__ == "__main__":
    main()
