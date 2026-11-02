"""rost5 / rost5_ref 改动（11-02，T15 小检查 2）的一致性核对。改动两处：
- 加选项 planxlate 时，按计划派在同一 (机型, 相位) 的就绪候选里取就绪最晚的一台（相同时取原顺序靠前的）；不加时与原来相同；
- rost5 的 R_LP 派工分支接上按计划派（原来只有打分规则的分支有，Python 参考实现两个分支都有；R_LP + planx 以前没有任何结果用过）。
(a) 新路径：随机计划（planx，到 2·S_max，每 (格, 站) 30% 概率一条，含少量本班后提前大修），偶数组加 planxlate、奇数组不加，
    规则 PB（T13 的手工权重）与松弛规则 fresh_noreset_u（rlp），rost5（新）与 Python 参考实现逐字段比较：
    Tmax = 3000 不判环（缺岗、提前大修次数）；Tcert = 20000 判环（缺岗、判环结果、周期、提前大修次数）。实例 m02-K1-00、m03-K12-00、m04-K1-01，
    每个 (实例, 规则) 4 组（台数 = 规则定出的台数 + 0..1）。
(b) 旧路径不变：新旧两个 rost5 在不加 planxlate 的路径上 eval 输出逐字段相同（去掉 sec）：
    按计划派（RH-MIP 的路径，PB 兜底 + planW）、松弛规则 rlp、T12 路径（fresh_rlp_w + planf + planW + 随机 RDD 基因型）、只算机器人（notools + 按计划派），
    Tmax = 3000 不判环与 60000 判环。每个 (实例, 路径) 2 组。
输出 results/T15/check_b.jsonl，打印不一致数。用法：python3 code/t15b_check.py <旧 rost5> <新 rost5>"""
import os, sys, json, random, subprocess, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import rost5_ref as ref  # noqa: E402
import t13_methods as TM  # noqa: E402
import t8_rules  # noqa: E402
from t11_job import write_geno  # noqa: E402

INSTS = ["m02-K1-00", "m03-K12-00", "m04-K1-01"]
RULES6 = os.path.join(ROOT, "results", "T6", "rules")


def c_json(r5, args):
    """运行 rost5，返回 JSON（去掉 sec）。"""
    d = json.loads(subprocess.run([r5] + [str(a) for a in args], capture_output=True, text=True, check=True).stdout)
    d.pop("sec", None); return d


def rand_plan(rng, I, end, W, path):
    """随机计划文件：每 (格, 站) 30% 概率一条 (s, D, k, i, r)，r = 1（本班后提前大修）只在 D < W 时、20% 概率。"""
    with open(path, "w") as f:
        f.write(f"end {end}\n")
        for D in range(end):
            for s in range(len(I.sites)):
                if rng.random() < 0.3:
                    k = rng.choice([k for k in range(len(I.T)) if s in I.elig[k]])
                    f.write(f"{s} {D} {k} {rng.randrange(I.L[k])} {int(D < W and rng.random() < 0.2)}\n")


def main():
    """(a) 新路径 C++ 与 Python；(b) 旧路径新旧程序。"""
    old, new = sys.argv[1], sys.argv[2]
    rng = random.Random(20261102); out = []; bad = 0
    with tempfile.TemporaryDirectory(prefix="t15bchk_") as tmpd:
        pf = os.path.join(tmpd, "p.planx")
        for n in INSTS:
            I = load_train(n); S = S_max(I); W = 2 * S; _, txt = inst_paths(n, "train")
            rlp_opts = t8_rules.rule_opts(n, "rlp"); R = ref.parse_rlp(I, rlp_opts[0].split("=", 1)[1])
            pb = c_json(new, ["size", txt, 365, 60000, 60000] + TM.opts("pb", n, S))
            # (a) 新路径
            for rule in ("pb", "rlp"):
                base = pb if rule == "pb" else c_json(new, ["rlp", txt, rlp_opts[0].split("=", 1)[1], 365, 60000, 60000])
                for c in range(4):
                    late = c % 2 == 0
                    rand_plan(rng, I, W, W, pf); X = ref.parse_planx(pf); X["late"] = late; lt = ["planxlate"] if late else []
                    Nk = [x + rng.randint(0, 1) for x in base["Nk"]]; Mj = [x + rng.randint(0, 1) for x in base["Mj"]]
                    if rule == "pb":
                        pol = dict(w=[float(x) for x in TM.W_PB] + [0.0] * 6, planx=X, plan_until=W); ev = TM.W_PB + ["|", f"planx={pf}", f"planW={W}"] + lt
                    else:
                        pol = dict(rlp=R, planx=X, plan_until=W); ev = ["|"] + rlp_opts + [f"planx={pf}", f"planW={W}"] + lt
                    diff = []
                    for Tm, cert in ((3000, False), (20000, True)):
                        cc = c_json(new, ["eval", txt, Tm, int(cert)] + Nk + ["|"] + Mj + ["|"] + ev)
                        py = ref.run(I, Nk, Mj, Tm, cert, **pol)
                        diff += [(Tm, k, py[k], cc[k]) for k in (("gaps", "cert", "period", "resets") if cert else ("gaps", "resets")) if py[k] != cc[k]]
                    bad += bool(diff)
                    out.append(dict(part="new_path", inst=n, rule=rule, case=c, late=int(late), cpp=dict(gaps=cc["gaps"], cert=cc["cert"], period=cc["period"]), same=int(not diff), diff=diff))
                    print("新路径", n, rule, c, "planxlate" if late else "", "一致" if not diff else diff, flush=True)
            # (b) 旧路径
            gp = os.path.join(tmpd, "g.geno")
            for path in ("rhmip", "rlp", "t12", "notools"):
                for c in range(2):
                    Nk = [x + rng.randint(0, 1) for x in pb["Nk"]]; Mj = [x + rng.randint(0, 1) for x in pb["Mj"]]
                    if path in ("rhmip", "notools"):
                        rand_plan(rng, I, 4 * S, W, pf); ev = TM.W_PB + ["|", f"planx={pf}", f"planW={W}"] + (["notools"] if path == "notools" else [])
                        if path == "notools": Mj = []
                    elif path == "rlp": ev = ["|"] + rlp_opts
                    else:
                        o12 = [f"rlp={RULES6}/{n}__fresh_rlp_w.rlp", f"planf={RULES6}/{n}__fresh_rlp_w.planf", f"planW={W}"]
                        write_geno(gp, W, 0, 0.4, [], 1)     # 全关的基因型，先数窗口决策数（与 t18_check 相同）
                        ndev = c_json(new, ["eval", txt, W, 0] + Nk + ["|"] + Mj + ["|", "|"] + o12 + [f"rdd={gp}"])["ndec_w"]
                        write_geno(gp, W, ndev, 0.4, [rng.random() for _ in range(2 * max(1, ndev // 4))], 1)
                        ev = ["|"] + o12 + ["resetopt", f"rdd={gp}"]
                    diffs = []
                    for Tm, cert in ((3000, 0), (60000, 1)):
                        eo = c_json(old, ["eval", txt, Tm, cert] + Nk + ["|"] + Mj + ["|"] + ev); en = c_json(new, ["eval", txt, Tm, cert] + Nk + ["|"] + Mj + ["|"] + ev)
                        if eo != en: diffs.append((Tm, {k: (eo.get(k), en.get(k)) for k in set(eo) | set(en) if eo.get(k) != en.get(k)}))
                    bad += bool(diffs)
                    out.append(dict(part="old_path", inst=n, path=path, case=c, same=int(not diffs), diff=diffs))
                    print("旧路径", n, path, c, "相同" if not diffs else diffs, flush=True)
    with open(os.path.join(ROOT, "results", "T15", "check_b.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(json.dumps(dict(cases=len(out), mismatch=bad)))


if __name__ == "__main__":
    main()
