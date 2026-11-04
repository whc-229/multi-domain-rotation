"""T19 (b) 一致性核对（rost5 / rost5_ref 11-05 的改动：devfirst 偏离优先）。
(a) 新路径：C++（新）与 Python 参考实现逐字段比较。默认规则 = (a) 的接入构造（code/t19b_job.py 的 conn_spec：稳态规则 + 计划 + planxlate + 记账相位 + planW），
    加 devfirst 与随机 RDD 基因型（偏离窗口 W = ⌊0.5·S_max⌋，ρ ∈ {0.1, 0.4}，resetopt）；台数为构造台数（第 0 组）或随机 + 0..1。
    每组两次运行：3000 格不判环（缺岗、提前大修次数、窗口决策数）与 2 万格判环（缺岗、判环结果、周期）。5 个实例 × 4 组 = 20 组。
(b) 旧路径：不加 devfirst 时新旧程序 eval 输出逐字段相同（去掉 sec），3000 格不判环与 6 万格判环：
    t19_check 的 6 条旧路径（5 列计划 + PB、松弛规则、T12 路径、只算机器人、稳态初值、线性规则 w + curdef），
    另加"接入构造 + RDD 基因型、不加 devfirst"（计划与偏离同时出现、按原来的计划优先）；3 个实例 × 7 条路径 × 2 组 = 42 组。
输出 results/T19/check_b.jsonl（不动 T19 (a) 的 check.jsonl），打印不一致数。用法：python3 code/t19b_check.py <旧 rost5> <新 rost5>"""
import os, sys, json, random, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import rost5_ref as ref  # noqa: E402
import t13_methods as TM  # noqa: E402
from t11_job import write_geno  # noqa: E402
from t15b_check import rand_plan  # noqa: E402
from t18_check import rand_geno  # noqa: E402
from t19_check import run_c, ev_args  # noqa: E402
import t19b_job as B  # noqa: E402

NEW_INSTS = ["m04-K1-00", "m04-K12-03", "m06-K12-04", "m08-K1-02", "m10-K12-05"]
OLD_INSTS = ["m04-K1-01", "m06-K1-00", "m08-K12-03"]


def ndev_of(r5, txt, W, Nk, Mj, opts, tmpd):
    """窗口决策数（全关的基因型跑到 W）。"""
    gp = os.path.join(tmpd, "alloff.geno"); write_geno(gp, W, 0, 0.4, [], 1)
    return run_c(r5, ev_args(txt, W, 0, Nk, Mj, None, opts + [f"rdd={gp}"]))["ndec_w"]


def main():
    """(a) 新路径 C++ 与 Python；(b) 旧路径新旧程序。"""
    old, new = sys.argv[1], sys.argv[2]
    rng = random.Random(20261105); out = []; bad = 0
    with tempfile.TemporaryDirectory(prefix="t19bchk_") as tmpd:
        # (a) 新路径
        for n in NEW_INSTS:
            sp = B.conn_spec(n); I, txt, S = sp["I"], sp["txt"], sp["S"]; W = int(B.MULT * S)
            nd = ndev_of(new, txt, W, sp["Nk"], sp["Mj"], sp["opts"], tmpd)
            for c in range(4):
                gp, D, rho = rand_geno(rng, tmpd, f"{n}_{c}", W, nd)
                Nk = list(sp["Nk"]) if c == 0 else [x + rng.randint(0, 1) for x in sp["Nk"]]
                Mj = list(sp["Mj"]) if c == 0 else [x + rng.randint(0, 1) for x in sp["Mj"]]
                diff = []
                for Tm, cert in ((3000, 0), (20000, 1)):
                    cc = run_c(new, ev_args(txt, Tm, cert, Nk, Mj, None, sp["opts"] + [f"rdd={gp}"]))
                    py = ref.run(I, Nk, Mj, Tm, bool(cert), rdd=D, **sp["py"])
                    keys = ("gaps", "cert", "period", "resets") if cert else ("gaps", "resets", "ndec_w")
                    diff += [(Tm, k, py[k], cc[k]) for k in keys if py[k] != cc[k]]
                bad += bool(diff)
                out.append(dict(part="new_path", inst=n, case=c, rho=rho, ndev=nd, n_active=sum(u >= 0 for u in D["dev"]),
                                cpp=dict(gaps=cc["gaps"], cert=cc["cert"], period=cc["period"], resets=cc["resets"]), same=int(not diff), diff=diff))
                print("新路径", n, c, f"ρ={rho}", "一致" if not diff else diff, flush=True)
        # (b) 旧路径
        pf = os.path.join(tmpd, "old.planx"); gp = os.path.join(tmpd, "g.geno")
        for n in OLD_INSTS:
            I = load_train(n); S = S_max(I); W2 = 2 * S; _, txt = inst_paths(n, "train")
            pb = run_c(new, ["size", txt, 365, 60000, 60000] + TM.opts("pb", n, S))
            st = json.load(open(os.path.join(ROOT, "results", "T3", "jobs", f"{n}__steady.json")))
            cs = B.conn_spec(n)
            for path in ("rhmip", "rlp", "t12", "notools", "steady", "wcurdef", "conn_rdd"):
                for c in range(2):
                    base = st if path == "steady" else (dict(Nk=cs["Nk"], Mj=cs["Mj"]) if path == "conn_rdd" else pb)
                    Nk = [x + rng.randint(0, 1) for x in base["Nk"]]; Mj = [x + rng.randint(0, 1) for x in base["Mj"]]
                    w = None
                    if path in ("rhmip", "notools"):
                        rand_plan(rng, I, 4 * S, W2, pf); w = TM.W_PB; o = [f"planx={pf}", f"planW={W2}"] + (["notools"] if path == "notools" else [])
                        if path == "notools": Mj = []
                    elif path == "rlp":     # fresh_noreset_u 只有 m ≤ 6 有（T3b），m = 8 用 fresh_noreset
                        u = f"{ROOT}/results/T3/rlp/{n}__fresh_noreset_u.rlp"; o = [f"rlp={u if os.path.exists(u) else u.replace('_u.rlp', '.rlp')}"]
                    elif path == "steady": o = [f"rlp={ROOT}/results/T3/rlp/{n}__steady.rlp"]
                    elif path == "wcurdef": w = [round(x + rng.gauss(0, 0.3), 4) for x in TM.W_PB]; o = ["curdef"]
                    elif path == "conn_rdd":
                        o0 = [x for x in cs["opts"] if x != "devfirst"]; W = int(B.MULT * S)
                        nd = ndev_of(new, txt, W, Nk, Mj, o0, tmpd)
                        write_geno(gp, W, nd, 0.4, [rng.random() for _ in range(2 * max(1, nd // 4))], 1)
                        o = o0 + [f"rdd={gp}"]
                    else:
                        R6 = os.path.join(ROOT, "results", "T6", "rules")
                        o12 = [f"rlp={R6}/{n}__fresh_rlp_w.rlp", f"planf={R6}/{n}__fresh_rlp_w.planf", f"planW={W2}"]
                        nd = ndev_of(new, txt, W2, Nk, Mj, o12, tmpd)
                        write_geno(gp, W2, nd, 0.4, [rng.random() for _ in range(2 * max(1, nd // 4))], 1)
                        o = o12 + ["resetopt", f"rdd={gp}"]
                    diffs = []
                    for Tm, cert in ((3000, 0), (60000, 1)):
                        eo = run_c(old, ev_args(txt, Tm, cert, Nk, Mj, w, o)); en = run_c(new, ev_args(txt, Tm, cert, Nk, Mj, w, o))
                        if eo != en: diffs.append((Tm, {k: (eo.get(k), en.get(k)) for k in set(eo) | set(en) if eo.get(k) != en.get(k)}))
                    bad += bool(diffs)
                    out.append(dict(part="old_path", inst=n, path=path, case=c, same=int(not diffs), diff=diffs))
                    print("旧路径", n, path, c, "相同" if not diffs else diffs, flush=True)
    with open(os.path.join(ROOT, "results", "T19", "check_b.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(json.dumps(dict(cases=len(out), new_path=sum(1 for r in out if r["part"] == "new_path"), old_path=sum(1 for r in out if r["part"] == "old_path"), mismatch=bad)))


if __name__ == "__main__":
    main()
