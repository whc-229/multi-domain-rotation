"""T18 补扫更短窗口（11-02 审核：W ∈ {0.125, 0.25}·S_max）开跑前的核对。rost5 没有改动，这里只确认短窗口下新路径照常工作：
rost5 与 Python 参考实现逐字段比较，做法与 t18_check.py 的 (a) 相同（fresh_rlp_w + 严格窗口 planW = 2·S_max + 随机基因型 + resetopt，
台数 = 默认规则台数 + 0..1 随机增量；Tmax = 3000 不判环、Tcert = 20000 判环），另查最后一次提前大修 < 2·S_max。
实例与 t18_check 相同（m = 4、6、8、10 各 1 个），W ∈ {0.125, 0.25}·S_max，每个 (实例, W) 4 组。
输出 results/T18/check_short.jsonl（不动原来的 check.jsonl），打印不一致数。
用法：python3 code/t18s_check.py"""
import os, sys, json, random, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import rost5_ref as ref  # noqa: E402
from t11_job import write_geno  # noqa: E402
from t18_check import c_eval, base_fleet, rand_geno, RULES, INSTS  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")


def main():
    """逐个 (实例, W, 组) 比较 C++ 与 Python，写 check_short.jsonl。"""
    rng = random.Random(20261102); out = []; bad = 0
    with tempfile.TemporaryDirectory(prefix="t18schk_") as tmpd:
        alloff = os.path.join(tmpd, "alloff.geno")
        for n in INSTS:
            I = load_train(n); S = S_max(I); P2 = 2 * S; _, txt = inst_paths(n, "train")
            rp, pp = f"{RULES}/{n}__fresh_rlp_w.rlp", f"{RULES}/{n}__fresh_rlp_w.planf"
            R = ref.parse_rlp(I, rp); PF = ref.parse_planf(I, pp)
            opts0 = [f"rlp={rp}", f"planf={pp}", f"planW={P2}"]
            Nk0, Mj0 = base_fleet(n)
            for mult in (0.125, 0.25):
                W = int(mult * S); write_geno(alloff, W, 0, 0.4, [], 1)
                ndev = c_eval(R5, txt, W, 0, Nk0, Mj0, opts0 + [f"rdd={alloff}"])["ndec_w"]
                for c in range(4):
                    gp, D, rho = rand_geno(rng, tmpd, f"s{mult}_{c}", W, ndev)
                    Nk = [x + rng.randint(0, 1) for x in Nk0]; Mj = [x + rng.randint(0, 1) for x in Mj0]
                    ce = c_eval(R5, txt, 3000, 0, Nk, Mj, opts0 + [f"rdd={gp}"])
                    pe = ref.run(I, Nk, Mj, 3000, False, rlp=R, planf=PF, rdd=D, plan_until=P2)
                    cc = c_eval(R5, txt, 20000, 1, Nk, Mj, opts0 + [f"rdd={gp}"])
                    pc = ref.run(I, Nk, Mj, 20000, True, rlp=R, planf=PF, rdd=D, plan_until=P2)
                    a = dict(gaps=ce["gaps"], resets=ce["resets"], ndec=ce["ndec_w"], cgaps=cc["gaps"], cert=cc["cert"], period=cc["period"])
                    p = dict(gaps=pe["gaps"], resets=pe["resets"], ndec=pe["ndec_w"], cgaps=pc["gaps"], cert=pc["cert"], period=pc["period"])
                    lr = ce.get("last_reset", -1)
                    ok = a == p and lr < P2; bad += (not ok)
                    out.append(dict(part="short_path", inst=n, mult=mult, W=W, planW=P2, case=c, rho=rho, ndev=ndev,
                                    n_active=sum(u >= 0 for u in D["dev"]), cpp=a, py=p, last_reset=lr, same=int(ok)))
                    print(n, f"W={mult}S（{W} 格）", c, "一致" if ok else "不一致", a, "" if ok else p, flush=True)
    with open(os.path.join(ROOT, "results", "T18", "check_short.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(json.dumps(dict(cases=len(out), mismatch=bad)))


if __name__ == "__main__":
    main()
