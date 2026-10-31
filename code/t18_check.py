"""T18 的一致性核对（rost5 改动：resetopt 的"本班后提前大修"备选项只在 t < planW 时提供）。
(a) 新路径（偏离窗口 W ≠ 提前大修窗口 planW = 2·S_max）：rost5（新）与 Python 参考实现（code/rost5_ref.py，同步改动）逐字段比较，
    W ∈ {0.5, 1, 3, 4}·S_max，fresh_rlp_w 规则 + 严格窗口 + 随机基因型（ρ ∈ {0.1, 0.4}，q_on = 0.4，resetopt），
    台数 = 默认规则定出的台数 + 0..1 随机增量；
    Tmax = 3000 不判环：缺岗、提前大修次数、最后一次提前大修时刻、窗口内决策数；Tcert = 20000 判环：缺岗、判环结果、周期；
    另查最后一次提前大修 < 2·S_max（窗口外不提前大修），并记下旧程序在同一组上的最后一次提前大修（说明改动是否必要）。
(b) 旧路径不变：新旧两个 rost5 在 T12 路径（W = planW = 2·S_max）与 T11 路径（不设 planW；只有 m ≤ 6 有规则文件）上的 eval 输出逐字段相同。
输出 results/T18/check.jsonl，打印不一致数。
用法：python3 code/t18_check.py <旧 rost5> <新 rost5>"""
import os, sys, json, random, subprocess, tempfile, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import rost5_ref as ref  # noqa: E402
from t11_job import write_geno  # noqa: E402
import t8_rules  # noqa: E402

RULES = os.path.join(ROOT, "results", "T6", "rules")
INSTS = ["m04-K1-00", "m06-K12-01", "m08-K1-02", "m10-K12-03"]


def c_eval(r5, txt, T, cert, Nk, Mj, opts):
    """rost5 eval，返回 JSON。"""
    a = [r5, "eval", txt, str(T), str(cert)] + [str(x) for x in Nk] + ["|"] + [str(x) for x in Mj] + ["|", "|"] + opts
    return json.loads(subprocess.run(a, capture_output=True, text=True, check=True).stdout)


def base_fleet(n):
    """默认规则 fresh_rlp_w 的台数：m ≥ 8 取 T12 base（严格窗口），否则取 T6 的结果。"""
    bp = os.path.join(ROOT, "results", "T12", "jobs", f"{n}__base.json")
    b = json.load(open(bp)) if os.path.exists(bp) else json.load(open(os.path.join(ROOT, "results", "T6", "rlp", f"{n}__fresh_rlp_w.json")))
    return b["Nk"], b["Mj"]


def rand_geno(rng, tmpd, tag, W, ndev):
    """随机基因型：写文件并返回 (路径, Python 解码结果)。"""
    rho = rng.choice([0.1, 0.4]); K = max(1, math.ceil(rho * ndev / 0.4 - 1e-9))
    g = [rng.random() for _ in range(2 * K)]
    gp = os.path.join(tmpd, f"{tag}.geno"); write_geno(gp, W, ndev, 0.4, g, 1)
    return gp, ref.decode_geno(g, W, ndev, 0.4, True), rho


def main():
    old, new = sys.argv[1], sys.argv[2]
    rng = random.Random(20261031); out = []; bad = 0
    with tempfile.TemporaryDirectory(prefix="t18chk_") as tmpd:
        alloff = os.path.join(tmpd, "alloff.geno")
        for n in INSTS:
            I = load_train(n); S = S_max(I); P2 = 2 * S; _, txt = inst_paths(n, "train")
            rp, pp = f"{RULES}/{n}__fresh_rlp_w.rlp", f"{RULES}/{n}__fresh_rlp_w.planf"
            R = ref.parse_rlp(I, rp); PF = ref.parse_planf(I, pp)
            opts0 = [f"rlp={rp}", f"planf={pp}", f"planW={P2}"]
            Nk0, Mj0 = base_fleet(n)
            # (a) 新路径
            for mult in (0.5, 1, 3, 4):
                W = int(mult * S); write_geno(alloff, W, 0, 0.4, [], 1)
                ndev = c_eval(new, txt, W, 0, Nk0, Mj0, opts0 + [f"rdd={alloff}"])["ndec_w"]
                for c in range(2):
                    gp, D, rho = rand_geno(rng, tmpd, f"a{mult}_{c}", W, ndev)
                    Nk = [x + rng.randint(0, 1) for x in Nk0]; Mj = [x + rng.randint(0, 1) for x in Mj0]
                    ce = c_eval(new, txt, 3000, 0, Nk, Mj, opts0 + [f"rdd={gp}"])
                    pe = ref.run(I, Nk, Mj, 3000, False, rlp=R, planf=PF, rdd=D, plan_until=P2)
                    cc = c_eval(new, txt, 20000, 1, Nk, Mj, opts0 + [f"rdd={gp}"])
                    pc = ref.run(I, Nk, Mj, 20000, True, rlp=R, planf=PF, rdd=D, plan_until=P2)
                    a = dict(gaps=ce["gaps"], resets=ce["resets"], ndec=ce["ndec_w"], cgaps=cc["gaps"], cert=cc["cert"], period=cc["period"])
                    p = dict(gaps=pe["gaps"], resets=pe["resets"], ndec=pe["ndec_w"], cgaps=pc["gaps"], cert=pc["cert"], period=pc["period"])
                    lr = ce.get("last_reset", -1)
                    lr_old = c_eval(old, txt, 3000, 0, Nk, Mj, opts0 + [f"rdd={gp}"]).get("last_reset", -1)   # 旧程序（不限备选项）的最后一次提前大修，只记录
                    ok = a == p and lr < P2; bad += (not ok)
                    out.append(dict(part="new_path", inst=n, mult=mult, W=W, planW=P2, case=c, rho=rho, ndev=ndev,
                                    n_active=sum(u >= 0 for u in D["dev"]), cpp=a, py=p, last_reset=lr, last_reset_old=lr_old, same=int(ok)))
                    print("新路径", n, f"W={mult}S", c, "一致" if ok else "不一致", a, "last_reset", lr, "旧程序", lr_old, "" if ok else p, flush=True)
            # (b) 旧路径：新旧程序逐字段相同
            write_geno(alloff, P2, 0, 0.4, [], 1)
            ndev = c_eval(old, txt, P2, 0, Nk0, Mj0, opts0 + [f"rdd={alloff}"])["ndec_w"]
            t11 = t8_rules.rule_opts(n, "rlp")
            for path, opts in (("T12", opts0),) + ((("T11", [o for o in t11 if not o.startswith("hsfb=")]),) if I.meta["m"] <= 6 else ()):   # T3 的 fresh_noreset_u 只有 m ≤ 6
                for c in range(3):
                    gp, D, rho = rand_geno(rng, tmpd, f"b{path}_{c}", P2, ndev)
                    Nk = [x + rng.randint(0, 1) for x in Nk0]; Mj = [x + rng.randint(0, 1) for x in Mj0]
                    diffs = []
                    for T, cert in ((3000, 0), (60000, 1)):
                        eo = c_eval(old, txt, T, cert, Nk, Mj, opts + [f"rdd={gp}"]); en = c_eval(new, txt, T, cert, Nk, Mj, opts + [f"rdd={gp}"])
                        eo.pop("sec", None); en.pop("sec", None)
                        if eo != en: diffs.append((T, {k: (eo.get(k), en.get(k)) for k in set(eo) | set(en) if eo.get(k) != en.get(k)}))
                    ok = not diffs; bad += (not ok)
                    out.append(dict(part="old_path", path=path, inst=n, case=c, rho=rho, same=int(ok), diff=diffs))
                    print("旧路径", path, n, c, "相同" if ok else diffs, flush=True)
    os.makedirs(os.path.join(ROOT, "results", "T18"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "T18", "check.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(json.dumps(dict(cases=len(out), mismatch=bad)))


if __name__ == "__main__":
    main()
