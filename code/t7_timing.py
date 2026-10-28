"""T7：单次 RDD 解码（= 一次完整的定台数 + 证书，即进化里一次适应度评估，不含 60 万格复核）的用时。
每个 m ∈ {2, 3, 4, 6} 取 K1、K12 各 2 个训练实例，3 种默认规则，每种 3 个随机基因型（ρ = 0.1、q_on = 0.4）和全关基因型各一次。
用 rost5 size 计时（size 的 sec 含 near_miss 与复核，这里另跑一次 Tverify = 0 的 size，只取定台数部分）。
输出 results/T7/timing.jsonl 和按 (m, 规则) 的表。用法：python3 code/t7_timing.py"""
import os, sys, json, random, subprocess, time, math
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import t8_rules  # noqa: E402
from t7_check import write_geno, opts_eval, c_json  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")


def main():
    out = os.path.join(ROOT, "results", "T7", "timing.jsonl"); gdir = os.path.join(ROOT, "results", "T7", "geno"); os.makedirs(gdir, exist_ok=True)
    rows = []
    for m in (2, 3, 4, 6):
        for f in ("K1", "K12"):
            for i in (0, 1):
                n = f"m{m:02d}-{f}-{i:02d}"; I = load_train(n); _, txt = inst_paths(n, "train"); W = 2 * S_max(I)
                for rule in ("rlp", "es", "hand"):
                    opts = t8_rules.rule_opts(n, rule)
                    base = c_json(["size", txt, 365, 60000, 0] + opts)
                    if base["R"] < 0: rows.append(dict(inst=n, m=m, rule=rule, skip=1)); continue
                    off = os.path.join(gdir, f"t_{n}_{rule}_off.txt"); write_geno(off, [0.0, 0.0], W, 0, 0.4, False)
                    ndev = c_json(["eval", txt, W, 0] + base["Nk"] + ["|"] + base["Mj"] + ["|"] + opts_eval(opts) + [f"rdd={off}"])["ndec_w"]
                    K = math.ceil(0.1 * ndev / 0.4 - 1e-9); rng = random.Random(m * 100 + i)
                    for gi in range(4):
                        g = [0.0, 0.0] if gi == 0 else [rng.random() for _ in range(2 * K)]
                        gp = os.path.join(gdir, f"t_{n}_{rule}_g{gi}.txt"); write_geno(gp, g, W, ndev, 0.4, False)
                        t0 = time.time(); r = c_json(["size", txt, 365, 60000, 0] + opts + [f"rdd={gp}"]); sec = time.time() - t0
                        rows.append(dict(inst=n, m=m, rule=rule, g=gi, ndev=ndev, K=K, R=r["R"], status=r["status"], sec=round(sec, 3)))
    with open(out, "w") as fo:
        for r in rows: fo.write(json.dumps(r, ensure_ascii=False) + "\n")
    by = defaultdict(list)
    for r in rows:
        if "sec" in r: by[(r["m"], r["rule"])].append(r)
    print("| m | 规则 | 次数 | N_dev（平均） | 单次解码平均 s | 最长 s | 1216 次估计（分钟） |\n|---|---|---|---|---|---|---|")
    for (m, rule), rs in sorted(by.items()):
        a = sum(r["sec"] for r in rs) / len(rs)
        print(f"| {m} | {rule} | {len(rs)} | {sum(r['ndev'] for r in rs) / len(rs):.0f} | {a:.2f} | {max(r['sec'] for r in rs):.2f} | {1216 * a / 60:.0f} |")


if __name__ == "__main__":
    main()
