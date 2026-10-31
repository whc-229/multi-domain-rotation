"""T12b 的一致性核对（rost5 改动：es 命令加 planW= 选项，计划提前大修只在 t < planW 时执行）。
(1) 不加 planW 时新旧程序的 es 输出逐字段相同（除用时）：实例 m02-K1-03、m03-K12-00，2 代，种子 0、1，参数与 T6 的 es_lpgplan 相同；
(2) 加 planW = 2·S_max 时，es 给出的最终解（权重、台数）用 rost5 eval（同一组选项 + planW）重新判定：
    复核缺岗、判环结果、周期与 es 的输出相同，最后一次提前大修 < 2·S_max；
(3) ES 的评价路径（14 维权重 + lpg + T3 的 fresh_rlp.planf + curdef + planW = 2·S_max）上 rost5 eval 与 Python 参考实现逐字段相同：
    权重取 T6 各实例 ES 结果的权重，台数 = 其台数 + 0..1 随机增量；Tmax = 3000 不判环（缺岗、提前大修次数、最后一次提前大修）与 20000 判环（缺岗、判环、周期）。
输出 results/T12b/check.jsonl，打印不一致数。用法：python3 code/t12b_check.py <旧 rost5> <新 rost5>"""
import os, sys, json, random, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import rost5_ref as ref  # noqa: E402
from t13_check import lpz_of  # noqa: E402

T3R = os.path.join(ROOT, "results", "T3", "rlp")


def c_json(r5, args):
    """运行 rost5，解析一行 JSON。"""
    return json.loads(subprocess.run([r5] + [str(a) for a in args], capture_output=True, text=True, check=True).stdout)


def es_opts(n):
    """T6 es_lpgplan 的规则选项。"""
    return [f"lpg={T3R}/{n}.lpg.txt", f"planf={T3R}/{n}__fresh_rlp.planf"]


def main():
    old, new = sys.argv[1], sys.argv[2]
    rng = random.Random(1231); out = []; bad = 0
    for n in ("m02-K1-03", "m03-K12-00"):
        I = load_train(n); S = S_max(I); W = 2 * S; _, txt = inst_paths(n, "train")
        for seed in (0, 1):
            base = ["es", txt, "certs", 365, 2, 4, 12, seed, 60000, 600000] + es_opts(n)
            a = c_json(old, base); b = c_json(new, base)
            a.pop("sec", None); b.pop("sec", None)
            ok = a == b; bad += (not ok)
            out.append(dict(part="es_no_planW", inst=n, seed=seed, same=int(ok), R=b["R"]))
            print("es 不加 planW", n, seed, "相同" if ok else "不同", b["R"], flush=True)
            e = c_json(new, base + [f"planW={W}"])
            if e["R"] >= 0:
                v = c_json(new, ["eval", txt, 600000, 1] + e["Nk"] + ["|"] + e["Mj"] + ["|"] + e["w"] + ["|"] + es_opts(n) + ["curdef", f"planW={W}"])
                v3 = c_json(new, ["eval", txt, 3000, 0] + e["Nk"] + ["|"] + e["Mj"] + ["|"] + e["w"] + ["|"] + es_opts(n) + ["curdef", f"planW={W}"])
                same = (v["gaps"] == e["verify_gaps"] and v["cert"] == e["verify_cert"] and (e["verify_cert"] != 1 or v["period"] == e["verify_period"]))
                ok = same and v3["last_reset"] < W
            else: ok = True; v3 = {"last_reset": None}
            bad += (not ok)
            out.append(dict(part="es_planW_recheck", inst=n, seed=seed, R=e["R"], last_reset=v3["last_reset"], W=W, same=int(ok)))
            print("es 加 planW 后复判", n, seed, e["R"], "last_reset", v3["last_reset"], "W", W, "一致" if ok else "不一致", flush=True)
    # (3) ES 评价路径：C++ 与 Python
    for n in ("m02-K1-03", "m03-K12-00", "m04-K1-01", "m06-K12-02"):
        I = load_train(n); S = S_max(I); W = 2 * S; _, txt = inst_paths(n, "train")
        d = json.load(open(os.path.join(ROOT, "results", "T6", "es", f"{n}__es_lpgplan__s0.json")))
        lz = lpz_of(I, f"{T3R}/{n}.lpg.txt"); pf = ref.parse_planf(I, f"{T3R}/{n}__fresh_rlp.planf")
        for c in range(4):
            Nk = [x + rng.randint(0, 1) for x in d["Nk"]]; Mj = [x + rng.randint(0, 1) for x in d["Mj"]]
            ev = ["|"] + d["w"] + ["|"] + es_opts(n) + ["curdef", f"planW={W}"]
            ce = c_json(new, ["eval", txt, 3000, 0] + Nk + ["|"] + Mj + ev); cc = c_json(new, ["eval", txt, 20000, 1] + Nk + ["|"] + Mj + ev)
            pol = dict(w=[float(x) for x in d["w"]], curdef=True, lpz=lz, planf=pf, plan_until=W)
            pe = ref.run(I, Nk, Mj, 3000, False, **pol); pc = ref.run(I, Nk, Mj, 20000, True, **pol)
            a = dict(gaps=ce["gaps"], resets=ce["resets"], cgaps=cc["gaps"], cert=cc["cert"], period=cc["period"])
            p = dict(gaps=pe["gaps"], resets=pe["resets"], cgaps=pc["gaps"], cert=pc["cert"], period=pc["period"])
            ok = a == p and ce["last_reset"] < W; bad += (not ok)
            out.append(dict(part="es_path_py", inst=n, case=c, cpp=a, py=p, last_reset=ce["last_reset"], W=W, same=int(ok)))
            print("ES 路径", n, c, "一致" if ok else "不一致", a, "" if ok else p, flush=True)
    os.makedirs(os.path.join(ROOT, "results", "T12b"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "T12b", "check.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(json.dumps(dict(cases=len(out), mismatch=bad)))


if __name__ == "__main__":
    main()
