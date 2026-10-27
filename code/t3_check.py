"""T3 一致性核对。
(a) rost5 vs Python 参考实现（code/rost5_ref.py）：R_LP 规则，(实例, 变体, 场景) × 若干台数，比较缺岗数、是否判出循环、周期长度；
    另加基础打分策略在 K12 实例上的比对（检查岗位化的代码路径）。
(b) rost5 vs rost4（新旧版本）：K1 实例上 eval（多组权重与选项）和 es（同参数、同种子）输出逐字段相同。
用法：python3 code/t3_check.py <a|b> <规则文件目录> <实例名,...> [组数]"""
import os, sys, json, random, subprocess
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, PRIOR, load_train, inst_paths  # noqa: E402
import rost5_ref as ref  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")
R4 = os.path.join(PRIOR, "cpp", "rost4")


def c_eval(exe, txt, Tmax, cert, Nk, Mj, w, opts):
    args = [exe, "eval", txt, str(Tmax), str(int(cert))] + [str(x) for x in Nk] + ["|"] + [str(x) for x in Mj] + ["|"] + [repr(x) for x in w] + (["|"] + opts if opts else [])
    out = subprocess.run(args, capture_output=True, text=True)
    return json.loads(out.stdout)


def check_a(rdir, names, ngroups, Tc=60000):
    """R_LP：每个 (实例, 变体, 场景) 取 ngroups 组台数（LP 取整 + 0..2 的随机增量）。"""
    rng = random.Random(7); tot = bad = 0
    for n in names:
        I = load_train(n); _, txt = inst_paths(n, "train")
        for tag in ["fresh_rlp", "fresh_noreset", "steady", "fresh_rlp_c", "fresh_noreset_c"]:
            rf = os.path.join(rdir, f"{n}__{tag}.rlp"); R = ref.parse_rlp(I, rf)
            pf = os.path.join(rdir, f"{n}__{tag}.planf") if tag.startswith("fresh_rlp") else None
            planf = ref.parse_planf(I, pf) if pf else None
            for g in range(ngroups):
                Nk = [x + rng.randint(0, 2) for x in R["N0"]]; Mj = [x + rng.randint(0, 2) for x in R["M0"]]
                for Tmax, cert in [(3000, False), (Tc, True)]:
                    py = ref.run(I, Nk, Mj, Tmax, cert, rlp=R, planf=planf)
                    cc = c_eval(R5, txt, Tmax, cert, Nk, Mj, [0] * 8, [f"rlp={rf}"] + ([f"planf={pf}"] if pf else []))
                    ok = py["gaps"] == cc["gaps"] and (not cert or (py["cert"] == cc["cert"] and py["period"] == cc["period"]
                                                                     and (py["cert"] != 1 or cc["exact"] == 1)))
                    ok = ok and (py["resets"] == cc["resets"] or py["cert"] == -1)
                    tot += 1; bad += (not ok)
                    if not ok: print("不一致", n, tag, g, Tmax, cert, py, cc, flush=True)
    return tot, bad


def check_a_score(names, ngroups, Tc=20000):
    """基础打分策略（无选项）在给定实例上的 rost5 vs 参考实现。"""
    rng = np.random.default_rng(3); r2 = random.Random(3); tot = bad = 0
    from tools_model import tool_lb
    from mdrr import ToolInst
    for n in names:
        I = load_train(n); _, txt = inst_paths(n, "train")
        _, xN, xM = tool_lb(ToolInst(I), integer=True)
        for g in range(ngroups):
            Nk = [int(round(v)) + r2.randint(0, 3) for v in xN]; Mj = [int(round(v)) + r2.randint(0, 3) for v in xM]
            w = rng.normal(0, 1, 8).tolist() if g else [0.1, 0, 10, 1, 0, 0, 1, 0.1]
            for Tmax, cert in [(3000, False), (Tc, True)]:
                py = ref.run(I, Nk, Mj, Tmax, cert, w=w)
                cc = c_eval(R5, txt, Tmax, cert, Nk, Mj, w, [])
                ok = py["gaps"] == cc["gaps"] and (not cert or (py["cert"] == cc["cert"] and py["period"] == cc["period"]))
                tot += 1; bad += (not ok)
                if not ok: print("不一致(打分)", n, g, Tmax, cert, py, cc, flush=True)
    return tot, bad


def check_b(names, rdir, ngroups=6):
    """rost5 vs rost4：K1 实例。eval 用 13 维随机权重与多种选项；输出逐字段比较。"""
    rng = np.random.default_rng(11); r2 = random.Random(11); tot = bad = 0
    from tools_model import tool_lb
    from mdrr import ToolInst
    for n in names:
        I = load_train(n); assert I.meta["family"] == "K1"; _, txt = inst_paths(n, "train")
        _, xN, xM = tool_lb(ToolInst(I), integer=True)
        pf = os.path.join(rdir, f"{n}__fresh_rlp.planf")
        optsets = [[], ["curdef"], ["reset", "Tw=300"], ["reset", "hard", "Tw=200"], ["stagger"]]
        if os.path.exists(pf): optsets.append([f"planf={pf}"])
        for g in range(ngroups):
            Nk = [int(round(v)) + r2.randint(0, 3) for v in xN]; Mj = [int(round(v)) + r2.randint(0, 3) for v in xM]
            w = rng.normal(0, 1, 13).tolist()
            for opts in optsets:
                for Tmax, cert in [(3000, 0), (60000, 1)]:
                    a = c_eval(R4, txt, Tmax, cert, Nk, Mj, w, opts); b = c_eval(R5, txt, Tmax, cert, Nk, Mj, w, opts)
                    tot += 1; bad += (a != b)
                    if a != b: print("新旧不一致 eval", n, opts, a, b, flush=True)
    return tot, bad


def check_b_es(names, rdir, seeds=(0, 1)):
    """es：curdef 与 lpgplan（lpg 文件由 planf 同源生成），小代数，输出除 sec 外逐字段相同。"""
    tot = bad = 0
    for n in names:
        _, txt = inst_paths(n, "train")
        cfgs = [["curdef"]]
        lpg = os.path.join(rdir, f"{n}.lpg.txt"); pf = os.path.join(rdir, f"{n}__fresh_rlp.planf")
        if os.path.exists(lpg): cfgs.append([f"lpg={lpg}", f"planf={pf}"])
        for sd in seeds:
            for opts in cfgs:
                args = ["es", txt, "certs", "365", "3", "4", "12", str(sd), "60000", "600000"] + opts
                a = json.loads(subprocess.run([R4] + args, capture_output=True, text=True).stdout)
                b = json.loads(subprocess.run([R5] + args, capture_output=True, text=True).stdout)
                a.pop("sec"); b.pop("sec")
                tot += 1; bad += (a != b)
                if a != b: print("新旧不一致 es", n, opts, sd, a, b, flush=True)
    return tot, bad


if __name__ == "__main__":
    mode, rdir, names = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
    ng = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    if mode == "a":
        print(json.dumps(dict(check="rlp_vs_ref", groups_mismatch=check_a(rdir, names, ng))), flush=True)
    elif mode == "ascore":
        print(json.dumps(dict(check="score_vs_ref", groups_mismatch=check_a_score(names, ng))), flush=True)
    elif mode == "b":
        print(json.dumps(dict(check="rost5_vs_rost4_eval", groups_mismatch=check_b(names, rdir, ng))), flush=True)
    elif mode == "bes":
        print(json.dumps(dict(check="rost5_vs_rost4_es", groups_mismatch=check_b_es(names, rdir))), flush=True)
