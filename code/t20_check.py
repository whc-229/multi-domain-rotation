"""T20-0 一致性核对（rost5 / rost5_ref 11-06 的改动：第一次缺岗 fg、共用外壳的补救 rescue、定台数限时 sizetl）。只用训练实例。
(A) 新路径：C++（新）与 Python 参考实现逐字段比较，含新字段 fg（第一次缺岗的时刻、站、原因）。
    8 条解码路径（第五批用到的全部）：pb、esdw（PB 权重加扰动 + curdef）、gphh（随机 GP 树）、rk、d5（带回溯修复）、
    rhmip（随机计划 planx + PB）、rddhand（PB + RDD 基因型 + resetopt）、ours（接入构造 + planx + devfirst + RDD 基因型；
    接入构造用 code/t20_ours.py 在训练实例上新做的严格窗口版本，产物在 results/T20/check_train/pre/）；
    2 个实例 × 2 组 = 每条路径 4 组，共 32 组。台数：第 0 组为默认台数 + 0..1，第 1 组随机减 1–3 个资源（让缺岗出现，fg 才有内容）。
    每组两次运行：3000 格不判环（缺岗、提前大修次数、窗口决策数、fg）与 2 万格判环（缺岗、判环结果、周期、fg）。
(B) 旧路径：不加 rescue / sizetl 时新旧程序输出逐字段相同（eval 去掉新字段 fg；所有命令去掉含 sec 的用时字段与 traj 里的秒数）：
    eval（(A) 的 8 条路径各 1 组）、size（fifo、pb、gphh）、fd（rk、esdw，限 40 代）、rddfd（rddhand、ours，限 40 代）。
(C) 补救的功能检查（新程序加 rescue）：T13 里外壳失败的 3 个作业——
    FIFO、PB 在 m04-K12-00（4 档起点定不出 / 复核 2 次缺岗）与 D5 在 m04-K1-00（起点定不出，另加 sizetl）；
    补救后的台数再用 rost5 judge 独立判定，字段与 size 输出相同。
(D) 消融 shorteval 用的 finite 选项与 judge 的 rescue（新程序）：训练实例 m04-K1-00 上，ours 的规则 + finite，rddfd 限 40 代（360 格）；
    最终解用 judge（finite、360 格）重判逐字段相同；再按长期口径 judge，不可行时 judge 加 rescue，补救后的台数用 judge（不加 rescue）独立判定，字段相同。
输出 results/T20/check0.jsonl，打印不一致数。用法：python3 code/t20_check.py <旧 rost5> <新 rost5>"""
import os, sys, json, random, subprocess, tempfile, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import rost5_ref as ref  # noqa: E402
import t13_methods as TM  # noqa: E402
from t11_job import write_geno  # noqa: E402
from t13_check import rand_tree  # noqa: E402
import t19b_job as B  # noqa: E402
import t20_ours  # noqa: E402

INSTS = ["m02-K12-00", "m04-K1-00"]
PATHS = ["pb", "esdw", "gphh", "rk", "d5", "rhmip", "rddhand", "ours"]
TIME_KEYS = ("sec",)


def c_json(r5, args):
    """运行 rost5，解析一行 JSON。"""
    r = subprocess.run([r5] + [str(a) for a in args], capture_output=True, text=True, check=True)
    return json.loads(r.stdout)


def strip_time(d):
    """去掉用时字段（键名含 sec）和 traj 里的秒数，供新旧程序比较。"""
    out = {k: v for k, v in d.items() if "sec" not in k}
    if "traj" in out: out["traj"] = [x[:3] + x[4:] for x in out["traj"]]
    return out


def base_fleet(r5, n, I, S, txt, meth, tf):
    """路径的默认台数（全 0 基因 / 基础权重用 rost5 size 定台数；ours 用接入构造台数）。"""
    W = 2 * S
    if meth == "ours":
        pre = t20_ours.prepare("check_train", n, "train"); return pre["conn"]["N"], pre["sigma_Mj"]
    if meth in ("pb", "rhmip", "rddhand"): o = TM.opts("pb", n, S)
    elif meth == "esdw": o = [TM.wstr(TM.W_PB)] + TM.opts("esdw", n, S)
    elif meth == "gphh": open(tf, "w").write("(- (* 10 flex) deficit)"); o = TM.opts("gphh", n, S, tf)
    else: o = TM.opts("rk", n, S) + ["pd=rk", f"W={W}"]     # d5 借用 rk 定出的台数（D5 在 m04-K1-00 上起点定不出，T13）
    b = c_json(r5, ["size", txt, 365, 20000, 20000] + o)
    assert b["R"] >= 0, (meth, n)
    return b["Nk"], b["Mj"]


def case(rng, r5, n, I, S, txt, meth, Nk, Mj, tmpd, tf):
    """生成一组随机输入：返回 (Python 参数, C++ eval 第 3 段权重, 第 4 段选项)。"""
    W = 2 * S; pol = {}; ev_w = []; ev_o = []
    if meth == "pb":
        pol["w"] = [float(x) for x in TM.W_PB] + [0.0] * 6; ev_w = TM.W_PB
    elif meth == "esdw":
        w = [x + 0.5 * rng.gauss(0, 1) * (abs(x) + 0.2) for x in TM.W_PB]
        pol.update(w=w + [0.0] * 6, curdef=True); ev_w = w; ev_o = ["curdef"]
    elif meth == "gphh":
        tr = rand_tree(rng); open(tf, "w").write(tr)
        pol.update(w=[float(x) for x in TM.W_TOOL] + [0.0] * 6, tree=ref.parse_tree(tr)); ev_w = TM.W_TOOL; ev_o = [f"tree=@{tf}"]
    elif meth in ("rk", "d5"):
        nd = c_json(r5, ["eval", txt, W, 0] + Nk + ["|"] + Mj + ["|", "|", "pd=rk", f"pdW={W}"])["ndec_w"]
        g = [rng.random() for _ in range(nd)]
        pol.update(w=[float(x) for x in TM.W_TOOL] + [0.0] * 6, curdef=True, pd=dict(mode=1 if meth == "rk" else 2, W=W, g=g, repair=(meth == "d5")))
        ev_w = TM.W_TOOL; ev_o = [f"pd={meth}", f"pdW={W}", "pdg=" + ",".join(repr(x) for x in g)]
    elif meth == "rhmip":
        pf = os.path.join(tmpd, "p.planx")
        with open(pf, "w") as f:
            f.write(f"end {4 * S}\n")
            for D in range(4 * S):
                for s2 in range(len(I.sites)):
                    if rng.random() < 0.3:
                        for _ in range(rng.randint(1, 2)):
                            k = rng.choice([k for k in range(len(I.T)) if s2 in I.elig[k]]); i = rng.randrange(I.L[k]); r = int(D < W and rng.random() < 0.2)
                            f.write(f"{s2} {D} {k} {i} {r}\n")
        pol.update(w=[float(x) for x in TM.W_PB] + [0.0] * 6, planx=ref.parse_planx(pf), plan_until=W); ev_w = TM.W_PB; ev_o = [f"planx={pf}", f"planW={W}"]
    elif meth == "rddhand":
        nd = c_json(r5, ["eval", txt, W, 0] + Nk + ["|"] + Mj + ["|"] + list(TM.W_PB) + ["|", f"rdd={_alloff(tmpd, W)}"])["ndec_w"]
        genes = [rng.random() for _ in range(2 * max(1, nd // 4))]; gp = os.path.join(tmpd, "g.geno"); write_geno(gp, W, nd, 0.4, genes, 1)
        pol.update(w=[float(x) for x in TM.W_PB] + [0.0] * 6, rdd=ref.decode_geno(genes, W, nd, 0.4, True)); ev_w = TM.W_PB; ev_o = [f"rdd={gp}"]
    elif meth == "ours":
        t20_ours.prepare("check_train", n, "train"); sp = t20_ours.spec("check_train", n, I, S, "ours"); Wd = int(B.MULT * S)
        o0 = [x for x in sp["opts"] if x != "devfirst"]
        nd = c_json(r5, ["eval", txt, Wd, 0] + Nk + ["|"] + Mj + ["|", "|"] + o0 + [f"rdd={_alloff(tmpd, Wd)}"])["ndec_w"]
        genes = [rng.random() for _ in range(2 * max(1, nd // 4))]; gp = os.path.join(tmpd, "g.geno"); write_geno(gp, Wd, nd, 0.4, genes, 1)
        pol.update(rdd=ref.decode_geno(genes, Wd, nd, 0.4, True), **sp["py"]); ev_o = sp["opts"] + [f"rdd={gp}"]
    return pol, ev_w, ev_o


def _alloff(tmpd, W):
    """全关的基因型（只为数窗口决策数）。"""
    gp = os.path.join(tmpd, "alloff.geno"); write_geno(gp, W, 0, 0.4, [], 1); return gp


def part_a(rng, new, tmpd, out):
    """新路径：C++（新）与 Python 逐字段（含 fg）。"""
    bad = 0; tf = os.path.join(tmpd, "tree.txt")
    for meth in PATHS:
        for n in INSTS:
            I = load_train(n); S = S_max(I); _, txt = inst_paths(n, "train")
            N0, M0 = base_fleet(new, n, I, S, txt, meth, tf)
            for c in range(2):
                Nk = [x + rng.randint(0, 1) for x in N0]; Mj = [x + rng.randint(0, 1) for x in M0]
                if c == 1:                     # 随机减 1–3 个资源，让缺岗出现
                    for _ in range(rng.randint(1, 3)):
                        pos = [(0, k) for k in range(len(Nk)) if Nk[k] > 0] + [(1, j) for j in range(len(Mj)) if Mj[j] > 0]
                        a, b = rng.choice(pos)
                        if a == 0: Nk[b] -= 1
                        else: Mj[b] -= 1
                pol, ev_w, ev_o = case(rng, new, n, I, S, txt, meth, Nk, Mj, tmpd, tf)
                diff = []; rec = dict(part="new_path", method=meth, inst=n, case=c, Nk=Nk, Mj=Mj)
                for Tm, cert in ((3000, False), (20000, True)):
                    py = ref.run(I, Nk, Mj, Tm, cert, **pol)
                    cc = c_json(new, ["eval", txt, Tm, int(cert)] + Nk + ["|"] + Mj + ["|"] + list(ev_w) + ["|"] + ev_o)
                    keys = ("gaps", "cert", "period", "fg") if cert else ("gaps", "resets", "ndec_w", "fg")
                    diff += [(Tm, k, py[k], cc[k]) for k in keys if py[k] != cc[k]]
                    rec[f"c{Tm}"] = {k: cc[k] for k in keys}
                rec.update(same=int(not diff), diff=diff); bad += bool(diff); out.append(rec)
                print("新路径", meth, n, c, "一致" if not diff else diff, rec["c20000"], flush=True)
    return bad


def part_b(rng, old, new, tmpd, out):
    """旧路径：不加新选项时新旧程序输出相同。"""
    bad = 0; tf = os.path.join(tmpd, "tree.txt")
    n = "m04-K1-00"; I = load_train(n); S = S_max(I); W = 2 * S; _, txt = inst_paths(n, "train")
    # eval：8 条路径各 1 组
    for meth in PATHS:
        N0, M0 = base_fleet(new, n, I, S, txt, meth, tf)
        Nk = [x + rng.randint(-1, 1) for x in N0]; Mj = [max(0, x + rng.randint(-1, 1)) for x in M0]; Nk = [max(0, x) for x in Nk]
        _, ev_w, ev_o = case(rng, new, n, I, S, txt, meth, Nk, Mj, tmpd, tf)
        for Tm, cert in ((3000, 0), (60000, 1)):
            a = ["eval", txt, Tm, cert] + Nk + ["|"] + Mj + ["|"] + list(ev_w) + ["|"] + ev_o
            eo = c_json(old, a); en = c_json(new, a); en.pop("fg")
            same = eo == en; bad += not same
            out.append(dict(part="old_path", cmd="eval", method=meth, T=Tm, same=int(same), diff=None if same else [eo, en]))
            print("旧路径 eval", meth, Tm, "相同" if same else (eo, en), flush=True)
    # size / fd / rddfd
    open(tf, "w").write("(- (* 10 flex) deficit)")
    runs = [("size", "fifo", ["size", txt, 365, 60000, 600000] + TM.opts("fifo", n, S)),
            ("size", "pb", ["size", txt, 365, 60000, 600000] + TM.opts("pb", n, S)),
            ("size", "gphh", ["size", txt, 365, 60000, 600000] + TM.opts("gphh", n, S, tf))]
    fdo = [f"W={W}", "stall=30", "swap=0.5", "fdsel=best"]
    for meth in ("rk", "esdw"):
        runs.append(("fd", meth, ["fd", txt, 365, 60000, 600000, 6 * S, 40, 10 ** 6, 4, 12, 0, TM.SEARCH[meth], 0.1, 0.4] + TM.opts(meth, n, S) + fdo))
    runs.append(("rddfd", "rddhand", ["rddfd", txt, 365, 60000, 600000, 6 * S, 40, 10 ** 6, 4, 12, 0, 0.1, 0.4] + TM.opts("rddhand", n, S) + fdo + ["resetopt"]))
    t20_ours.prepare("check_train", n, "train"); sp = t20_ours.spec("check_train", n, I, S, "ours")
    runs.append(("rddfd", "ours", ["rddfd", txt, 365, 60000, 600000, 6 * S, 40, 10 ** 6, 4, 12, 0, 0.1, 0.4] + sp["opts"] +
                 [f"W={int(B.MULT * S)}", "resetopt", "stall=30", "swap=0.5", "fdsel=best"]))
    for cmd, meth, a in runs:
        eo = strip_time(c_json(old, a)); en = strip_time(c_json(new, a))
        same = eo == en; bad += not same
        out.append(dict(part="old_path", cmd=cmd, method=meth, R=en.get("R"), same=int(same),
                        diff=None if same else {k: (eo.get(k), en.get(k)) for k in set(eo) | set(en) if eo.get(k) != en.get(k)}))
        print("旧路径", cmd, meth, "R", en.get("R"), "相同" if same else out[-1]["diff"], flush=True)
    return bad


def part_c(new, out):
    """补救的功能检查：T13 外壳失败的 3 个作业，加 rescue 后的结果；补救后的台数用 judge 独立判定。"""
    bad = 0
    cases = [("m04-K12-00", "fifo", []), ("m04-K12-00", "pb", []), ("m04-K1-00", "d5", ["sizetl=20"])]
    for n, meth, extra in cases:
        I = load_train(n); S = S_max(I); W = 2 * S; _, txt = inst_paths(n, "train")
        o = TM.opts(meth, n, S) + ([f"pd={meth}", f"W={W}"] if meth == "d5" else [])
        t0 = time.time(); old_like = c_json(new, ["size", txt, 365, 60000, 600000] + o + (extra if meth == "d5" else []))
        t1 = time.time(); r = c_json(new, ["size", txt, 365, 60000, 600000] + o + ["rescue"] + extra); t2 = time.time()
        rec = dict(part="rescue", inst=n, method=meth, extra=extra, without=dict(R=old_like["R"], status=old_like["status"], verify_gaps=old_like["verify_gaps"], factor=old_like["start_factor"], sec=round(t1 - t0, 1)),
                   with_rescue={k: r.get(k) for k in ("R", "status", "Nk", "Mj", "start_factor", "verify_gaps", "verify_cert", "exact", "rescue_add", "size_timeouts")}, sec=round(t2 - t1, 1))
        if r["R"] >= 0:
            j = c_json(new, ["judge", txt, 60000, 600000, ",".join(map(str, r["Nk"])), ",".join(map(str, r["Mj"]))] + o)
            same = (j["R"] == r["R"] and j["status"] == r["status"] and j["verify_gaps"] == r["verify_gaps"] and j["verify_cert"] == r["verify_cert"]
                    and j["exact"] == r["exact"] and j["resets"] == r["resets"] and j["last_reset"] == r["last_reset"])
            rec.update(judge_feasible=j["feasible"], judge_same=int(same)); bad += not same
        out.append(rec); print("补救", n, meth, json.dumps(rec, ensure_ascii=False), flush=True)
    return bad


def part_d(new, tmpd, out):
    """finite 选项与 judge 的 rescue。"""
    n = "m04-K1-00"; I = load_train(n); S = S_max(I); _, txt = inst_paths(n, "train"); Wd = int(B.MULT * S)
    t20_ours.prepare("check_train", n, "train"); sp = t20_ours.spec("check_train", n, I, S, "ours")
    r = c_json(new, ["rddfd", txt, 360, 360, 360, 360, 40, 10 ** 6, 4, 12, 0, 0.1, 0.4] + sp["opts"] + ["rescue", f"W={Wd}", "resetopt", "stall=30", "swap=0.5", "fdsel=best", "finite"])
    gp = os.path.join(tmpd, "fin.geno"); write_geno(gp, r["W"], r["ndev"], r["qon"], r["genes"], 1)
    fl = [",".join(map(str, r["Nk"])), ",".join(map(str, r["Mj"]))]; tail = sp["opts"] + [f"W={Wd}", "resetopt", f"rdd={gp}"]
    j = c_json(new, ["judge", txt, 360, 360] + fl + tail + ["finite"])
    same1 = j["R"] == r["R"] and j["feasible"] == 1 and j["status"] == r["status"] == "finite" and j["verify_gaps"] == r["verify_gaps"] == 0
    L = c_json(new, ["judge", txt, 60000, 600000] + fl + tail)
    rec = dict(part="finite", inst=n, R_finite=r["R"], base_R=r["base_R"], gens=r["gens"], judge_finite_same=int(same1), long_feasible=L["feasible"], long=L)
    bad = int(not same1)
    if L["feasible"] != 1:
        Rr = c_json(new, ["judge", txt, 60000, 600000] + fl + tail + ["rescue"])
        J2 = c_json(new, ["judge", txt, 60000, 600000, ",".join(map(str, Rr["Nk"])), ",".join(map(str, Rr["Mj"]))] + tail)
        same2 = all(J2[k] == Rr[k] for k in ("R", "feasible", "status", "verify_gaps", "verify_cert", "verify_period", "exact", "resets", "last_reset"))
        rec.update(rescue=Rr, rescue_judge_same=int(same2)); bad += not same2
    out.append(rec); print("finite", json.dumps(rec, ensure_ascii=False), flush=True)
    return bad


def main():
    """(A) 新路径、(B) 旧路径、(C) 补救功能检查、(D) finite 与 judge 的 rescue。"""
    old, new = sys.argv[1], sys.argv[2]
    rng = random.Random(20261106); out = []
    with tempfile.TemporaryDirectory(prefix="t20chk_") as tmpd:
        ba = part_a(rng, new, tmpd, out)
        bb = part_b(rng, old, new, tmpd, out)
        bc = part_c(new, out)
        bd = part_d(new, tmpd, out)
    os.makedirs(os.path.join(ROOT, "results", "T20"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "T20", "check0.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(json.dumps(dict(new_path=sum(r["part"] == "new_path" for r in out), new_path_mismatch=ba,
                          old_path=sum(r["part"] == "old_path" for r in out), old_path_mismatch=bb, rescue_cases=3, rescue_judge_mismatch=bc, finite_mismatch=bd)))


if __name__ == "__main__":
    main()
