"""T20 补足阶段停止条件（11-06 审核）的核对。只用训练实例。
改动：rost5 加 rescue 时，补足至少补 60 次；之后只要 365 格缺岗在最近 5 次补台里刷新过最低值就继续，否则本档失败。size 命令另加 hardtl（总限时）。
(F1) C++ 与 Python 参考实现（rost5_ref.size_res, fillpat=True）逐步相同：
     每次补足加的资源下标与加后的缺岗（C++ sizelog 的 "fill add=… gaps=…" 行），以及 size_res 的最终结果（C++ 未补救时）。
     为了走到"越过下限后按缺岗是否刷新决定继续"的分支，用只为核对的 filladd=3、10 把下限调小（Python 同样 max_add），另有默认 60 一组；
     起点 = PB 定出的台数 × 0.4（向下取整），规则：PB、2 条随机 GP 树；实例 m02-K12-00、m04-K1-00。
(F2) 新旧程序都加 rescue 时，补足在 60 次以内的（新程序输出 fill_max ≤ 60）输出逐字段相同（去掉用时与新字段 fill_max）：
     size（FIFO、PB、GPHH、σ 定台数、只用接入构造）；fd（RK、ES-DW、D5，10 代）；rddfd（RDD-hand、本方法，10 代）；实例 m04-K1-00、m04-K12-00。
(F3) 训练实例重测 报告/T20_0_冒烟.md 第 5 节的表：PB 与 RK 默认排序加 rescue 定台数，96 个训练实例 × 2，旧程序与新程序各跑一次；
     新程序加 hardtl = 该 m 的预算。fill_max ≤ 60 的必须与旧程序逐字段相同；汇总"定不出"的变化。
输出 results/T20/check_fill.jsonl（F1、F2、F3 的比较），results/T20/probe/fill2_old.jsonl、fill2_new.jsonl（F3 的原始结果）。
用法：python3 code/t20_check_fill.py <旧 rost5> <新 rost5> [并行数]"""
import os, sys, json, random, re, subprocess, tempfile
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths, train_names  # noqa: E402
from freshlb import S_max  # noqa: E402
import rost5_ref as ref  # noqa: E402
import t13_methods as TM  # noqa: E402
from t13_check import rand_tree  # noqa: E402
import t19b_job as B  # noqa: E402
import t20_ours  # noqa: E402
from t20_check import c_json, strip_time  # noqa: E402
from t20_job import TLIM  # noqa: E402


def parse_sizelog(path):
    """C++ sizelog → (补足轨迹 [(下标, 缺岗)], 补足后台数, 回删后台数)。"""
    tr, filled, pruned = [], None, None
    for ln in open(path):
        m = re.match(r"fill add=(\d+) gaps=(\d+)", ln)
        if m: tr.append((int(m.group(1)), int(m.group(2)))); continue
        if ln.startswith("filled"): filled = ln.split(None, 1)[1].strip()
        if ln.startswith("pruned"): pruned = ln.split(None, 1)[1].strip()
    return tr, filled, pruned


def f1(rng, new, tmpd, out):
    """C++ 与 Python 的补足 / 回删 / 最终结果逐步比较。"""
    bad = 0
    for n in ("m02-K12-00", "m04-K1-00"):
        I = load_train(n); S = S_max(I); _, txt = inst_paths(n, "train")
        b = c_json(new, ["size", txt, 365, 20000, 20000] + TM.opts("pb", n, S))
        N0 = [int(x * 0.4) for x in b["Nk"]]; M0 = [int(x * 0.4) for x in b["Mj"]]
        rules = [("pb", [o for o in TM.opts("pb", n, S) if o != "facloop"], dict(w=[float(x) for x in TM.W_PB] + [0.0] * 6))]
        for t in range(2):
            tr = rand_tree(rng); tf = os.path.join(tmpd, f"tree{t}.txt"); open(tf, "w").write(tr)
            rules.append((f"gp{t}", [TM.wstr(TM.W_TOOL), f"tree=@{tf}"], dict(w=[float(x) for x in TM.W_TOOL] + [0.0] * 6, tree=ref.parse_tree(tr))))
        for rname, opts, pol in rules:
            for fa in (3, 10, 60):
                lg = os.path.join(tmpd, "size.log")
                cc = c_json(new, ["size", txt, 365, 20000, 20000] + opts + ["N0=" + ",".join(map(str, N0)), "M0=" + ",".join(map(str, M0)),
                                                                           "rescue", f"filladd={fa}", f"sizelog={lg}"])
                ctr, cfill, cprune = parse_sizelog(lg)
                ptr = []
                R, Nk, Mj, st, per = ref.size_res(I, pol, 365, 20000, N0, M0, max_add=fa, fillpat=True, trace=ptr)
                diff = []
                if ctr != ptr: diff.append(("trace", len(ctr), len(ptr)))
                if cc["rescue_add"] == 0 and (cc["R"], cc["Nk"], cc["Mj"], cc["status"]) != (R, Nk, Mj, st): diff.append(("final", cc["R"], R))
                # 补足的结局：filled = 补到 0 缺岗；stall = 越过下限后最近 5 次没刷新最低值而停（只为说明覆盖了哪个分支）
                stop = "filled" if cfill else "stall"
                rec = dict(part="F1", inst=n, rule=rname, filladd=fa, adds=len(ctr), stop=stop, R_c=cc["R"], R_py=R, same=int(not diff), diff=diff)
                bad += bool(diff); out.append(rec); print("F1", json.dumps(rec, ensure_ascii=False), flush=True)
    return bad


def f2(old, new, tmpd, out):
    """新旧程序都加 rescue：fill_max ≤ 60 时逐字段相同。"""
    bad = 0; tf = os.path.join(tmpd, "tree.txt"); open(tf, "w").write("(- (* 10 flex) deficit)")
    for n in ("m04-K1-00", "m04-K12-00"):
        I = load_train(n); S = S_max(I); W = 2 * S; _, txt = inst_paths(n, "train")
        pre = t20_ours.prepare("check_train", n, "train"); art = os.path.join(ROOT, "results", "T20", "check_train", "pre", n)
        sp_c = t20_ours.spec("check_train", n, I, S, "conn_only"); sp_o = t20_ours.spec("check_train", n, I, S, "ours")
        sz = ["size", txt, 365, 60000, 600000]
        runs = [("size", "fifo", sz + TM.opts("fifo", n, S)), ("size", "pb", sz + TM.opts("pb", n, S)), ("size", "gphh", sz + TM.opts("gphh", n, S, tf)),
                ("size", "sigma", sz + [f"rlp={os.path.join(art, 'steady.rlp')}"]), ("size", "conn_only", sz + sp_c["opts"])]
        fdo = [f"W={W}", "stall=30", "swap=0.5", "fdsel=best"]
        for meth in ("rk", "esdw", "d5"):
            runs.append(("fd", meth, ["fd", txt, 365, 60000, 600000, 6 * S, 10, 10 ** 6, 4, 12, 0, TM.SEARCH[meth], 0.1, 0.4] + TM.opts(meth, n, S) + fdo))
        runs.append(("rddfd", "rddhand", ["rddfd", txt, 365, 60000, 600000, 6 * S, 10, 10 ** 6, 4, 12, 0, 0.1, 0.4] + TM.opts("rddhand", n, S) + fdo + ["resetopt"]))
        runs.append(("rddfd", "ours", ["rddfd", txt, 365, 60000, 600000, 6 * S, 10, 10 ** 6, 4, 12, 0, 0.1, 0.4] + sp_o["opts"] +
                     [f"W={int(B.MULT * S)}", "resetopt", "stall=30", "swap=0.5", "fdsel=best"]))
        for cmd, meth, a in runs:
            a = a + ["rescue"]
            eo = strip_time(c_json(old, a)); en = strip_time(c_json(new, a)); fm = en.pop("fill_max")
            same = eo == en; need = fm <= 60; bad += need and not same
            rec = dict(part="F2", inst=n, cmd=cmd, method=meth, R_old=eo.get("R"), R_new=en.get("R"), fill_max=fm, same=int(same),
                       diff=None if same else {k: (eo.get(k), en.get(k)) for k in set(eo) | set(en) if eo.get(k) != en.get(k)})
            out.append(rec); print("F2", json.dumps(rec, ensure_ascii=False)[:400], flush=True)
    return bad


def probe_one(a):
    """F3 的一个 (训练实例, 方法, 程序)：rost5 size 加 rescue（新程序另加 hardtl = 预算）。"""
    r5, tag, n, meth = a
    I = load_train(n); S = S_max(I); _, txt = inst_paths(n, "train")
    o = TM.opts(meth, n, S) + (["pd=rk", f"W={2 * S}"] if meth == "rk" else []) + ["rescue"] + ([f"hardtl={TLIM[I.meta['m']]}"] if tag == "new" else [])
    r = c_json(r5, ["size", txt, "365", "60000", "600000"] + o)
    return tag, dict(inst=n, method=meth, **r)


def f3(old, new, out, P):
    """训练实例重测：旧 / 新各跑一次，fill_max ≤ 60 的逐字段相同；写原始结果。"""
    jobs = [(r5, tag, n, m) for n in train_names([2, 3, 4, 6, 8, 10, 15, 20]) for m in ("pb", "rk") for r5, tag in ((old, "old"), (new, "new"))]
    res = {"old": {}, "new": {}}
    d = os.path.join(ROOT, "results", "T20", "probe")
    with open(os.path.join(d, "fill2_old.jsonl"), "w") as fo, open(os.path.join(d, "fill2_new.jsonl"), "w") as fn, ThreadPoolExecutor(P) as ex:
        for tag, r in ex.map(probe_one, jobs):
            res[tag][(r["inst"], r["method"])] = r
            (fo if tag == "old" else fn).write(json.dumps(r, ensure_ascii=False) + "\n"); (fo if tag == "old" else fn).flush()
            print("F3", tag, r["inst"], r["method"], r["R"], r["status"], r.get("fill_max"), r["sec"], flush=True)
    bad = 0
    for key in res["new"]:
        eo = strip_time(res["old"][key]); en = strip_time(res["new"][key]); fm = en.pop("fill_max"); en.pop("size_timeouts", None); eo.pop("size_timeouts", None)
        same = eo == en; need = fm <= 60; bad += need and not same
        out.append(dict(part="F3", inst=key[0], method=key[1], R_old=eo["R"], R_new=en["R"], fill_max=fm,
                        timeouts=res["new"][key].get("size_timeouts"), same=int(same), need_same=int(need)))
    return bad


def main():
    """F1、F2、F3，写 results/T20/check_fill.jsonl，打印各部分的不一致数。"""
    old, new = sys.argv[1], sys.argv[2]; P = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    rng = random.Random(20261106); out = []
    with tempfile.TemporaryDirectory(prefix="t20fill_") as tmpd:
        b1 = f1(rng, new, tmpd, out)
        b2 = f2(old, new, tmpd, out)
    b3 = f3(old, new, out, P)
    with open(os.path.join(ROOT, "results", "T20", "check_fill.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(json.dumps(dict(F1=sum(r["part"] == "F1" for r in out), F1_mismatch=b1, F2=sum(r["part"] == "F2" for r in out), F2_mismatch_within60=b2,
                          F3=sum(r["part"] == "F3" for r in out), F3_mismatch_within60=b3)))


if __name__ == "__main__":
    main()
