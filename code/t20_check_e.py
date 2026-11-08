"""T20 (e)（搜索型方法定台数失败时的搜索起点，11-06 晚二审核）的核对与训练实例量测。只用训练实例。
(E1) 不触发 (e) 时与现装程序逐字段相同（去掉用时与新字段 e_fixed / e_gaps0 / e_found_gen；新程序必须报 e_fixed = −1）：
     size（FIFO、PB、GPHH、σ、只用接入构造）、fd（RK、ES-DW、D5，10 代）、rddfd（RDD-hand、本方法，10 代），都加 rescue；实例 m04-K1-00、m04-K12-00。
(E1b) 触发 (e) 的功能检查：D5 在 m04-K1-00 上加 sizetl=20 起点定不出（T13、T20-0 已知），fd 限时 300 秒；
     (e) 必须触发（e_fixed ≥ 0）；找到解时最终台数用 rost5 judge 独立判定，字段相同且可行。
(E2) 量测（审核要求，不是装上的门槛）：报告/T20_0_冒烟.md 第 7.4 节里 PB / RK 定不出的 41 个 (训练实例, 起点)——
     RK 起点定不出的跑 RK-EA，PB 起点定不出的跑 RDD-hand（它的默认规则是 PB）；各 1 个种子（0），按 T20 预算与 T20 作业相同的参数；
     记 (e) 是否触发、能否在预算内找到可行解；找到时用 judge 独立判定。
输出 results/T20/check_e.jsonl（E1、E1b）与 results/T20/probe/e_train.jsonl（E2）。
用法：python3 code/t20_check_e.py <现装 rost5> <新 rost5> [E2 并行数，0 = 不跑 E2] [e2only]
11-09 修改：原 E1b 用例在新补足规则下起点已能定出（base_R = 69），(e) 不触发，作废（结果仍在 check_e.jsonl，ok = 0）；
     改由 E2 中预算最小的 3 个作业做触发检查（E1c，见报告），所以 E2 改为按预算从小到大派；加 e2only 参数只跑 E2。"""
import os, sys, json, subprocess, tempfile
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import t13_methods as TM  # noqa: E402
import t19b_job as B  # noqa: E402
import t20_ours  # noqa: E402
from t20_check import c_json, strip_time  # noqa: E402
from t20_job import TLIM  # noqa: E402

E_KEYS = ("e_fixed", "e_gaps0", "e_found_gen")


def judge_same(r5, txt, rec, opts):
    """最终台数用 judge 独立判定：返回 (可行, 与 fd 输出的判定字段是否相同)。"""
    if rec["R"] < 0: return None, None
    j = c_json(r5, ["judge", txt, 60000, 600000, ",".join(map(str, rec["Nk"])), ",".join(map(str, rec["Mj"]))] + opts)   # 同 t13_job.checks 的写法
    same = all(j[k] == rec[k] for k in ("R", "status", "verify_gaps", "verify_cert", "verify_period", "exact"))
    return int(j["feasible"]), int(same)


def e1(old, new, tmpd, out):
    """不触发 (e)：新旧逐字段相同。"""
    bad = 0; tf = os.path.join(tmpd, "tree.txt"); open(tf, "w").write("(- (* 10 flex) deficit)")
    for n in ("m04-K1-00", "m04-K12-00"):
        I = load_train(n); S = S_max(I); W = 2 * S; _, txt = inst_paths(n, "train")
        t20_ours.prepare("check_train", n, "train"); art = os.path.join(ROOT, "results", "T20", "check_train", "pre", n)
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
            eo = strip_time(c_json(old, a)); en = strip_time(c_json(new, a)); ev = {k: en.pop(k) for k in E_KEYS if k in en}
            trig = ev.get("e_fixed", -1) >= 0
            same = eo == en; ok = same and not trig; bad += not ok
            rec = dict(part="E1", inst=n, cmd=cmd, method=meth, R_old=eo.get("R"), R_new=en.get("R"), e=ev, same=int(same), ok=int(ok),
                       diff=None if same else {k: (eo.get(k), en.get(k)) for k in set(eo) | set(en) if eo.get(k) != en.get(k)})
            out.append(rec); print("E1", json.dumps(rec, ensure_ascii=False)[:300], flush=True)
    return bad


def e1b(new, out):
    """触发 (e) 的功能检查：D5 在 m04-K1-00 上起点定不出。"""
    n = "m04-K1-00"; I = load_train(n); S = S_max(I); W = 2 * S; _, txt = inst_paths(n, "train")
    o = TM.opts("d5", n, S) + ["rescue", "sizetl=20"]
    r = c_json(new, ["fd", txt, 365, 60000, 600000, 6 * S, 10 ** 9, 300, 4, 12, 0, TM.SEARCH["d5"], 0.1, 0.4] + o + [f"W={W}", "stall=30", "swap=0.5", "fdsel=best"])
    trig = r["e_fixed"] >= 0
    feas, same = judge_same(new, txt, r, [x for x in o if not x.startswith("sizetl=")] + [f"W={W}", "pd=d5", "pdg=" + ",".join(repr(x) for x in r["genes"])])
    ok = trig and (r["R"] < 0 or (feas == 1 and same == 1))
    rec = dict(part="E1b", inst=n, method="d5", base_R=r["base_R"], e_fixed=r["e_fixed"], e_gaps0=r["e_gaps0"], e_found_gen=r["e_found_gen"], R=r["R"], status=r["status"],
               gens=r["gens"], sec=r["sec"], judge_feasible=feas, judge_same=same, ok=int(ok))
    out.append(rec); print("E1b", json.dumps(rec, ensure_ascii=False), flush=True)
    return int(not ok)


def e2_one(a):
    """E2 的一个作业：训练实例上按 T20 预算跑 RK-EA（rk 起点定不出）或 RDD-hand（pb 起点定不出）。"""
    r5, n, start = a
    I = load_train(n); S = S_max(I); W = 2 * S; _, txt = inst_paths(n, "train"); TL = TLIM[I.meta["m"]]
    meth = "rk" if start == "rk" else "rddhand"
    opts = TM.opts(meth, n, S) + ["rescue"]
    fdo = [f"W={W}", "stall=30", "swap=0.5", "fdsel=best"]
    if meth == "rddhand":
        args = ["rddfd", txt, 365, 60000, 600000, 6 * S, 10 ** 9, TL, 4, 12, 0, 0.1, 0.4] + opts + fdo + ["resetopt"]
    else:
        args = ["fd", txt, 365, 60000, 600000, 6 * S, 10 ** 9, TL, 4, 12, 0, TM.SEARCH[meth], 0.1, 0.4] + opts + fdo
    r = c_json(r5, args)
    feas = same = None
    if r["R"] >= 0:
        with tempfile.TemporaryDirectory(prefix="t20e_") as d:
            if meth == "rddhand":
                from t11_job import write_geno
                gp = os.path.join(d, "g.geno"); write_geno(gp, r["W"], r["ndev"], r["qon"], r["genes"], 1)
                jo = opts + [f"W={W}", "resetopt", f"rdd={gp}"]
            else:
                jo = opts + [f"W={W}", "pd=rk", "pdg=" + ",".join(repr(x) for x in r["genes"])]
            feas, same = judge_same(r5, txt, r, jo)
    return dict(inst=n, start=start, method=meth, TL=TL, base_R=r["base_R"], f0_src=r["f0_src"], e_fixed=r["e_fixed"], e_gaps0=r["e_gaps0"],
                e_found_gen=r["e_found_gen"], R=r["R"], status=r["status"], gens=r["gens"], sec=r["sec"], judge_feasible=feas, judge_same=same)


def e2(new, P):
    """41 个 (训练实例, 起点)，按预算从小到大派（前 3 个兼作 (e) 的触发检查 E1c）。"""
    fails = [json.loads(l) for l in open(os.path.join(ROOT, "results", "T20", "probe", "fill2_new.jsonl"))]
    jobs = sorted([(new, r["inst"], r["method"]) for r in fails if r["R"] < 0], key=lambda a: (TLIM[int(a[1][1:3])], a[1], a[2]))
    print("E2", len(jobs), "个作业", flush=True)
    path = os.path.join(ROOT, "results", "T20", "probe", "e_train.jsonl")
    done = set()
    if os.path.exists(path): done = {(json.loads(l)["inst"], json.loads(l)["start"]) for l in open(path)}
    jobs = [a for a in jobs if (a[1], a[2]) not in done]
    with open(path, "a") as f, ThreadPoolExecutor(P) as ex:
        for r in ex.map(e2_one, jobs):
            f.write(json.dumps(r, ensure_ascii=False) + "\n"); f.flush(); print("E2", json.dumps(r, ensure_ascii=False), flush=True)


def main():
    """E1、E1b 写 check_e.jsonl 并打印不一致数；再跑 E2（可续跑）。"""
    old, new = sys.argv[1], sys.argv[2]; P = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    if len(sys.argv) > 4 and sys.argv[4] == "e2only":
        e2(new, P); print("E2 结束", flush=True); return
    out = []
    with tempfile.TemporaryDirectory(prefix="t20e_") as tmpd:
        b1 = e1(old, new, tmpd, out)
    b1b = e1b(new, out)
    with open(os.path.join(ROOT, "results", "T20", "check_e.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(json.dumps(dict(E1=sum(r["part"] == "E1" for r in out), E1_bad=b1, E1b_bad=b1b)), flush=True)
    if P > 0: e2(new, P)
    print("E2 结束", flush=True)


if __name__ == "__main__":
    main()
