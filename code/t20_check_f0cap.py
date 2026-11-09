"""T20 后备上限（11-09 下午审核 (b)：rescue 时原有后备最多用预算的一半，剩下的留给 (e)；实现为"后备开始时剩余预算的一半"）的核对。只用训练实例。
(G1) 不触发 (e)、后备不跑时与现装程序逐字段相同：同 code/t20_check_e.py 的 E1（20 组，去掉用时与新字段；新程序必须报 e_fixed = −1、f0_capped = 0）。
(G2) 后备在一半预算内定出的情形：E2 里 A 类（后备定出）的 5 个训练作业，按 T20 预算设限时、跑 10 代（结果只由代数决定），
     新旧程序逐字段相同；新程序 f0_capped 必须为 0（若为 1，说明这一例后备超过了一半预算，单独记下，不算不一致）。
(G3) 上限起作用的情形：E2 里 D 类（后备用完预算、(e) 没触发）的 m10-K12-03 RK-EA，按 T20 预算完整跑；
     必须 f0_capped = 1、后备用时 ≤ 一半预算 + 60 秒（剩余的一半更小，这是宽松的上界）、(e) 触发（e_fixed ≥ 0）、总用时 ≤ 预算 + 60 秒；找到解时用 judge 独立判定。
输出 results/T20/check_f0cap.jsonl，最后一行是汇总 {"G1_bad", "G2_bad", "G3_bad"}。
用法：python3 code/t20_check_f0cap.py <现装 rost5> <新 rost5>"""
import os, sys, json, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import t13_methods as TM  # noqa: E402
import t20_check_e as CE  # noqa: E402
from t20_check import c_json, strip_time  # noqa: E402
from t20_job import TLIM  # noqa: E402

NEW_KEYS = ("e_fixed", "e_gaps0", "e_found_gen", "f0_capped")
G2_CASES = [("m06-K12-05", "rk"), ("m08-K12-04", "rk"), ("m10-K1-00", "rk"), ("m15-K1-04", "rddhand"), ("m15-K12-00", "rddhand")]
G3_CASE = ("m10-K12-03", "rk")


def fd_args(meth, n, gens, tl):
    """与 E2 / T20 作业相同的 fd / rddfd 参数（种子 0，加 rescue）；gens、限时可指定。"""
    I = load_train(n); S = S_max(I); W = 2 * S; _, txt = inst_paths(n, "train")
    opts = TM.opts(meth, n, S) + ["rescue"]; fdo = [f"W={W}", "stall=30", "swap=0.5", "fdsel=best"]
    if meth == "rddhand":
        return txt, opts, W, ["rddfd", txt, 365, 60000, 600000, 6 * S, gens, tl, 4, 12, 0, 0.1, 0.4] + opts + fdo + ["resetopt"]
    return txt, opts, W, ["fd", txt, 365, 60000, 600000, 6 * S, gens, tl, 4, 12, 0, TM.SEARCH[meth], 0.1, 0.4] + opts + fdo


def g1(old, new, out):
    """E1 原样重做；新字段多了 f0_capped，必须为 0。"""
    CE.E_KEYS = NEW_KEYS
    with tempfile.TemporaryDirectory(prefix="t20g_") as tmpd:
        bad = CE.e1(old, new, tmpd, out)
    for r in out:
        if r["part"] == "E1":
            r["part"] = "G1"
            if r["e"].get("f0_capped", 0) != 0: r["ok"] = 0; bad += 1
    return bad


def g2(old, new, out):
    """后备在一半预算内定出：新旧逐字段相同。"""
    bad = 0
    for n, meth in G2_CASES:
        TL = TLIM[int(n[1:3])]; _, _, _, a = fd_args(meth, n, 10, TL)
        eo = strip_time(c_json(old, a)); en = strip_time(c_json(new, a)); ev = {k: en.pop(k) for k in NEW_KEYS if k in en}
        for k in ("e_fixed", "e_gaps0", "e_found_gen"): eo.pop(k, None)
        capped = ev.get("f0_capped", 0) == 1; same = eo == en
        ok = same or capped; bad += not ok
        rec = dict(part="G2", inst=n, method=meth, TL=TL, f0_src_old=eo.get("f0_src"), f0_src_new=en.get("f0_src"), R_old=eo.get("R"), R_new=en.get("R"),
                   e=ev, same=int(same), capped=int(capped), ok=int(ok),
                   diff=None if same else {k: (eo.get(k), en.get(k)) for k in set(eo) | set(en) if eo.get(k) != en.get(k)})
        out.append(rec); print("G2", json.dumps(rec, ensure_ascii=False)[:400], flush=True)
    return bad


def g3(new, out):
    """上限起作用：后备停在一半预算，(e) 触发；找到解时 judge。"""
    n, meth = G3_CASE; TL = TLIM[int(n[1:3])]; txt, opts, W, a = fd_args(meth, n, 10 ** 9, TL)
    r = c_json(new, a)
    feas = same = None
    if r["R"] >= 0:
        feas, same = CE.judge_same(new, txt, r, opts + [f"W={W}", "pd=rk", "pdg=" + ",".join(repr(x) for x in r["genes"])])
    ok = (r["f0_capped"] == 1 and r["sec_f0"] <= 0.5 * TL + 60 and r["e_fixed"] >= 0 and r["sec"] <= TL + 60
          and (r["R"] < 0 or (feas == 1 and same == 1)))
    rec = dict(part="G3", inst=n, method=meth, TL=TL, base_R=r["base_R"], f0_src=r["f0_src"], sec_f0=r["sec_f0"], f0_capped=r["f0_capped"],
               e_fixed=r["e_fixed"], e_gaps0=r["e_gaps0"], e_found_gen=r["e_found_gen"], R=r["R"], status=r["status"], gens=r["gens"], sec=r["sec"],
               judge_feasible=feas, judge_same=same, ok=int(ok))
    out.append(rec); print("G3", json.dumps(rec, ensure_ascii=False), flush=True)
    return int(not ok)


def main():
    """G1、G2、G3 依次跑，写 check_f0cap.jsonl。"""
    old, new = sys.argv[1], sys.argv[2]; out = []
    b1 = g1(old, new, out); b2 = g2(old, new, out); b3 = g3(new, out)
    summ = dict(G1=sum(r["part"] == "G1" for r in out), G1_bad=b1, G2_bad=b2, G3_bad=b3)
    with open(os.path.join(ROOT, "results", "T20", "check_f0cap.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
        f.write(json.dumps(summ) + "\n")
    print(json.dumps(summ), flush=True)


if __name__ == "__main__":
    main()
