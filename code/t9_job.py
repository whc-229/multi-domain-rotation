"""T9 单个作业：整数机队松弛 IFR / IFR^Rw。结果写进 results/T9/jobs/<名>.json（原子写入，已存在就跳过）。
模型：freshlb.fresh_lb(integer=True)，即全新启动时间展开 LP 里只把各型台数 N_k、各种工具数 M_j 设为整数，流量仍连续。H = 2·S_max。
- fr：不提前大修（IFR）；
- frRw：提前大修弧只在值岗开始 D < W = 2·S_max 时存在（IFR^Rw）。H = 2·S_max 时窗口覆盖整个时域，
  所以 IFR^Rw(2S) 就是 freshlb 的 reset=True 整数机队模型（越过 H 的值岗不进任何约束，有没有都一样）。
有效性：任何可行排班截到 H 后满足 LP 的全部约束，且台数本来就是整数，所以最优值（以及 MIP 对偶界）是下界。
求解：建模用 freshlb 原代码，只把 scipy 的 milp 换成 highspy（code/t4_job._milp_highs：单线程、mip_rel_gap = 0、限时）。
下界 = ⌈对偶界 − 1e-6⌉。
用法：python3 code/t9_job.py <实例名> <fr|frRw> [限时秒，默认 3600]"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train  # noqa: E402
import freshlb  # noqa: E402
from t4_job import _milp_highs  # noqa: E402

OUT = os.environ.get("T9_OUT", os.path.join(ROOT, "results", "T9", "jobs"))


def ifr(I, reset, tl, mult=2):
    """整数机队松弛：返回 dict(H, obj, bound, LB, status, N, M, sec)。"""
    H = int(mult * freshlb.S_max(I))
    old = freshlb.milp; freshlb.milp = _milp_highs(tl)
    try:
        t0 = time.time()
        r = freshlb.fresh_lb(I, H, tools=True, reset=reset, integer=True, time_limit=tl)
        sec = time.time() - t0
    finally:
        freshlb.milp = old
    b = r.get("bound")
    return dict(H=H, obj=r.get("obj"), bound=b, LB=int(np.ceil(b - 1e-6)) if b is not None and np.isfinite(b) else None,
                status=r.get("status"), N=r.get("N"), M=r.get("M"), n=r.get("n"), sec=round(sec, 1))


def main():
    inst, kind = sys.argv[1], sys.argv[2]
    tl = float(sys.argv[3]) if len(sys.argv) > 3 else 3600
    os.makedirs(OUT, exist_ok=True)
    name = f"{inst}__ifr_{kind}"
    path = os.path.join(OUT, name + ".json")
    if os.path.exists(path): print("已完成，跳过", name); return
    I = load_train(inst)
    rec = dict(inst=inst, m=I.meta["m"], family=I.meta["family"], kind=kind, tl=tl)
    rec.update(ifr(I, kind == "frRw", tl))
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps(rec, ensure_ascii=False))


if __name__ == "__main__":
    main()
