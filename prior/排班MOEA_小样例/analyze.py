"""分析 runs.pkl -> summary.md
HV：每个实例把所有运行的可行非支配解合在一起确定归一化范围（理想点=最小值，参考点=最大值×1.1 后归一化到 1.1）。
同时报告：可行解出现的代数、最小资源数、"资源数 ≤ 阈值"时的最好均衡度。
"""
import pickle, sys
from math import erf, sqrt
import numpy as np

fn = sys.argv[1] if len(sys.argv) > 1 else "runs.pkl"
res = sum((pickle.load(open(f, "rb")) for f in fn.split("+")), [])   # 支持 a.pkl+b.pkl 合并
fn = fn.split("+")[0]
insts = list(dict.fromkeys(r["inst"] for r in res))
vars_ = list(dict.fromkeys(r["var"] for r in res))
base = vars_[0]


def hv2(P, ref):
    P = P[(P < ref).all(1)]
    if not len(P):
        return 0.0
    P = P[np.argsort(P[:, 0])]
    h, prev = 0.0, ref[1]
    for f1, f2 in P:
        if f2 < prev:
            h += (ref[0] - f1) * (prev - f2); prev = f2
    return h


def wilcoxon(a, b):
    d = np.asarray(a) - np.asarray(b); d = d[d != 0]; n = len(d)
    if n == 0:
        return 1.0
    ad = np.abs(d); rk = np.argsort(np.argsort(ad)) + 1.0
    for v in np.unique(ad):
        rk[ad == v] = rk[ad == v].mean()
    W = rk[d > 0].sum(); mu = n * (n + 1) / 4; sd = np.sqrt(n * (n + 1) * (2 * n + 1) / 24)
    z = (abs(W - mu) - .5) / sd
    return 2 * (1 - .5 * (1 + erf(z / sqrt(2))))


out = ["# 排班 MOEA 对比结果（自动生成）\n",
       "HV 越大越好（每个实例统一归一化）；p 为与基线的配对 Wilcoxon（同种子配对）。"
       "“首个可行代”只在随机解大多不可行的实例上有意义。\n"]
summary_rank = {v: [] for v in vars_}
for I in insts:
    rs = [r for r in res if r["inst"] == I]
    allP = np.vstack([r["ndF"] for r in rs if len(r["ndF"])])
    lo, hi = allP.min(0), allP.max(0)
    ref = np.full(2, 1.1)
    for r in rs:
        r["hv"] = hv2((r["ndF"] - lo) / np.maximum(hi - lo, 1e-9), ref) if len(r["ndF"]) else 0.0
        r["minf1"] = r["ndF"][:, 0].min() if len(r["ndF"]) else np.nan
    seeds = sorted(set(r["seed"] for r in rs))
    get = lambda v, key: np.array([next(r[key] for r in rs if r["var"] == v and r["seed"] == s) for s in seeds], float)
    out.append(f"## {I}\n")
    out.append("| 变体 | HV 均值 | HV 标准差 | p(vs 基线) | 最小资源数（均值） | 首个可行代（均值） | 单次用时 s |")
    out.append("|---|---|---|---|---|---|---|")
    hvs = {v: get(v, "hv") for v in vars_}
    order = sorted(vars_, key=lambda v: -hvs[v].mean())
    for k, v in enumerate(order):
        summary_rank[v].append(k + 1)
    for v in vars_:
        ff = [r["first_feas"] for r in rs if r["var"] == v]
        ffm = np.mean([f for f in ff if f is not None]) if any(f is not None for f in ff) else np.nan
        nofeas = sum(f is None for f in ff)
        p = "—" if v == base else f"{wilcoxon(hvs[v], hvs[base]):.3f}"
        out.append(f"| {v} | {hvs[v].mean():.4f} | {hvs[v].std():.4f} | {p} | {np.nanmean(get(v, 'minf1')):.1f} | "
                   f"{ffm:.1f}{'（'+str(nofeas)+' 次无可行解）' if nofeas else ''} | "
                   f"{np.mean([r['secs'] for r in rs if r['var'] == v]):.1f} |")
    # 纯规则参照：全 0 基因 = 每班都复用（原解码按编号取第一个；排序解码取最闲的）
    try:
        import roster_sim as RS
        rules = []
        for d, nm in [("orig", "按编号复用"), ("rank", "最闲优先")]:
            inst = RS.INSTANCES[I](); inst.decoder = d
            o = RS.evaluate(np.zeros(inst.n_var), inst)
            rules.append((nm, o["f"] if o["feasible"] else None))
        txt = []
        for nm, f in rules:
            if f is None:
                txt.append(f"{nm}：不可行"); continue
            cover = {v: np.mean([((r["ndF"] <= f).all(1)).any() for r in rs if r["var"] == v]) for v in vars_}
            best = max(cover, key=cover.get)
            txt.append(f"{nm}：(资源数 {f[0]:.0f}, 均衡度 {f[1]:.1f})，被 MOEA 前沿弱支配的运行比例最高为 {cover[best]:.0%}（{best}）")
        out.append("\n纯规则参照（每班都复用）：" + "；".join(txt))
    except Exception as e:
        out.append(f"\n（规则参照计算失败：{e}）")
    out.append("")
out.append("## 跨实例 HV 平均排名（越小越好）\n")
out.append("| 变体 | 平均排名 |"); out.append("|---|---|")
for v in sorted(vars_, key=lambda v: np.mean(summary_rank[v])):
    out.append(f"| {v} | {np.mean(summary_rank[v]):.2f} |")
open(fn.replace(".pkl", "_summary.md"), "w", encoding="utf-8").write("\n".join(out) + "\n")
print("\n".join(out))
