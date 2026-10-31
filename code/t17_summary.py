"""T17 汇总：时间界的松弛一侧。
1. LP 下界随时域 H 的收敛：m ∈ {2, 3, 4, 6, 8, 10} 的 72 个训练实例，H ∈ {0.5, 1, 1.5, 2, 3, 4}·S_max，
   LP_fr(H)（不提前大修）与 LP_fr^Rw(H)（提前大修窗口 [0, min(2·S_max, H))）。
   来源：0.5、1、1.5、3 倍 = T17（results/T17/jobs/）；2 倍 = T2 的 fr__m2 / frR__m2（H = 2·S_max 时 LP_fr^Rw = LP_fr^R）；
   4 倍 = T2 的 fr__m4、T6 的 frRw__m4（正式运行失败时取 retry1）。
   LP 值取原始目标与证书对偶目标中较小者（与 T2 汇总相同）；对偶证书违反量 = max(路径、工具、等待、M)。
   H₀ = LP(H) ≥ 99.5%·LP(4S) 的最小 H（倍数）；另记 ⌈LP(H)⌉ = ⌈LP(4S)⌉ 的最小 H（整数下界意义上的收敛）。
2. 对偶价格的时间分布（H = 4·S_max）：覆盖对偶 π_{s,t}（results/T2_duals 的 fr__m4、results/T6_duals 的 frRw__m4），
   对偶目标 = Σ_{s,t} π_{s,t}·K_s；按每 0.5·S_max 一段统计这一段的贡献占总量的比例（净值；π 可以为负，另记负值部分的比例）。
   4S 的 LP 来自关 crossover 重跑（retry1）的实例不计入：内点对偶不是顶点，π 有成对抵消的大正负值。
3. 自然周期：已有全部可行解（T3 / T3b、T6 / T6b 的规则与 ES、T8、T11、T12 / T12b 的 RDD）的周期，
   可行 = 严格证明（verify_cert = 1 且 exact = 1）或 60 万格复核 0 缺岗；严格证明的解周期 = verify_period（格；1 格 = U 小时），
   只复核通过、判不出循环的记为"> 判环上限"（60 万格复核里没有回到出现过的状态）。按 m 分组。
输出 results/T17_horizon.jsonl（逐实例：LP 曲线、H₀、价格分布、周期列表）与 markdown 表（打印）。
用法：python3 code/t17_summary.py"""
import os, sys, json, glob, math
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, U, train_names, load_train  # noqa: E402

MULTS = (0.5, 1.0, 1.5, 2.0, 3.0, 4.0)
TVER = 600000                       # 复核时域（格）
YEAR = 365 * 24 / U                 # 一年的格数


def jl(p):
    """读 json；不存在返回 None。"""
    return json.load(open(p)) if os.path.exists(p) else None


def pick(base):
    """<base>.json 求解成功就用它，否则用 <base>__retry1.json；返回 (记录, 文件名主干) 或 (None, None)。"""
    for b in (base, base + "__retry1"):
        d = jl(b + ".json")
        if d is not None and d.get("obj") is not None: return d, os.path.basename(b)
    return None, None


def lp_value(d):
    """LP 值（原始目标与对偶目标的较小者）与证书违反量。"""
    v = min(d.get("obj_full", d["obj"]), d.get("dual_obj", d.get("obj_full", d["obj"])))
    viol = max(d.get("path_viol", math.nan), d.get("tool_viol", math.nan), d.get("wait_viol", math.nan), d.get("M_viol", math.nan))
    return v, viol


def lp_record(n, kind, mult):
    """(kind ∈ fr / frRw, 倍数) 的 LP 结果：(值, 违反量, 来源, 对偶文件) 或 None。"""
    T2J, T6J, T17J = (os.path.join(ROOT, "results", x, "jobs") for x in ("T2", "T6", "T17"))
    if mult in (2.0, 4.0):
        if kind == "fr": d, b = pick(os.path.join(T2J, f"{n}__fr__m{mult:g}")); src, dd = "T2", "T2_duals"
        elif mult == 2.0: d, b = pick(os.path.join(T2J, f"{n}__frR__m2")); src, dd = "T2 frR", "T2_duals"
        else: d, b = pick(os.path.join(T6J, f"{n}__frRw__m4")); src, dd = "T6", "T6_duals"
    else:
        d, b = pick(os.path.join(T17J, f"{n}__{kind}__m{mult:g}")); src, dd = "T17", "T17_duals"
    if d is None: return None
    v, viol = lp_value(d)
    return dict(LP=v, viol=viol, src=src, dual=os.path.join(ROOT, "results", dd, b + ".npz"), sec=d.get("wall"), H=d.get("H"), retry=int(b.endswith("__retry1")))


def price_shares(I, npz, S):
    """H = 4·S_max 的覆盖对偶：每 0.5·S_max 一段的净贡献占对偶目标的比例，以及负值部分占正值部分的比例。"""
    z = np.load(npz); pi = z["pi"]; H = pi.shape[1]
    K = np.array([I.sites[s]["K"] for s in range(pi.shape[0])], float)[:, None]
    c = pi * K; tot = c.sum(); seg = max(1, S // 2)
    shares = [float(c[:, a:a + seg].sum() / tot) for a in range(0, H, seg)]
    return dict(shares=shares[:8] if len(shares) > 8 else shares, tail=float(sum(shares[8:])) if len(shares) > 8 else 0.0,
                neg_ratio=float(-c[c < 0].sum() / c[c > 0].sum()), total=float(tot))


def feasible(d):
    """可行（严格证明或 60 万格复核 0 缺岗）时返回 (R, 周期或 None)；否则 None。"""
    if not d or d.get("R") is None or d["R"] < 0: return None
    proven = (d.get("verify_cert") == 1 and d.get("exact") == 1) or (d.get("proven") == 1 and d.get("verify_period", -1) > 0)
    if not (proven or d.get("verify_gaps") == 0): return None
    return d["R"], (d["verify_period"] if proven and d.get("verify_period", -1) > 0 else None)


def all_periods():
    """已有可行解的周期：{实例: [(来源, 周期或 None), …]}。"""
    src = [("T3", "T3/jobs/*.json"), ("T3b", "T3b/jobs/*.json"), ("T6", "T6/rlp/*.json"), ("T6", "T6/es/*.json"),
           ("T8", "T8/jobs/*.json"), ("T11", "T11/jobs/*.json"), ("T12", "T12/jobs/*.json"), ("T12b", "T12b/jobs/*.json")]
    out = defaultdict(list)
    for tag, g in src:
        for p in glob.glob(os.path.join(ROOT, "results", g)):
            d = json.load(open(p)); f = feasible(d)
            if f is None: continue
            out[os.path.basename(p).split("__")[0]].append((tag, f[1]))
    return out


def q(xs, a):
    """分位数（xs 非空）。"""
    return float(np.quantile(xs, a))


def main():
    names = sorted(train_names([2, 3, 4, 6, 8, 10]), key=lambda n: (int(n[1:3]), n))
    per = all_periods(); rows = []
    for n in names:
        I = load_train(n)
        from freshlb import S_max
        S = S_max(I)
        row = dict(inst=n, m=I.meta["m"], family=I.meta["family"], S=S)
        for kind in ("fr", "frRw"):
            recs = {m: lp_record(n, kind, m) for m in MULTS}
            row[f"{kind}_LP"] = {f"{m:g}": (r["LP"] if r else None) for m, r in recs.items()}
            row[f"{kind}_viol"] = max((r["viol"] for r in recs.values() if r), default=None)
            row[f"{kind}_sec"] = {f"{m:g}": (r["sec"] if r else None) for m, r in recs.items()}
            L4 = recs[4.0]["LP"] if recs[4.0] else None
            if L4:
                ratio = {f"{m:g}": (r["LP"] / L4 if r else None) for m, r in recs.items()}
                row[f"{kind}_ratio"] = ratio
                ok = [m for m in MULTS if recs[m] and recs[m]["LP"] >= 0.995 * L4 - 1e-9]
                row[f"{kind}_H0"] = min(ok) if ok else None
                c4 = math.ceil(L4 - 1e-6)
                okc = [m for m in MULTS if recs[m] and math.ceil(recs[m]["LP"] - 1e-6) >= c4]
                row[f"{kind}_H0_ceil"] = min(okc) if okc else None
                row[f"{kind}_complete"] = int(all(recs[m] for m in MULTS))
                row[f"{kind}_monotone"] = int(all(recs[a]["LP"] <= recs[b]["LP"] + 1e-6 for a, b in zip(MULTS, MULTS[1:]) if recs[a] and recs[b]))
                # 关 crossover 的重跑给的是内点对偶（不是顶点），π 有成对抵消的大正负值，时间分布没有意义，不进价格统计
                if recs[4.0]["retry"]: row[f"{kind}_price_skip"] = 1
                elif os.path.exists(recs[4.0]["dual"]): row[f"{kind}_price"] = price_shares(I, recs[4.0]["dual"], S)
        row["periods"] = per.get(n, [])
        rows.append(row)
    tmp = os.path.join(ROOT, "results", "T17_horizon.jsonl.tmp")
    with open(tmp, "w") as f:
        for r in rows: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(tmp, tmp[:-4])

    by = defaultdict(list)
    for r in rows: by[(r["m"], r["family"])].append(r)
    L = []
    for kind, title in (("fr", "LP_fr（不提前大修）"), ("frRw", "LP_fr^Rw（启动窗口内提前大修）")):
        L += [f"### {title}：LP(H) / LP(4S)，按 (m, 族) 平均\n",
              "| m | 族 | 完整 | " + " | ".join(f"{m:g}S" for m in MULTS) + " | H₀（99.5%）中位 / 最大 | H₀ ≤ 2S 的实例 | ⌈LP⌉ 收敛 H 中位 / 最大 |",
              "|---|---|---|" + "---|" * len(MULTS) + "---|---|---|"]
        for (m, fam), xs in sorted(by.items()):
            ok = [x for x in xs if x.get(f"{kind}_complete")]
            cols = []
            for mm in MULTS:
                v = [x[f"{kind}_ratio"][f"{mm:g}"] for x in ok]
                cols.append(f"{100 * sum(v) / len(v):.1f}%" if v else "—")
            h0 = [x[f"{kind}_H0"] for x in ok if x.get(f"{kind}_H0") is not None]
            hc = [x[f"{kind}_H0_ceil"] for x in ok if x.get(f"{kind}_H0_ceil") is not None]
            tail = (f" | {np.median(h0):g} / {max(h0):g} | {sum(h <= 2 for h in h0)}/{len(h0)} | {np.median(hc):g} / {max(hc):g} |" if h0 and hc else " | — | — | — |")
            L.append(f"| {m} | {fam} | {len(ok)}/{len(xs)} | " + " | ".join(cols) + tail)
        vio = [x[f"{kind}_viol"] for x in rows if x.get(f"{kind}_viol") is not None]
        nm = sum(1 for x in rows if x.get(f"{kind}_monotone") == 0)
        L += ["", f"对偶证书最大违反量 {max(vio):.1e}（{len(vio)} 个实例的全部 H）；LP(H) 不随 H 单调的实例 {nm} 个。", ""]
    for kind, title in (("fr", "LP_fr(4S)"), ("frRw", "LP_fr^Rw(4S)")):
        L += [f"### 对偶价格的时间分布（{title}）：每 0.5·S_max 一段占对偶目标的比例，按 (m, 族) 平均\n",
              "| m | 族 | " + " | ".join(f"[{a / 2:g}, {(a + 1) / 2:g})S" for a in range(8)) + " | 负值 / 正值 |", "|---|---|" + "---|" * 9]
        for (m, fam), xs in sorted(by.items()):
            ps = [x[f"{kind}_price"] for x in xs if x.get(f"{kind}_price")]
            if not ps: continue
            cols = [f"{100 * np.mean([p['shares'][a] for p in ps if len(p['shares']) > a]):.1f}%" for a in range(8)]
            L.append(f"| {m} | {fam} | " + " | ".join(cols) + f" | {100 * np.mean([p['neg_ratio'] for p in ps]):.1f}% |")
        sk = [x["inst"] for x in rows if x.get(f"{kind}_price_skip")]
        L += ["", f"净值占比（π 可为负，各段相加为 100%）；最后一列 = 负贡献总量 / 正贡献总量。4S 的 LP 来自关 crossover 重跑的实例不计入（内点对偶）：{', '.join(sk) or '无'}。", ""]
    # 自然周期（按 m；含 m = 15、20 的 T12 解）
    allp = all_periods(); bym = defaultdict(list)
    for n, ps in allp.items(): bym[int(n[1:3])].extend(ps)
    L += ["### 自然周期（已有全部可行解，按 m）\n",
          f"| m | 可行解 | 判出周期 | 周期中位数（格 / 年） | 25%–75%（年） | 最长（年） | 判不出（> {TVER} 格 ≈ {TVER / YEAR:.0f} 年）|", "|---|---|---|---|---|---|---|"]
    for m in sorted(bym):
        ps = bym[m]; pp = [p for _, p in ps if p is not None]
        if pp:
            L.append(f"| {m} | {len(ps)} | {len(pp)} | {q(pp, .5):.0f} / {q(pp, .5) / YEAR:.1f} | {q(pp, .25) / YEAR:.1f}–{q(pp, .75) / YEAR:.1f} | {max(pp) / YEAR:.1f} | {len(ps) - len(pp)} |")
        else:
            L.append(f"| {m} | {len(ps)} | 0 | — | — | — | {len(ps)} |")
    print("\n".join(L))


if __name__ == "__main__":
    main()
