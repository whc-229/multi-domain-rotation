"""T18 汇总：时间界的编码一侧（RDD 偏离窗口 W 的扫描）。
来源：W = 0 → 默认规则本身（m = 4、6：T12b 的严格窗口 base；m = 8、10：T12 的 base）；
      W = 2·S_max → T16 的 β = 0（设置与 T18 相同，限时相同）；W = 0.5、1、3、4 倍 → T18（results/T18/jobs/）。
可行 = 严格证明或 60 万格复核 0 缺岗；只在可行解上比资源数。
1. 每组 (m, 族) 的平均资源数随 W 的变化（实例 × 种子平均）；
2. W₀（按实例）：每个 W 取 3 个种子里最少的资源数 R(W)，W₀ = R(W) 达到全部 W 中最小值的最小 W；另给"差 1 个以内"的版本；
3. W₀ 与 T17 的 H₀（LP_fr^Rw 达到 4S 值 99.5% 的最小 H）逐实例对照；
4. 偏离落点：最终解里激活的偏离的时刻 t，按每 0.5·S_max 一段统计占比（W = 2S 的 T16 作业没有记，这里用 t18_job.dev_positions 补算）。
输出 results/T18_window.jsonl（逐作业）与 markdown 表（打印）。
用法：python3 code/t18_summary.py"""
import os, sys, json
from collections import defaultdict, Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names, load_train  # noqa: E402
from freshlb import S_max  # noqa: E402
from t18_job import dev_positions  # noqa: E402

MULTS = (0.0, 0.5, 1.0, 2.0, 3.0, 4.0)


def jl(p):
    """读 json；不存在返回 None。"""
    return json.load(open(p)) if os.path.exists(p) else None


def feas(d):
    """可行时返回资源数，否则 None。"""
    if not d or d.get("R") is None or d["R"] < 0: return None
    return d["R"] if ((d.get("verify_cert") == 1 and d.get("exact") == 1) or d.get("verify_gaps") == 0) else None


def job(n, mult, seed):
    """(实例, 倍数, 种子) 的结果记录（W = 0 不分种子）。"""
    if mult == 0.0:
        return jl(os.path.join(ROOT, "results", "T12b" if int(n[1:3]) <= 6 else "T12", "jobs", f"{n}__base.json"))
    if mult == 2.0: return jl(os.path.join(ROOT, "results", "T16", "jobs", f"{n}__b0__s{seed}.json"))
    return jl(os.path.join(ROOT, "results", "T18", "jobs", f"{n}__w{mult:g}__s{seed}.json"))


def main():
    H0 = {}
    if os.path.exists(os.path.join(ROOT, "results", "T17_horizon.jsonl")):
        H0 = {r["inst"]: r.get("frRw_H0") for r in map(json.loads, open(os.path.join(ROOT, "results", "T17_horizon.jsonl")))}
    rows = []; inst_rows = []
    for n in sorted(train_names([4, 6, 8, 10]), key=lambda n: (int(n[1:3]), n)):
        I = load_train(n); S = S_max(I)
        best = {}
        for mult in MULTS:
            Rs = []
            for seed in ((None,) if mult == 0.0 else (0, 1, 2)):
                d = job(n, mult, seed)
                if d is None: continue
                R = feas(d)
                pos = d.get("dev_pos")
                if mult == 2.0 and pos is None and R is not None and d.get("genes"):
                    pos, _ = dev_positions(I, n, d, 2 * S)
                rows.append(dict(inst=n, m=I.meta["m"], family=I.meta["family"], S=S, mult=mult, seed=seed, R=R, status=d.get("status"),
                                 proven=int(d.get("verify_cert") == 1 and d.get("exact") == 1), n_active=d.get("n_active"), ndev=d.get("ndev"), K=d.get("K"),
                                 window_ok=int(d.get("last_reset", -1) < 2 * S), wall=d.get("wall"), check_judge=d.get("check_judge"),
                                 check_dump=(d.get("check_dump_same"), d.get("check_dump_n")), check_ndec=d.get("check_ndec"),
                                 dev_t=[p[2] for p in pos] if pos else []))
                if R is not None: Rs.append(R)
            best[mult] = min(Rs) if Rs else None
        vals = {m: r for m, r in best.items() if r is not None}
        w0 = w0_1 = None
        if len(vals) == len(MULTS):
            lo = min(vals.values())
            w0 = min(m for m in MULTS if vals[m] == lo); w0_1 = min(m for m in MULTS if vals[m] <= lo + 1)
        inst_rows.append(dict(inst=n, m=I.meta["m"], family=I.meta["family"], best={f"{m:g}": best[m] for m in MULTS}, W0=w0, W0_1=w0_1, H0=H0.get(n)))
    tmp = os.path.join(ROOT, "results", "T18_window.jsonl.tmp")
    with open(tmp, "w") as f:
        for r in rows: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(tmp, tmp[:-4])

    by = defaultdict(list)
    for r in rows: by[(r["m"], r["family"])].append(r)
    L = ["### 平均资源数随偏离窗口 W 的变化（可行解；实例 × 种子平均；括号内为可行 / 作业数）\n",
         "| m | 族 | " + " | ".join(f"W = {m:g}S" for m in MULTS) + " |", "|---|---|" + "---|" * len(MULTS)]
    for (m, fam), xs in sorted(by.items()):
        cols = []
        for mm in MULTS:
            ys = [x for x in xs if x["mult"] == mm]; ok = [x["R"] for x in ys if x["R"] is not None]
            cols.append(f"{np.mean(ok):.1f}（{len(ok)}/{len(ys)}）" if ok else f"—（0/{len(ys)}）")
        L.append(f"| {m} | {fam} | " + " | ".join(cols) + " |")
    byi = defaultdict(list)
    for r in inst_rows: byi[(r["m"], r["family"])].append(r)
    L += ["", "### W₀（按实例；每个 W 取 3 个种子里的最少资源数）\n",
          "| m | 族 | 完整实例 | W₀ 的分布（最好值） | W₀ 的分布（差 1 个以内） | 比 W = 0 少（个，最好种子，按 W） |", "|---|---|---|---|---|---|"]
    for (m, fam), xs in sorted(byi.items()):
        ok = [x for x in xs if x["W0"] is not None]
        d0 = Counter(f"{x['W0']:g}S" for x in ok); d1 = Counter(f"{x['W0_1']:g}S" for x in ok)
        imp = []
        for mm in MULTS[1:]:
            v = [x["best"]["0"] - x["best"][f"{mm:g}"] for x in ok]
            imp.append(f"{mm:g}S: {np.mean(v):.1f}" if v else "—")
        L.append(f"| {m} | {fam} | {len(ok)}/{len(xs)} | " + ", ".join(f"{k} ×{v}" for k, v in sorted(d0.items(), key=lambda kv: float(kv[0][:-1]))) + " | "
                 + ", ".join(f"{k} ×{v}" for k, v in sorted(d1.items(), key=lambda kv: float(kv[0][:-1]))) + " | " + "；".join(imp) + " |")
    pairs = [(x["W0"], x["H0"]) for x in inst_rows if x["W0"] is not None and x["H0"] is not None]
    if pairs:
        L += ["", "### W₀ 与 T17 的 H₀（LP_fr^Rw）逐实例对照：行 = W₀，列 = H₀，格内为实例数\n"]
        hs = sorted({h for _, h in pairs}); ws = sorted({w for w, _ in pairs})
        L += ["| W₀ \\ H₀ | " + " | ".join(f"{h:g}S" for h in hs) + " |", "|---|" + "---|" * len(hs)]
        for w in ws: L.append(f"| {w:g}S | " + " | ".join(str(sum(1 for a, b in pairs if a == w and b == h)) for h in hs) + " |")
        L.append(f"\nW₀ ≤ H₀ 的实例 {sum(a <= b for a, b in pairs)} / {len(pairs)}；W₀（差 1 个以内）≤ H₀ 的实例 "
                 f"{sum(x['W0_1'] <= x['H0'] for x in inst_rows if x['W0'] is not None and x['H0'] is not None)} / {len(pairs)}。")
    L += ["", "### 偏离落点：激活的偏离按时刻每 0.5·S_max 一段的占比（全部作业合并，按 m、W）\n",
          "| m | W | 偏离数 | " + " | ".join(f"[{a / 2:g}, {(a + 1) / 2:g})S" for a in range(8)) + " |", "|---|---|---|" + "---|" * 8]
    for m in (4, 6, 8, 10):
        for mm in MULTS[2:]:
            ts = [t / x["S"] for x in rows if x["m"] == m and x["mult"] == mm for t in x["dev_t"]]
            if not ts: continue
            h = np.histogram(ts, bins=np.arange(0, 4.01, 0.5))[0]
            L.append(f"| {m} | {mm:g}S | {len(ts)} | " + " | ".join(f"{100 * c / len(ts):.0f}%" if a < 2 * mm else "" for a, c in enumerate(h)) + " |")
    t18 = [r for r in rows if r["mult"] not in (0.0, 2.0)]
    cj = [r["check_judge"] for r in t18 if r["check_judge"] is not None]
    L += ["", f"核对（T18 作业）：最终解重新判定逐字段相同 {sum(cj)} / {len(cj)}；抽样缺岗与 eval 相同 {sum(r['check_dump'][0] or 0 for r in t18)} / {sum(r['check_dump'][1] or 0 for r in t18)}；"
          f"参考实现数出的窗口决策数与 eval 相同 {sum(r['check_ndec'] == 1 for r in t18)} / {sum(r['check_ndec'] is not None for r in t18)}；"
          f"提前大修都在 [0, 2S) 内 {sum(r['window_ok'] for r in rows if r['R'] is not None)} / {sum(r['R'] is not None for r in rows)}。"]
    print("\n".join(L))


if __name__ == "__main__":
    main()
