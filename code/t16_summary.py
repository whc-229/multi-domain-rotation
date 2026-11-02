"""T16 汇总：对偶引导偏离位置（β = 0 均匀、0.5、1）。
来源：results/T16/jobs/<实例>__b<β>__s<种子>.json（t16_job.py），m ∈ {4, 6, 8, 10} 的 48 个训练实例 × 3 个 β × 种子 0–2。
可行 = 严格证明或 60 万格复核 0 缺岗；只在可行解上比资源数。
1. 每组 (m, 族) 的平均资源数（实例 × 种子平均），可行 / 严格证明个数；
2. 对 β = 0 的配对比较：同一 (实例, 种子) 的资源数差（少 / 同 / 多、平均差），以及同一实例 3 个种子里最好值的差；
3. 偏离落在哪里：把最终解码出的每个激活偏离对应到它所在的窗口决策 (站 s, 时刻 t, 值岗时长 d)，
   用 LP_fr^Rw(4S) 的覆盖对偶算该决策的权重 w = Σ_{u∈[t, t+d)} max(π_{s,u}, 0)（与 t16_job 的定义相同），
   再与同一次运行里全部窗口决策的权重比：落在 w > 0 的决策上的占比（对照：全部决策中 w > 0 的占比）、
   落在权重前 20% 的占比（对照 20%）、平均分位（中位秩，对照 0.5）；
4. 时间分布：[0, 2·S_max) 按 0.25·S_max 分 8 段，对偶价格总量的占比（各实例归一化后平均）与激活偏离的占比（按 β）。
偏离位置用 Python 参考实现（rost5_ref，带同一个 cdf）解码，激活个数与 rost5 输出的 n_active 核对。
输出 results/T16_dual_guided.jsonl（逐作业）与 markdown 表（打印）。
用法：python3 code/t16_summary.py"""
import os, sys, json
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names, load_train  # noqa: E402
from freshlb import S_max  # noqa: E402
from t6_rerun import frrw_npz  # noqa: E402
from t12_job import RULES  # noqa: E402
import rost5_ref as ref  # noqa: E402

BETAS = (0.0, 0.5, 1.0)
NB = 8           # 时间分布的段数：[0, 2S) 每段 0.25·S_max


def feas(d):
    """可行时返回资源数，否则 None。"""
    if not d or d.get("R") is None or d["R"] < 0: return None
    return d["R"] if ((d.get("verify_cert") == 1 and d.get("exact") == 1) or d.get("verify_gaps") == 0) else None


def dec_weight(pi, s, t, dur):
    """窗口决策 (站 s, 时刻 t, 值岗时长 dur) 的对偶权重 Σ_{u∈[t, t+dur)} max(π_{s,u}, 0)。"""
    seg = pi[s, t:min(t + dur, pi.shape[1])]
    return float(np.clip(seg, 0, None).sum())


def mid_rank(ws, x):
    """x 在 ws 中的中位秩分位：(严格更小的个数 + 0.5·相等的个数) / 总数。"""
    ws = np.asarray(ws)
    return float(((ws < x).sum() + 0.5 * (ws == x).sum()) / len(ws))


def deviations(I, n, d, S, pi):
    """把最终解解码，返回 (激活偏离 [(t, s, dur, w, 分位, 是否前 20%)], 全部窗口决策的权重, 参考实现的激活个数)。"""
    D = ref.decode_geno(d["genes"], d["W"], d["ndev"], d["qon"], True, cdf=d.get("cdf") or None)
    R = ref.parse_rlp(I, f"{RULES}/{n}__fresh_rlp_w.rlp"); PF = ref.parse_planf(I, f"{RULES}/{n}__fresh_rlp_w.planf")
    tr = {}
    ref.run(I, d["Nk"], d["Mj"], d["W"], False, rlp=R, planf=PF, rdd=D, plan_until=2 * S, trace=tr)
    decs = tr.get("decs", [])
    ws = [dec_weight(pi, s, t, dur) for s, t, dur in decs]
    q80 = np.quantile(ws, 0.8) if ws else 0.0
    out = []
    for c, u in enumerate(D["dev"]):
        if u < 0 or c >= len(decs): continue
        s, t, dur = decs[c]; w = ws[c]
        out.append((t, s, dur, w, mid_rank(ws, w), int(w > q80)))
    return out, ws, sum(u >= 0 for u in D["dev"])


def main():
    rows = []; dual_hist = defaultdict(list)
    for n in sorted(train_names([4, 6, 8, 10]), key=lambda n: (int(n[1:3]), n)):
        I = load_train(n); S = S_max(I); m = I.meta["m"]
        pi = np.load(frrw_npz(n))["pi"]
        # 对偶价格总量在 [0, 2S) 上按 0.25·S 一段的占比（只取正价格）
        pp = np.clip(pi[:, :2 * S], 0, None).sum(axis=0)
        h = np.array([pp[int(b * S / 4):int((b + 1) * S / 4)].sum() for b in range(NB)])
        dual_hist[m].append(h / h.sum() if h.sum() > 0 else h)
        for beta in BETAS:
            for seed in (0, 1, 2):
                p = os.path.join(ROOT, "results", "T16", "jobs", f"{n}__b{beta:g}__s{seed}.json")
                if not os.path.exists(p): continue
                d = json.load(open(p)); R = feas(d)
                dev, ws, nact = deviations(I, n, d, S, pi) if (R is not None and d.get("genes")) else ([], [], None)
                rows.append(dict(inst=n, m=m, family=I.meta["family"], S=S, beta=beta, seed=seed, R=R, base_R=d.get("base_R"),
                                 proven=int(d.get("verify_cert") == 1 and d.get("exact") == 1), window_ok=d.get("window_ok"),
                                 n_active=d.get("n_active"), n_active_ref=nact, ndev=d.get("ndev"), K=d.get("K"), wall=d.get("wall"),
                                 check_judge=d.get("check_judge"), check_dump=(d.get("check_dump_same"), d.get("check_dump_n")),
                                 dec_pos_frac=float(np.mean([w > 0 for w in ws])) if ws else None,
                                 dev=[[t, s, dur, round(w, 6), round(pc, 4), top] for t, s, dur, w, pc, top in dev]))
    tmp = os.path.join(ROOT, "results", "T16_dual_guided.jsonl.tmp")
    with open(tmp, "w") as f:
        for r in rows: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(tmp, tmp[:-4])

    by = defaultdict(list)
    for r in rows: by[(r["m"], r["family"])].append(r)
    L = ["### 1. 平均资源数（可行解；实例 × 种子平均；括号内为可行 / 作业数，严格证明数）\n",
         "| m | 族 | 默认规则 | " + " | ".join(f"β = {b:g}" for b in BETAS) + " |", "|---|---|---|" + "---|" * len(BETAS)]
    for (m, fam), xs in sorted(by.items()):
        base = np.mean([x["base_R"] for x in xs if x["beta"] == 0 and x["base_R"] is not None and x["base_R"] >= 0])
        cols = []
        for b in BETAS:
            ys = [x for x in xs if x["beta"] == b]; ok = [x["R"] for x in ys if x["R"] is not None]
            cols.append(f"{np.mean(ok):.2f}（{len(ok)}/{len(ys)}，{sum(x['proven'] for x in ys)}）" if ok else f"—（0/{len(ys)}）")
        L.append(f"| {m} | {fam} | {base:.2f} | " + " | ".join(cols) + " |")

    def paired(xs, b):
        """同一 (实例, 种子) 上 β = b 与 β = 0 的差（b − 0）。"""
        d0 = {(x["inst"], x["seed"]): x["R"] for x in xs if x["beta"] == 0}
        return [x["R"] - d0[(x["inst"], x["seed"])] for x in xs if x["beta"] == b and x["R"] is not None and d0.get((x["inst"], x["seed"])) is not None]

    def best_diff(xs, b):
        """同一实例 3 个种子里最好值的差（b − 0）。"""
        bb = defaultdict(dict)
        for x in xs:
            if x["R"] is not None: bb[x["inst"]][x["beta"]] = min(bb[x["inst"]].get(x["beta"], 10 ** 9), x["R"])
        return [v[b] - v[0.0] for v in bb.values() if b in v and 0.0 in v]

    def fmt(ds):
        """少 / 同 / 多（对数）与平均差。"""
        return f"{sum(x < 0 for x in ds)} / {sum(x == 0 for x in ds)} / {sum(x > 0 for x in ds)}，{np.mean(ds):+.2f}" if ds else "—"

    L += ["", "### 2. 与 β = 0 配对比较（β 减 β = 0；少 / 同 / 多，平均差）\n",
          "| m | 族 | β = 0.5（同种子） | β = 1（同种子） | β = 0.5（3 种子最好） | β = 1（3 种子最好） |", "|---|---|---|---|---|---|"]
    for (m, fam), xs in sorted(by.items()):
        L.append(f"| {m} | {fam} | {fmt(paired(xs, 0.5))} | {fmt(paired(xs, 1.0))} | {fmt(best_diff(xs, 0.5))} | {fmt(best_diff(xs, 1.0))} |")
    for b in (0.5, 1.0):
        L.append(f"\n合计 β = {b:g}：同种子 {fmt(paired(rows, b))}；3 种子最好 {fmt(best_diff(rows, b))}。")

    L += ["", "### 3. 偏离落在对偶价格高的决策上吗（最终解里激活的偏离；对照 = 同一次运行全部窗口决策）\n",
          "| m | β | 作业 | 平均激活偏离数 | 落在 w > 0 的决策（对照） | 落在权重前 20%（对照 20%） | 平均分位（对照 0.5） |", "|---|---|---|---|---|---|---|"]
    for m in (4, 6, 8, 10):
        for b in BETAS:
            xs = [x for x in rows if x["m"] == m and x["beta"] == b and x["R"] is not None]
            dv = [e for x in xs for e in x["dev"]]
            if not dv: continue
            L.append(f"| {m} | {b:g} | {len(xs)} | {np.mean([len(x['dev']) for x in xs]):.1f} | {100 * np.mean([e[3] > 0 for e in dv]):.0f}%"
                     f"（{100 * np.mean([x['dec_pos_frac'] for x in xs]):.0f}%） | {100 * np.mean([e[5] for e in dv]):.0f}% | {np.mean([e[4] for e in dv]):.2f} |")

    L += ["", "### 4. 时间分布：[0, 2·S_max) 每 0.25·S_max 一段的占比（对偶价格总量 vs 激活偏离）\n",
          "| m | 来源 | " + " | ".join(f"[{b / 4:g}, {(b + 1) / 4:g})S" for b in range(NB)) + " |", "|---|---|" + "---|" * NB]
    for m in (4, 6, 8, 10):
        dh = np.mean(dual_hist[m], axis=0)
        L.append(f"| {m} | 对偶价格 | " + " | ".join(f"{100 * v:.0f}%" for v in dh) + " |")
        for b in BETAS:
            ts = [e[0] / x["S"] for x in rows if x["m"] == m and x["beta"] == b for e in x["dev"]]
            if not ts: continue
            h = np.histogram(ts, bins=np.linspace(0, 2, NB + 1))[0]
            L.append(f"| {m} | 偏离 β = {b:g}（{len(ts)}） | " + " | ".join(f"{100 * c / len(ts):.0f}%" for c in h) + " |")

    cj = [r["check_judge"] for r in rows if r["check_judge"] is not None]
    na = [(r["n_active"], r["n_active_ref"]) for r in rows if r["n_active_ref"] is not None]
    L += ["", f"核对：作业 {len(rows)} / 432；最终解重新判定逐字段相同 {sum(cj)} / {len(cj)}；抽样缺岗与 eval 相同 "
          f"{sum(r['check_dump'][0] or 0 for r in rows)} / {sum(r['check_dump'][1] or 0 for r in rows)}；"
          f"参考实现解码的激活偏离数与 rost5 相同 {sum(a == b for a, b in na)} / {len(na)}；"
          f"提前大修都在 [0, 2S) 内 {sum(r['window_ok'] == 1 for r in rows if r['R'] is not None)} / {sum(r['R'] is not None for r in rows)}。"]
    print("\n".join(L))


if __name__ == "__main__":
    main()
