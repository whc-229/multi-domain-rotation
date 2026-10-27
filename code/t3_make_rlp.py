"""T3：由 T2 的 LP 解生成 R_LP 规则文件（rost5 rlp=、Python 参考实现共用）和 LP 提前大修计划文件（planf）。

规则文件全部用整数，C++ 与 Python 做同样的整数运算，结果逐位一致：
- H：计划时域（全新启动 = LP 的 H；稳态 = 0）；steady：是否稳态初值。
- NLPq：LP 台数 × 1e6；N0 / M0：LP 台数、工具数向上取整（定台数起点）。
- Q_s、p_{k,s}：t ≥ H 后的计划增量。p/Q_s 是速率 LP 班率 y_{k,s} 的小分母有理近似，要求 Σ_k p_{k,s} d_k = K_s Q_s（与覆盖一致，lag 不漂移）。
  取 Q_s ∈ [1, 30] 中最大误差 ≤ 0.01 班/格的最小者；都达不到就取误差最小的。分母小，lag 的周期短，证书容易判出。
- CQ_{k,s} = ⌊40·Q_s / d_k⌋（lag 截断在 ±40 格对应的班数）；CP_k = ⌊40·L_k / d_k⌋（相位 lag 截断，单位 1/L_k）。
- F[k][s][t] = round(1e6 × Σ_{D ≤ t} LP 值岗流_{k,s}(D))，FP[k][i][t] 同理按相位，t ∈ [0, H)。
- W[k][i]（稳态初值）：相位 i 的就绪窗口 = 进入相位 i 的上一段（相位 i−1）在各允许站的平均步长，取整（≥ 1）。
planf：与 prior 的 make_planf.py 相同的格式与口径（#N LP 台数；k i D 流量，只保留 D < H/2 的提前大修流）。
诊断变体 *_c：缩放分母用 ⌈N_LP⌉（NLPq = 1e6·⌈N_LP⌉，planf 的 #N 同样取整），其余相同。
用法：python3 code/t3_make_rlp.py <实例名>    → results/T3/rlp/<实例名>__{fresh_rlp, fresh_noreset, steady, fresh_rlp_c, fresh_noreset_c}.rlp 与 .planf"""
import os, sys, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train  # noqa: E402

T2J = os.path.join(ROOT, "results", "T2", "jobs")
T2D = os.path.join(ROOT, "results", "T2_duals")
OUT = os.path.join(ROOT, "results", "T3", "rlp")
QMAX, QTOL, CLIP_SLOTS = 30, 0.01, 40


def rationalize(I, s, y):
    """站 s：返回 (Q, {k: p}, 最大误差)。y: {k: 班率}。"""
    ks = [k for k in range(len(I.T)) if s in I.elig[k]]
    K = I.sites[s]["K"]; best = None
    for Q in range(1, QMAX + 1):
        tgt = K * Q
        # 枚举除最后一型外的 p，最后一型由等式确定
        def rec(a, rem, cur):
            nonlocal best
            k = ks[a]; d = I.d[k]
            if a == len(ks) - 1:
                if rem % d: return
                pp = dict(cur); pp[k] = rem // d
                err = max(abs(pp[kk] / Q - y.get(kk, 0.0)) for kk in ks)
                key = (err > QTOL, Q if err <= QTOL else 0, err)
                if best is None or key < best[0]: best = (key, Q, pp, err)
                return
            for x in range(rem // d + 1):
                cur[k] = x; rec(a + 1, rem - x * d, cur)
            cur.pop(k, None)
        rec(0, tgt, {})
    return best[1], best[2], best[3]


def write_rlp(path, I, H, steady, NLP, N0, M0, flows, y_by_site):
    """写规则文件。flows: (k, i, s, D, r, x) 数组（稳态时为空）；y_by_site: {s: {k: 速率 LP 班率}}。返回有理近似的诊断。"""
    nk, ns = len(I.T), len(I.sites)
    lines = [f"H {H}", f"steady {int(steady)}",
             "NLPq " + " ".join(str(int(round(1e6 * v))) for v in NLP),
             "N0 " + " ".join(map(str, N0)), "M0 " + " ".join(map(str, M0))]
    Qs, P, errs = [], {}, []
    for s in range(ns):
        Q, pp, err = rationalize(I, s, y_by_site.get(s, {})); Qs.append(Q); errs.append(err)
        for k, v in pp.items(): P[(k, s)] = v
    lines.append("Q " + " ".join(map(str, Qs)))
    lines.append("CP " + " ".join(str(40 * I.L[k] // I.d[k]) for k in range(nk)))
    for k in range(nk):
        for s in I.elig[k]:
            lines.append(f"p {k} {s} {P[(k, s)]} {CLIP_SLOTS * Qs[s] // I.d[k]}")
    for k in range(nk):
        L = I.L[k]
        w = [max(1, int(round(np.mean([I.step(k, (i - 1) % L, s)[1] for s in I.elig[k]])))) for i in range(L)]
        lines.append(f"W {k} " + " ".join(map(str, w)))
    if H > 0:
        Fks = np.zeros((nk, ns, H)); Fki = {k: np.zeros((I.L[k], H)) for k in range(nk)}
        for (k, i, s, D, r, x) in flows:
            k, i, s, D = int(k), int(i), int(s), int(D)
            if D < H: Fks[k, s, D] += x; Fki[k][i, D] += x
        for k in range(nk):
            for s in I.elig[k]:
                lines.append(f"F {k} {s} " + " ".join(str(int(round(1e6 * v))) for v in np.cumsum(Fks[k, s])))
            for i in range(I.L[k]):
                lines.append(f"FP {k} {i} " + " ".join(str(int(round(1e6 * v))) for v in np.cumsum(Fki[k][i])))
    with open(path, "w") as f: f.write("\n".join(lines) + "\n")
    return dict(Q=Qs, max_err=max(errs))


def write_planf(path, H, N, flows):
    """与 make_planf.py 相同：#N LP 台数；k i D 流量（提前大修弧、D < H/2，按 (k,i,D) 合并）。"""
    agg = {}
    for (k, i, s, D, r, x) in flows:
        if int(r) == 1 and int(D) < H // 2:
            key = (int(k), int(i), int(D)); agg[key] = agg.get(key, 0.0) + float(x)
    with open(path, "w") as f:
        f.write("#N " + " ".join(repr(float(v)) for v in N) + "\n")
        for (k, i, D) in sorted(agg): f.write(f"{k} {i} {D} {repr(agg[(k, i, D)])}\n")
    return len(agg)


def write_lpg(path, H, N, M, flows):
    """与 make_planf.py 写的 .lpg.txt 相同：#N0 ⌊N⌋、#M0 ⌊M⌋；z k s = D < H/2 的值岗流在站 s 内的占比。"""
    dz = {}
    for (k, i, s, D, r, x) in flows:
        if int(D) < H // 2 and x > 1e-9: dz[(int(k), int(s))] = dz.get((int(k), int(s)), 0.0) + float(x)
    tot = {}
    for (k, s), x in dz.items(): tot[s] = tot.get(s, 0.0) + x
    with open(path, "w") as f:
        f.write("#N0 " + " ".join(str(int(math.floor(x + 1e-6))) for x in N) + "\n")
        f.write("#M0 " + " ".join(str(int(math.floor(x + 1e-6))) for x in M) + "\n")
        for (k, s), x in dz.items(): f.write(f"z {k} {s} {x / tot[s]:.4f}\n")


def make(name, out=OUT, t2j=T2J, t2d=T2D, mult=4):
    """生成三份规则文件和一份 planf。"""
    os.makedirs(out, exist_ok=True)
    I = load_train(name)
    rate = json.load(open(os.path.join(t2j, f"{name}__rate__m0.json")))
    y = {}
    for key, v in rate["y"].items():
        k, s = map(int, key.split(",")); y.setdefault(s, {})[k] = v
    ceil = lambda v: int(math.ceil(v - 1e-6))
    info = {}
    for tag, kind in [("fresh_rlp", "frR"), ("fresh_noreset", "fr")]:
        d = np.load(os.path.join(t2d, f"{name}__{kind}__m{mult:g}.npz"))
        H = int(d["H"]); N = d["N"]; M = d["M"]; fl = d["duty_flow"]
        info[tag] = write_rlp(os.path.join(out, f"{name}__{tag}.rlp"), I, H, False, N, [ceil(v) for v in N], [ceil(v) for v in M], fl, y)
        if kind == "frR":
            info["planf_entries"] = write_planf(os.path.join(out, f"{name}__fresh_rlp.planf"), H, N, fl)
            write_lpg(os.path.join(out, f"{name}.lpg.txt"), H, N, M, fl)     # lpgplan ES 参照用
        # 诊断变体（不是规格里的 R_LP）：缩放分母用 ⌈N_LP⌉ 而不是 N_LP，见 报告/T3_松弛规则.md
        N0c = [ceil(v) for v in N]
        write_rlp(os.path.join(out, f"{name}__{tag}_c.rlp"), I, H, False, N0c, N0c, [ceil(v) for v in M], fl, y)
        if kind == "frR": write_planf(os.path.join(out, f"{name}__fresh_rlp_c.planf"), H, [float(v) for v in N0c], fl)
    info["steady"] = write_rlp(os.path.join(out, f"{name}__steady.rlp"), I, 0, True, rate["N"],
                                    [ceil(v) for v in rate["N"]], [ceil(v) for v in rate["M"]], np.zeros((0, 6)), y)
    return info


if __name__ == "__main__":
    for n in sys.argv[1:]:
        print(n, json.dumps(make(n)))
