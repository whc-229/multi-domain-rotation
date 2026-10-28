"""T2 单个作业：一个 (训练实例, 下界种类, H 倍数)。结果写进自己的文件，已存在就跳过（续跑）。

种类：
- rate：速率下界 R_rate（prior 的 tool_lb，LP 与整数台数），并存 LP 的机型—站点班率 y[k,s]（T3 外推用）；
- fr：全新启动 LP_fr(H)（freshlb.fresh_lb 的模型，不许提前大修）；
- frR：同上，允许提前大修。
fr / frR 的求解：整体 LP（code/fresh_full.py：freshlb 建模，highspy 内点法 + crossover），然后
  1) 按 freshlb 的列、行枚举顺序重建索引，并与截下的矩阵逐块比对（确认每列、每行的含义）；
  2) 对偶证书：把整体 LP 的对偶代入列生成的定价（code/cg_fresh.py 的最长路），检查所有单机路径的约化成本 ≥ 0、
     工具弧 / 工具等待弧 / M_j 的约化成本 ≥ 0。都成立时 Σ π_{s,t} K_s 是不依赖求解器的有效下界；
  3) 保存覆盖对偶 π（槽权重）、耦合对偶、按 (k, i, s, D, 是否提前大修) 聚合的值岗流、N、M。
用法：python3 code/t2_job.py <实例名> <rate|fr|frR> <H 倍数>"""
import os, sys, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, U, ToolInst, load_train  # noqa: E402
import tools_model  # noqa: E402
from freshlb import S_max  # noqa: E402

OUT = os.environ.get("T2_OUT", os.path.join(ROOT, "results", "T2", "jobs"))        # 环境变量只用于开发期把试算结果写到临时目录
DUAL = os.environ.get("T2_DUAL", os.path.join(ROOT, "results", "T2_duals"))


def job_name(inst, kind, mult):
    return f"{inst}__{kind}__m{mult:g}"


def rate_job(I):
    """R_rate：tool_lb 的 LP 值与整数台数版本；截下 LP 解里的班率 y[k,s]。"""
    TI = ToolInst(I)
    cap = {}
    real = tools_model.milp
    def spy(c, **kw):
        r = real(c, **kw); cap["x"] = r.x; return r
    tools_model.milp = spy
    try:
        t0 = time.time()
        lp, N, M = tools_model.tool_lb(TI, integer=False)
        x = cap["x"]
        yv = [(k, s) for k in range(len(I.T)) for s in I.elig[k]]
        y = {f"{k},{s}": float(x[a]) for a, (k, s) in enumerate(yv)}
        lpi, Ni, Mi = tools_model.tool_lb(TI, integer=True)
        sec = time.time() - t0
    finally:
        tools_model.milp = real
    return dict(R_rate=float(lp), R_rate_int=float(round(lpi)), N=[float(v) for v in N], M=[float(v) for v in M],
                N_int=[int(round(v)) for v in Ni], M_int=[int(round(v)) for v in Mi], y=y, sec=round(sec, 2))


def fresh_index(I, H, reset):
    """按 freshlb.fresh_lb 的循环顺序重建：机器人弧列的含义、行的顺序。返回 dict。"""
    TI = ToolInst(I)
    nk, ns, nt = len(I.T), len(I.sites), TI.nt
    v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]; c0 = max(v0)
    nrn = sum(I.L[k] * (H - v0[k]) for k in range(nk))
    cov = [(s, t) for s in range(ns) for t in range(c0 + I.lead(s), H)]
    duty_key = {}; arcs = []          # arcs: (列, k, i, v, s, r)，r = 1 为提前大修弧；等待弧 s = −1
    col = 0
    for k in range(nk):
        L = I.L[k]
        for i in range(L):
            for v in range(v0[k], H):
                arcs.append((col, k, i, v, -1, 0)); col += 1
                for s in I.elig[k]:
                    ld, nx = I.step(k, i, s)
                    opts = [0] + ([1] if (reset and i != L - 1) else [])
                    D = v + ld
                    for r in opts:
                        arcs.append((col, k, i, v, s, r))
                        if D < H and (k, s, D) not in duty_key: duty_key[(k, s, D)] = len(duty_key)
                        col += 1
    return dict(nrn=nrn, ntn=nt * H, cov=cov, duty_key=duty_key, arcs=arcs, n_rob=col, v0=v0, c0=c0)


def verify_index(I, H, idx, A):
    """用重建的索引生成覆盖块、耦合块里机器人弧的部分，与截下的矩阵比对。返回不一致的非零元个数。"""
    from scipy.sparse import coo_matrix
    TI = ToolInst(I)
    nA = idx["nrn"] + idx["ntn"]; ncov = len(idx["cov"]); nq = len(idx["duty_key"])
    covpos = {st: a for a, st in enumerate(idx["cov"])}
    R, C, V = [], [], []
    for (col, k, i, v, s, r) in idx["arcs"]:
        if s < 0: continue
        D = v + I.lead(s)
        for q in range(I.d[k]):
            a = covpos.get((s, D + q))
            if a is not None: R.append(nA + a); C.append(col); V.append(1.0)
        if D < H: R.append(nA + ncov + idx["duty_key"][(k, s, D)]); C.append(col); V.append(-float(TI.R[k]["q"]))
    E = coo_matrix((V, (R, C)), shape=(nA + ncov + nq, idx["n_rob"])).tocsr()
    got = A.tocsr()[:nA + ncov + nq, :idx["n_rob"]].tocsr()
    got = got[nA:, :]; E = E[nA:, :]
    d = (got - E); d.eliminate_zeros()
    return int(d.nnz), int(A.shape[0] == nA + ncov + nq)


def certificate(I, H, reset, idx, y_full, obj):
    """对偶证书：把整体 LP 的对偶映射到列生成主问题的行，检查所有列的约化成本。"""
    from cg_fresh import FreshCG
    cg = FreshCG(I, H, reset)
    nA = idx["nrn"] + idx["ntn"]; ncov = len(idx["cov"])
    y = np.zeros(cg.nrow)
    assert list(cg.cov_row) == idx["cov"]
    y[:cg.nc] = y_full[nA:nA + ncov]
    for key, row in cg.duty_key.items(): y[row] = y_full[nA + ncov + idx["duty_key"][key]]
    y[cg.tool_row0:] = y_full[idx["nrn"]:nA]
    pi, mu = cg.split_duals(y)
    path_viol = []
    for k in range(cg.nk):
        W = cg.weights(k, pi, mu); val, path = cg.price(k, W)
        path_viol.append(float(val - cg.cost[k]))
    phi = y[cg.tool_row0:].reshape(cg.nt, H)
    tool_viol = max((y[cg.duty_key[(k, s, D)]] - (phi[j, a] - (phi[j, b] if b < H else 0.0))
                     for (_, k, s, D, j, a, b) in cg.tool_arc), default=0.0)
    wait_viol = float(np.max(phi[:, 1:] - phi[:, :-1])) if H > 1 else 0.0
    M_viol = float(np.max(phi[:, 0] - 1.0))
    dual_obj = float(y[:cg.nc] @ cg.rhs[:cg.nc])
    return dict(path_viol=max(path_viol), path_viol_k=path_viol, tool_viol=float(tool_viol), wait_viol=wait_viol, M_viol=M_viol,
                dual_obj=dual_obj, dual_gap=dual_obj - obj, pi=pi, mu=mu, cg=cg)


def fresh_job(I, H, reset, dual_path, perm=None, crossover="on", tl=6 * 3600):
    """整体 LP + 对偶证书 + 保存对偶与值岗流。perm / crossover / tl 只在重跑（--tag）时改动，见 main。"""
    from fresh_full import fresh_full
    keep = {}
    t0 = time.time()
    r = fresh_full(I, H, reset=reset, solver="ipm", time_limit=tl, keep=keep, perm=perm, crossover=crossover)
    sec = time.time() - t0
    out = dict(obj=r["obj"], N=r.get("N"), M=r.get("M"), lp_status=keep.get("status"), n=r.get("n"), sec=round(sec, 1))
    if r["obj"] is None or keep.get("row_dual") is None: return out
    obj = float(keep["c"] @ keep["x"])
    out["obj_full"] = obj
    idx = fresh_index(I, H, reset)
    bad, rows_ok = verify_index(I, H, idx, keep["A"])
    out["index_mismatch"] = bad; out["rows_ok"] = rows_ok
    t1 = time.time()
    cert = certificate(I, H, reset, idx, keep["row_dual"], obj)
    out.update({k: cert[k] for k in ("path_viol", "path_viol_k", "tool_viol", "wait_viol", "M_viol", "dual_obj", "dual_gap")})
    out["cert_sec"] = round(time.time() - t1, 1)
    # 值岗流
    x = keep["x"]; fl = {}
    for (col, k, i, v, s, rr) in idx["arcs"]:
        if s < 0 or x[col] <= 1e-9: continue
        key = (k, i, s, v + I.lead(s), rr); fl[key] = fl.get(key, 0.0) + float(x[col])
    keys = sorted(fl)
    duty_flow = np.array([list(kk) + [fl[kk]] for kk in keys]) if keys else np.zeros((0, 6))
    cg = cert["cg"]
    yv = np.array([[k, s, D, float(cert["mu"][(k, s)][D])] for (k, s, D) in cg.duty_key]) if cg.duty_key else np.zeros((0, 4))
    np.savez_compressed(dual_path, pi=cert["pi"], cov_first=np.array(cg.cov_first), mu=yv, duty_flow=duty_flow,
                        N=np.array(r["N"]), M=np.array(r["M"]), H=H, reset=int(reset), obj=obj, v0=np.array(idx["v0"]), c0=idx["c0"])
    out["n_flow_arcs"] = len(keys)
    out["resets_flow"] = float(sum(fl[kk] for kk in keys if kk[4] == 1))
    return out


def main():
    """用法：t2_job.py <实例> <rate|fr|frR> <倍数> [--tag 名 --perm 种子 --crossover on|off --tl 秒]。
    --tag 用于失败作业的重跑：结果写到 <原名>__<tag>.json（不覆盖原结果），汇总时取求解成功的那一份。"""
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("inst"); ap.add_argument("kind"); ap.add_argument("mult", type=float)
    ap.add_argument("--tag"); ap.add_argument("--perm", type=int); ap.add_argument("--crossover", default="on")
    ap.add_argument("--tl", type=float, default=6 * 3600)
    a = ap.parse_args()
    inst, kind, mult = a.inst, a.kind, a.mult
    os.makedirs(OUT, exist_ok=True); os.makedirs(DUAL, exist_ok=True)
    name = job_name(inst, kind, mult) + (f"__{a.tag}" if a.tag else "")
    path = os.path.join(OUT, name + ".json")
    if os.path.exists(path):
        print("已完成，跳过", name); return
    I = load_train(inst)
    S = S_max(I); H = int(mult * S)
    rec = dict(inst=inst, m=I.meta["m"], family=I.meta["family"], kind=kind, mult=mult, S=S, H=H)
    t0 = time.time()
    if kind == "rate":
        rec.update(rate_job(I))
    else:
        rec.update(fresh_job(I, H, kind == "frR", os.path.join(DUAL, name + ".npz"), a.perm, a.crossover, a.tl))
        if a.tag: rec.update(tag=a.tag, perm=a.perm, crossover=a.crossover, tl=a.tl)
    rec["wall"] = round(time.time() - t0, 1)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)      # 原子写入：中断时不会留下半个结果文件
    print(json.dumps({k: v for k, v in rec.items() if k not in ("y", "path_viol_k")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
