"""T10：m = 20 上 fresh_rlp_w 规则失败的诊断（只用训练实例）。输出 results/T10/diag.json，并打印 markdown 表。

依赖：
- results/T10/sizelog/<实例>.log：rost5 rlp … sizelog=… 的定台数过程记录（4 个起点倍数依次写在同一个文件里）；
- results/T6/rules/：T6 生成的 fresh_rlp_w 规则文件与 planf；results/T6_duals/：LP_fr^Rw(4S) 的值岗流；
- results/T6/rlp_diag/<实例>__fresh_rlp_w__trial60.json：T6b 不设停止条件的补台结果（补台最多的机型从这里看）。
三部分：
1. 缺岗在哪：每个起点倍数在证书阶段最后一组台数下，用 rost5 eval（gaplog）跑 6 万格，记下每个缺岗的 (时刻, 站, 就绪候选数)；
   按时刻分成 窗口内 [0, W)、窗口后到 LP 时域 [W, H)、H 之后；给出缺岗站可派的机型，及这些机型在 LP 值岗流里的份额。
2. 补台最多的机型：在 LP 计划里承担哪些站、能不能去缺岗站；并核对"给它加 1 台"对 6 万格缺岗数没有影响（中性补台）。
3. t ≥ H 后外推比例的有理近似：所有 m ≥ 8 训练实例的 Q_s 与最大误差（来自规则生成时的记录 <实例>.w.done）。
用法：python3 code/t10_diag.py"""
import os, sys, json, re, subprocess, tempfile
from collections import Counter, defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths, train_names  # noqa: E402
from t6_rerun import frrw_npz  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")
OUT = os.path.join(ROOT, "results", "T10")
RULES = os.path.join(ROOT, "results", "T6", "rules")
INSTS = ["m20-K1-01", "m20-K1-02", "m20-K12-00"]
FACS = [100, 115, 130, 150]


def parse_sizelog(path):
    """把定台数过程记录拆成 4 个起点倍数的段。每段返回 dict(fill_ok, pruned, certs=[(a, gaps, Nk, Mj)], end)。"""
    segs, cur = [], None
    for line in open(path):
        line = line.strip()
        if cur is None: cur = dict(fill_ok=False, pruned=None, certs=[], end=None)
        if line.startswith("fail_fill"): cur["end"] = "fail_fill"; segs.append(cur); cur = None; continue
        if line.startswith("pruned"):
            cur["fill_ok"] = True; cur["pruned"] = [json.loads(x) for x in re.findall(r"\[[^\]]*\]", line)]
        elif line.startswith("cert a="):
            a = int(re.search(r"a=(\d+)", line).group(1)); g = int(re.search(r"gaps=(\d+)", line).group(1))
            Nk, Mj = [json.loads(x) for x in re.findall(r"\[[^\]]*\]", line)]
            cur["certs"].append((a, g, Nk, Mj))
        elif line.startswith("fail_cert_stall"): cur["end"] = "fail_cert_stall"; segs.append(cur); cur = None
    return segs


def eval_gaps(name, Nk, Mj, T=60000, extra=()):
    """rost5 eval（判环）跑 T 格，返回 (eval 输出, 缺岗列表 [(t, s, 候选数)])。"""
    _, txt = inst_paths(name, "train")
    with tempfile.NamedTemporaryFile("r", suffix=".gap", dir=OUT, delete=False) as f: gp = f.name
    args = [R5, "eval", txt, str(T), "1"] + [str(x) for x in Nk] + ["|"] + [str(x) for x in Mj] + ["|", "|",
            f"rlp={RULES}/{name}__fresh_rlp_w.rlp", f"planf={RULES}/{name}__fresh_rlp_w.planf", f"gaplog={gp}"] + list(extra)
    r = json.loads(subprocess.run(args, capture_output=True, text=True, check=True).stdout)
    gaps = [tuple(int(x) for x in l.split()[:3]) for l in open(gp)]
    os.remove(gp)
    return r, gaps


def lp_shares(I, name):
    """LP_fr^Rw(4S) 值岗流：返回 (H, W, flow[k][s] 按时段 [0,W)、[W,H) 的班数, 每站总班数)。"""
    d = np.load(frrw_npz(name)); H, W = int(d["H"]), int(d["W"])
    fl = defaultdict(lambda: [0.0, 0.0])
    for (k, i, s, D, r, x) in d["duty_flow"]:
        fl[(int(k), int(s))][0 if D < W else 1] += float(x)
    return H, W, fl


def main():
    os.makedirs(OUT, exist_ok=True)
    rep = {}
    for n in INSTS:
        I = load_train(n); H, W, fl = lp_shares(I, n)
        segs = parse_sizelog(os.path.join(OUT, "sizelog", f"{n}.log"))
        rate = json.load(open(os.path.join(ROOT, "results", "T6", "jobs", f"{n}__rate.json")))["y_g"]
        out = dict(H=H, W=W, segs=[])
        for fac, sg in zip(FACS, segs):
            e = dict(fac=fac, end=sg["end"], fill_ok=sg["fill_ok"], cert_steps=len(sg["certs"]),
                     cert_gaps=[c[1] for c in sg["certs"]])
            if sg["certs"]:
                _, g, Nk, Mj = sg["certs"][-1]
                r, gaps = eval_gaps(n, Nk, Mj)
                e.update(Nk=Nk, Mj=Mj, eval=r, gaps=gaps)
                # 中性补台核对：最后一步每种资源各加 1 后的 6 万格缺岗数
                trial = []
                for idx in range(len(Nk) + len(Mj)):
                    N2, M2 = list(Nk), list(Mj)
                    if idx < len(Nk): N2[idx] += 1
                    else: M2[idx - len(Nk)] += 1
                    trial.append(eval_gaps(n, N2, M2, extra=())[0]["gaps"])
                e["trial_add_gaps"] = trial
            out["segs"].append(e)
        # 缺岗站可派的机型与 LP 份额
        sites = sorted({g[1] for e in out["segs"] for g in e.get("gaps", [])})
        out["gap_sites"] = {s: dict(types=[k for k in range(len(I.T)) if s in I.elig[k]],
                                   lp_win={k: round(fl[(k, s)][0], 2) for k in range(len(I.T)) if s in I.elig[k]},
                                   lp_after={k: round(fl[(k, s)][1], 2) for k in range(len(I.T)) if s in I.elig[k]},
                                   rate={k: round(rate.get(f"{k},{s}", 0.0), 4) for k in range(len(I.T)) if s in I.elig[k]})
                            for s in sites}
        # 补台最多的机型（T6b trial60 的结果相对起点）
        t60 = json.load(open(os.path.join(ROOT, "results", "T6", "rlp_diag", f"{n}__fresh_rlp_w__trial60.json")))
        inc = [a - b for a, b in zip(t60["Nk"], t60["N0"])]
        kmax = int(np.argmax(inc))
        out["trial60"] = dict(Nk=t60["Nk"], N0=t60["N0"], Mj=t60["Mj"], M0=t60["M0"], inc=inc, kmax=kmax,
                              kmax_sites={s: [round(fl[(kmax, s)][0], 2), round(fl[(kmax, s)][1], 2)] for s in I.elig[kmax]},
                              kmax_reaches_gap_sites=[s for s in sites if s in I.elig[kmax]])
        rep[n] = out
    # 有理近似（所有 m ≥ 8）
    qa = {}
    for n in train_names([8, 10, 15, 20]):
        p = os.path.join(RULES, f"{n}.w.done")
        if os.path.exists(p): d = json.loads(open(p).read()); qa[n] = dict(Q=d["Q"], max_err=d["max_err"])
    rep["_rational"] = qa
    with open(os.path.join(OUT, "diag.json.tmp"), "w") as f: json.dump(rep, f, ensure_ascii=False, indent=1)
    os.replace(os.path.join(OUT, "diag.json.tmp"), os.path.join(OUT, "diag.json"))
    # 打印
    for n in INSTS:
        o = rep[n]; print(f"\n## {n}  H={o['H']} W={o['W']}")
        for e in o["segs"]:
            print(f"- 起点 ×{e['fac']/100:.2f}：{e['end']}；证书阶段 {e['cert_steps']} 步，6 万格缺岗 {e['cert_gaps']}")
            if "gaps" in e:
                print(f"  最后台数 Nk={e['Nk']} Mj={e['Mj']}：缺岗 {[(t, s, c) for t, s, c in e['gaps']]}；各资源加 1 后缺岗 {e['trial_add_gaps']}")
        for s, v in o["gap_sites"].items(): print(f"  站 {s}：可派机型 {v['types']}，LP 窗口内班数 {v['lp_win']}，[W,H) 班数 {v['lp_after']}，速率 LP 班率 {v['rate']}")
        t = o["trial60"]; print(f"  trial60：补台 {t['inc']}，最多的是机型 {t['kmax']}；它在 LP 里承担的站（窗口内, [W,H)）{t['kmax_sites']}；能去的缺岗站 {t['kmax_reaches_gap_sites']}")
    errs = [v["max_err"] for v in qa.values()]
    print(f"\n有理近似：{len(qa)} 个 m ≥ 8 实例，Q_s 取值 {sorted(Counter(q for v in qa.values() for q in v['Q']).items())}，最大误差 {max(errs):.2e}")


if __name__ == "__main__":
    main()
