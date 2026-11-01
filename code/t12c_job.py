"""T12c 单个作业（一个 m ∈ {2, 3} 的训练实例）：核对作为上界的 ES 解能否用全精度权重复现。
起因（T12b）：rost5 es 输出的权重 w 只有 4 位小数，内部用的是全精度；用 w 重新判定可能与 es 自己的判定不同。
T4 / T4b 的上界（以及 T8 的"ES 饱和值"、T9 不提前大修的上界）里有一部分只来自 ES 解。

候选解与 T4 / T4b 汇总相同：T3 的 R_LP 各变体与 es（10 代）、T4 的 esub（30 代）、T4b 的 es100（100 代，种子 0–4），
可行 = 严格证明或 60 万格复核 0 缺岗；不允许提前大修的场景只取提前大修次数为 0 的解。
对每个场景（不允许 / 允许提前大修）和每个上界口径（T4：T3 + T4；T4b：再加 es100），按资源数从小到大逐个看候选解：
- 非 ES 解（R_LP 规则文件，确定性复现）：直接采用；
- ES 解：用当前 rost5（输出全精度权重 w17）按原命令、原种子补跑，加 certadd=why（这些作业在 T6 改补台规则之前跑的；
  试算 4 个作业逐字段复现，不加时 1 / 4 不同）。
  先核对补跑与原作业相同（R、状态、台数、复核结果、提前大修次数、4 位小数的 w、逐代最好值）；
  再用 w17 跑 rost5 eval（60 万格判环，同一组选项）：缺岗、判环、周期、exact、提前大修次数、值岗次数与原作业相同，算"能复现"。
  补跑与原作业不同时（不应出现），退而用原作业的 4 位小数 w 重新判定，相同也算能复现，并单独记下。
  另记下用 4 位小数 w 重新判定是否相同（只记录）。
- 第一个能复现的解的资源数就是新上界；同资源数时先看非 ES 解，再按 T3 → T4 → T4b、种子从小到大。
每个 ES 解的补跑结果缓存在 results/T12c/es/<原文件名>（原子写入），两个场景、两个口径共用。
结果写进 results/T12c/jobs/<实例>.json（原子写入，已存在就跳过；运行时加 fcntl 锁）。
用法：python3 code/t12c_job.py <实例>"""
import os, sys, json, glob, time, subprocess, fcntl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")
OUT = os.path.join(ROOT, "results", "T12c", "jobs")
ESD = os.path.join(ROOT, "results", "T12c", "es")
RLP_KINDS = ("fresh_rlp", "fresh_noreset", "fresh_rlp_c", "fresh_noreset_c")
ES_KINDS = ("es_curdef", "es_lpgplan", "esub_curdef", "esub_lpgplan", "es100_curdef", "es100_lpgplan")
SAME_KEYS = ("R", "status", "Nk", "Mj", "verify_gaps", "verify_cert", "verify_period", "exact", "resets", "w", "hist")
SRC_ORDER = {"T3": 0, "T4": 1, "T4b": 2}


def feasible(d):
    """严格证明或 60 万格复核 0 缺岗。"""
    return d.get("R") is not None and d["R"] >= 0 and (d.get("proven") == 1 or d.get("verify_gaps") == 0)


def candidates(n):
    """该实例的全部可行候选解：[(R, 是否 ES, 来源, 种子, 文件, 记录)]。"""
    out = []
    for src, pat in (("T3", f"results/T3/jobs/{n}__*.json"), ("T4", f"results/T4/jobs/{n}__esub_*.json"), ("T4b", f"results/T4b/jobs/{n}__es100_*.json")):
        for f in sorted(glob.glob(os.path.join(ROOT, pat))):
            d = json.load(open(f))
            if d.get("kind") not in RLP_KINDS + ES_KINDS or not feasible(d): continue
            out.append((d["R"], d["kind"] in ES_KINDS, src, d.get("seed", -1), f, d))
    return out


def rejudge(d, w):
    """用给定权重对 d 的台数跑 rost5 eval（60 万格判环，与 es 的复核相同的选项），返回 eval 的输出。"""
    a = d["cmd"].split()
    ev = [R5, "eval", a[2], "600000", "1"] + [str(x) for x in d["Nk"]] + ["|"] + [str(x) for x in d["Mj"]] + ["|"] + [repr(float(x)) for x in w] + ["|"] + a[11:] + ["curdef"]
    return json.loads(subprocess.run(ev, capture_output=True, text=True, check=True).stdout)


def same_verdict(d, v):
    """重新判定的结果与原作业的复核是否相同。"""
    return (v["gaps"] == d["verify_gaps"] and v["cert"] == d["verify_cert"] and v["resets"] == d["resets"] and v["duties"] == d["duties"]
            and (d["verify_cert"] != 1 or (v["period"] == d["verify_period"] and v["exact"] == d["exact"])))


def check_es(f, d):
    """补跑一个 ES 解并重新判定（有缓存就直接读），返回记录。"""
    cp = os.path.join(ESD, os.path.basename(f))
    if os.path.exists(cp): return json.load(open(cp))
    t0 = time.time()
    r = json.loads(subprocess.run(d["cmd"].split() + ["certadd=why"], capture_output=True, text=True, check=True).stdout)
    diff = [k for k in SAME_KEYS if d.get(k) != r.get(k)]
    v17 = rejudge(d, r["w17"]) if not diff else None
    v4 = rejudge(d, d["w"])
    rec = dict(file=os.path.relpath(f, ROOT), inst=d["inst"], kind=d["kind"], seed=d.get("seed"), R=d["R"], resets=d["resets"],
               same_as_orig=int(not diff), diff=diff, w17=r.get("w17"),
               rejudge_w17=None if v17 is None else int(same_verdict(d, v17)), rejudge_w4=int(same_verdict(d, v4)),
               v17={k: v17[k] for k in ("gaps", "cert", "period", "exact", "resets", "duties")} if v17 else None,
               v4={k: v4[k] for k in ("gaps", "cert", "period", "exact", "resets", "duties")},
               orig={k: d[k] for k in ("verify_gaps", "verify_cert", "verify_period", "exact", "resets", "duties")}, wall=round(time.time() - t0, 1))
    rec["reproducible"] = int(rec["rejudge_w17"] == 1 if not diff else rec["rejudge_w4"] == 1)
    rec["via"] = "w17" if not diff else "w4"
    tmp = cp + ".tmp"
    with open(tmp, "w") as fh: fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, cp)
    return rec


def walk(cands):
    """按资源数从小到大找第一个能复现的解，返回 (上界, 来源, 走过的 ES 解)。"""
    order = sorted(cands, key=lambda c: (c[0], c[1], SRC_ORDER[c[2]], c[5]["kind"], c[3]))
    tried = []
    for R, is_es, src, seed, f, d in order:
        if not is_es: return R, dict(kind=d["kind"], src=src, file=os.path.relpath(f, ROOT)), tried
        rec = check_es(f, d)
        tried.append(dict(file=rec["file"], R=R, reproducible=rec["reproducible"], via=rec["via"]))
        if rec["reproducible"]: return R, dict(kind=d["kind"], src=src, seed=seed, file=rec["file"], via=rec["via"]), tried
    return None, None, tried


def main():
    """两个场景 × 两个口径（T4、T4b）逐个走候选解，原子写入。"""
    n = sys.argv[1]
    os.makedirs(OUT, exist_ok=True); os.makedirs(ESD, exist_ok=True)
    path = os.path.join(OUT, f"{n}.json")
    if os.path.exists(path): print("已完成，跳过", n); return
    lk = open(path + ".lock", "w")
    try: fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError: print("另一个进程正在做", n); sys.exit(4)
    if os.path.exists(path): print("已完成，跳过", n); return
    t0 = time.time(); C = candidates(n); res = dict(inst=n)
    for rs in (0, 1):
        cs = [c for c in C if rs == 1 or c[5].get("resets", 0) == 0]
        for scope, keep in (("T4", ("T3", "T4")), ("T4b", ("T3", "T4", "T4b"))):
            sub = [c for c in cs if c[2] in keep]
            old = min((c[0] for c in sub), default=None)
            ub, src, tried = walk(sub)
            res[f"r{rs}_{scope}"] = dict(UB_old=old, UB_new=ub, src=src, tried=tried)
    res["wall"] = round(time.time() - t0, 1)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh: fh.write(json.dumps(res, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: (v if not isinstance(v, dict) else (v["UB_old"], v["UB_new"], len(v["tried"]))) for k, v in res.items()}, ensure_ascii=False))


if __name__ == "__main__":
    main()
