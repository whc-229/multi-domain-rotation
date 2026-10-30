"""T9 汇总：整数机队松弛 IFR / IFR^Rw 的下界。输出 results/T9_intfleet.jsonl 和 markdown 表（打印）。
对照：
- ⌈LP(2S)⌉：同一模型的 LP 松弛（T2 的 fr m2、frR m2；H = 2S 时 LP_fr^Rw = LP_fr^R，见 t9_job.py）；
- 原最好下界：不提前大修 = max(⌈LP_fr(4S)⌉, T4 的 ILP 下界)；启动期提前大修 = ⌈LP_fr^Rw(4S)⌉（T6，新口径的主下界）；
- 已知最好上界：
  * 不提前大修 = T4 / T4b 的最好上界、T8 的 RDD 与默认规则可行解、T6 的 fresh_noreset_u 中的最小者；
  * 启动期提前大修 = 上一行（不提前大修的解在这个场景下也可行）与 T6 / T6b 的 fresh_rlp_w（只在窗口内提前大修）中的最小者。
    T4 / T4b 的"允许提前大修"上界来自提前大修不限时刻的解，在"只许启动期提前大修"的场景下不一定可行，只作参照列 UB_any，不进差距。
新下界 = max(原最好下界, IFR 下界)。差距 = 上界 − 下界。
整数台数对比：IFR 有可行解时，逐型比较 IFR 解的 N_k、M_j 与 LP(2S) 的小数值。
用法：python3 code/t9_summary.py"""
import os, sys, json, glob, math
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names, load_train  # noqa: E402

J = os.path.join(ROOT, "results", "T9", "jobs")


def jl(p):
    """读 json；不存在返回 None。"""
    return json.load(open(p)) if os.path.exists(p) else None


def ceil(x):
    """向上取整（容差 1e-6）。"""
    return int(math.ceil(x - 1e-6))


def t8_best(n):
    """T8 里该实例所有 RDD 可行解（严格证明或 60 万格复核 0 缺岗）的最少资源数；没有返回 None。"""
    best = None
    for f in glob.glob(os.path.join(ROOT, "results", "T8", "jobs", f"{n}__*.json")):
        d = json.load(open(f))
        for R, ok in ((d.get("R"), d.get("proven") == 1 or d.get("verify_gaps") == 0),
                      (d.get("base_R"), d.get("base_proven") == 1 or d.get("base_verify_gaps") == 0)):
            if R is not None and R >= 0 and ok: best = R if best is None else min(best, R)
    return best


def feas_json(p):
    """规则结果文件可行（严格证明或 60 万格复核 0 缺岗）时返回资源数，否则 None。"""
    d = jl(p)
    if not d or d.get("R") is None or d["R"] < 0: return None
    return d["R"] if (d.get("proven") == 1 or d.get("verify_gaps") == 0) else None


def ub_noreset(n, T4b):
    """不提前大修场景的已知最好上界（见文件头）。"""
    c = [T4b[(n, 0)]["UB_new"] if (n, 0) in T4b else None, t8_best(n),
         feas_json(os.path.join(ROOT, "results", "T6", "rlp", f"{n}__fresh_noreset_u.json"))]
    c = [x for x in c if x is not None]
    return min(c) if c else None


def rows():
    """逐 (实例, 场景) 的记录。"""
    T4 = {(r["inst"], r["reset"]): r for r in map(json.loads, open(os.path.join(ROOT, "results", "T4_exact.jsonl")))}
    T4b = {(r["inst"], r["reset"]): r for r in map(json.loads, open(os.path.join(ROOT, "results", "T4b_gap.jsonl")))}
    out = []
    for n in sorted(train_names(m_list=[2, 3, 4, 6]), key=lambda n: (int(n[1:3]), n)):
        for kind, rs in (("fr", 0), ("frRw", 1)):
            d = jl(os.path.join(J, f"{n}__ifr_{kind}.json"))
            if d is None: continue
            lp2 = jl(os.path.join(ROOT, "results", "T2", "jobs", f"{n}__{'fr' if rs == 0 else 'frR'}__m2.json"))
            if rs == 0:
                old_lb = T4[(n, 0)]["LB"] if (n, 0) in T4 else ceil(json.load(open(os.path.join(ROOT, "results", "T2", "jobs", f"{n}__fr__m4.json")))["obj"])
            else:
                w = jl(os.path.join(ROOT, "results", "T6", "jobs", f"{n}__frRw__m4.json"))
                old_lb = ceil(w["obj"]) if w and w.get("obj") is not None else None
            ub = ub_noreset(n, T4b)
            ub_any = T4b[(n, 1)]["UB_new"] if rs == 1 and (n, 1) in T4b else None
            if rs == 1:
                w = feas_json(os.path.join(ROOT, "results", "T6", "rlp", f"{n}__fresh_rlp_w.json"))
                if w is not None: ub = w if ub is None else min(ub, w)
            lb_raw = d.get("LB")
            # MIP 的松弛就是同一个 LP，所以有效下界 = max(MIP 对偶界, ⌈LP(2S)⌉)；对偶界低于 LP 说明限时内根节点 LP 都没解完
            lb_i = max(x for x in (lb_raw, ceil(lp2["obj"]) if lp2 else None) if x is not None) if (lb_raw is not None or lp2) else None
            new_lb = max(x for x in (old_lb, lb_i) if x is not None) if (old_lb is not None or lb_i is not None) else None
            rec = dict(inst=n, m=d["m"], family=d["family"], reset=rs, LP2=lp2["obj"] if lp2 else None, LP2_ceil=ceil(lp2["obj"]) if lp2 else None,
                       IFR_LB=lb_i, IFR_LB_raw=lb_raw, root_unsolved=int(lb_raw is not None and lp2 is not None and d.get("bound") is not None and d["bound"] < lp2["obj"] - 1e-6), IFR_bound=d.get("bound"), IFR_obj=d.get("obj"), IFR_status=d.get("status"), IFR_sec=d.get("sec"),
                       old_LB=old_lb, new_LB=new_lb, UB=ub, UB_any=ub_any, IFR_N=d.get("N"), IFR_M=d.get("M"),
                       LP2_N=lp2.get("N") if lp2 else None, LP2_M=lp2.get("M") if lp2 else None)
            out.append(rec)
    return out


def mean(xs):
    """忽略 None 的平均。"""
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def fm(x, nd=2):
    """格式化。"""
    return "—" if x is None else f"{x:.{nd}f}"


def main():
    rs = rows()
    with open(os.path.join(ROOT, "results", "T9_intfleet.jsonl.tmp"), "w") as f:
        for r in rs: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(os.path.join(ROOT, "results", "T9_intfleet.jsonl.tmp"), os.path.join(ROOT, "results", "T9_intfleet.jsonl"))
    L = ["| m | 族 | 场景 | 实例 | ⌈LP(2S)⌉ | IFR 有效下界 | IFR > ⌈LP(2S)⌉ 的实例 | 根节点未解完 | 原最好下界 | 新下界 | IFR 抬高原下界的实例 | 上界 | 原差距 | 新差距 | 解到最优 | 有整数可行解 | 平均用时 s |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    by = defaultdict(list)
    for r in rs: by[(r["m"], r["family"], r["reset"])].append(r)
    for (m, f, s), xs in sorted(by.items()):
        ok = [x for x in xs if x["UB"] is not None and x["new_LB"] is not None and x["old_LB"] is not None]
        L.append(f"| {m} | {f} | {'启动期提前大修' if s else '不提前大修'} | {len(xs)} | {fm(mean([x['LP2_ceil'] for x in xs]))} | {fm(mean([x['IFR_LB'] for x in xs]))} | "
                 f"{sum((x['IFR_LB'] or 0) > (x['LP2_ceil'] or 0) for x in xs)} | {sum(x['root_unsolved'] for x in xs)} | {fm(mean([x['old_LB'] for x in xs]))} | {fm(mean([x['new_LB'] for x in xs]))} | "
                 f"{sum((x['IFR_LB'] or 0) > (x['old_LB'] or 0) for x in xs)} | {fm(mean([x['UB'] for x in ok]))} | {fm(mean([x['UB'] - x['old_LB'] for x in ok]))} | "
                 f"{fm(mean([x['UB'] - x['new_LB'] for x in ok]))} | {sum(x['IFR_status'] == 'Optimal' for x in xs)} | {sum(x['IFR_N'] is not None for x in xs)} | "
                 f"{fm(mean([x['IFR_sec'] for x in xs]), 0)} |")
    print("\n".join(L))
    # 整数台数与小数台数
    print("\n### 整数解与 LP(2S) 的台数（有整数可行解的实例）\n")
    for r in rs:
        if r["IFR_N"] is None or r["LP2_N"] is None: continue
        I = load_train(r["inst"])
        dn = [(I.tn[k], round(r["IFR_N"][k] - r["LP2_N"][k], 2)) for k in range(len(I.tn))]
        dm = [(j, round(r["IFR_M"][j] - r["LP2_M"][j], 2)) for j in range(len(r["IFR_M"]))]
        print(f"- {r['inst']} {'R' if r['reset'] else ''}：IFR 目标 {fm(r['IFR_obj'])}、对偶界 {fm(r['IFR_bound'])}；机器人 Δ {dn}；工具 Δ {dm}")


if __name__ == "__main__":
    main()
