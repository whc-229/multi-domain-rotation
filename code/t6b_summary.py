"""T6b 汇总：松弛规则 fresh_rlp_w 与 ES lpgplan（新补台）在各规模上离 ⌈LP_fr^Rw(4S)⌉ 的差距。
- fresh_rlp_w：m ≤ 6 取 T6 的结果，m ∈ {8, 10, 15, 20} 是 T6b 新跑的（同一脚本、同一目录 results/T6/rlp/）；
- ES lpgplan（新补台）：m ≤ 6 取 T6，m ∈ {8, 10} 是 T6b 新跑的（results/T6/es/，种子 0–2）；m ≥ 15 不跑；
- 可行 = 严格证明或 60 万格复核 0 缺岗；差距只在可行解上算，(R − ⌈LP⌉) / ⌈LP⌉。
输出 results/T6b_large.jsonl（逐实例）并打印按 (m, 族) 的 markdown 表。
用法：python3 code/t6b_summary.py"""
import os, sys, json, math
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names  # noqa: E402

RLP = os.path.join(ROOT, "results", "T6", "rlp"); ES = os.path.join(ROOT, "results", "T6", "es")


def jl(p):
    """读 json；不存在返回 None。"""
    return json.load(open(p)) if os.path.exists(p) else None


def feas(d):
    """可行时返回资源数，否则 None。"""
    if not d or d.get("R") is None or d["R"] < 0: return None
    return d["R"] if (d.get("proven") == 1 or d.get("verify_gaps") == 0) else None


def rows():
    """逐实例记录：LP 下界、rlp_w 与 3 个 ES 种子的资源数、证明情况、用时。"""
    lp = {r["inst"]: r["LP_frRw"] for r in map(json.loads, open(os.path.join(ROOT, "results", "T6_bounds.jsonl")))}
    out = []
    for n in sorted(train_names(), key=lambda n: (int(n[1:3]), n)):
        w = jl(os.path.join(RLP, f"{n}__fresh_rlp_w.json"))
        es = [jl(os.path.join(ES, f"{n}__es_lpgplan__s{s}.json")) for s in range(3)]
        if w is None and not any(es): continue
        lb = int(math.ceil(lp[n] - 1e-6))
        out.append(dict(inst=n, m=int(n[1:3]), family=n.split("-")[1], LB=lb,
                        rlp_w=feas(w), rlp_w_ran=int(w is not None), rlp_w_proven=int(bool(w and w.get("proven") == 1)), rlp_w_sec=w and w.get("sec"),
                        es=[feas(d) for d in es], es_ran=[int(d is not None) for d in es], es_proven=[int(bool(d and d.get("proven") == 1)) for d in es],
                        es_sec=[d and d.get("sec") for d in es]))
    return out


def mean(xs):
    """忽略 None 的平均；全空返回 None。"""
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def fm(x, nd=1):
    """格式化。"""
    return "—" if x is None else f"{x:.{nd}f}"


def main():
    """写 jsonl，打印表。"""
    rs = rows()
    tmp = os.path.join(ROOT, "results", "T6b_large.jsonl.tmp")
    with open(tmp, "w") as f:
        for r in rs: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(tmp, os.path.join(ROOT, "results", "T6b_large.jsonl"))
    by = defaultdict(list)
    for r in rs: by[(r["m"], r["family"])].append(r)
    L = ["| m | 族 | 实例 | ⌈LP_fr^Rw⌉ | fresh_rlp_w 离下界 % | 可行 / 证明 / 跑完 | 用时 s | ES 离下界 %（种子平均） | ES 可行 / 证明 / 跑完（种子） | ES 用时 s/种子 | 同实例 rlp_w − ES 平均（个） |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for (m, f), xs in sorted(by.items()):
        gw = mean([100 * (x["rlp_w"] - x["LB"]) / x["LB"] for x in xs if x["rlp_w"] is not None])
        ges = mean([100 * (e - x["LB"]) / x["LB"] for x in xs for e in x["es"] if e is not None])
        nes = sum(sum(x["es_ran"]) for x in xs)
        # 同实例配对：rlp_w 与该实例可行 ES 种子的平均
        dd = mean([x["rlp_w"] - mean(x["es"]) for x in xs if x["rlp_w"] is not None and mean(x["es"]) is not None])
        L.append(f"| {m} | {f} | {len(xs)} | {fm(mean([x['LB'] for x in xs]), 1)} | {fm(gw)} | "
                 f"{sum(x['rlp_w'] is not None for x in xs)} / {sum(x['rlp_w_proven'] for x in xs)} / {sum(x['rlp_w_ran'] for x in xs)} | "
                 f"{fm(mean([x['rlp_w_sec'] for x in xs]))} | {fm(ges)} | "
                 + (f"{sum(e is not None for x in xs for e in x['es'])} / {sum(sum(x['es_proven']) for x in xs)} / {nes}" if nes else "—")
                 + f" | {fm(mean([s for x in xs for s in x['es_sec']]), 0)} | {fm(dd, 2)} |")
    print("\n".join(L))


if __name__ == "__main__":
    main()
