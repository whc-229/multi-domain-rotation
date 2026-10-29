"""T6 汇总：下界口径修正与规则重跑。输出 results/T6_bounds.jsonl、results/T6_rlp.jsonl 和 markdown 表（打印到标准输出）。
口径：
- 启动代价（不提前大修）= (LP_fr − R_rate_g) / R_rate_g；启动期提前大修 = (LP_fr^Rw − R_rate_g) / R_rate_g；LP 都用 H = 4·S_max（T2 的 fr、T6 的 frRw）。
- 旧口径（T2）= (LP_fr − R_rate) / R_rate，R_rate 是 tool_lb（连续小时）。差异拆成两部分：
  机器人网格化 = R_rate_g_oldtools − R_rate；工具占用改为 freshlb 口径 = R_rate_g − R_rate_g_oldtools。
- 稳态提前大修本来能省多少（只作参照）= R_rate_g − R_rate_g^R。
- 规则：可行 = 严格证明或 60 万格复核 0 缺岗；比较只在两边都可行的实例上做。
用法：python3 code/t6_summary.py"""
import os, sys, json, glob
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names  # noqa: E402

T2J = os.path.join(ROOT, "results", "T2", "jobs")
T6J = os.path.join(ROOT, "results", "T6", "jobs")


def jl(p):
    """读一个 json 文件；不存在返回 None。"""
    return json.load(open(p)) if os.path.exists(p) else None


def t2_obj(n, kind, mult=4):
    """T2 的 LP 值：正式运行失败时取重跑（--tag）成功的那一份。"""
    for p in [os.path.join(T2J, f"{n}__{kind}__m{mult}.json")] + sorted(glob.glob(os.path.join(T2J, f"{n}__{kind}__m{mult}__*.json"))):
        d = jl(p)
        if d and d.get("obj") is not None: return d
    return None


def feas(d):
    """可行时返回资源数。"""
    if not d or d.get("R") is None or d["R"] < 0: return None
    return d["R"] if (d.get("proven") == 1 or d.get("verify_gaps") == 0) else None


def mean(xs):
    """忽略 None 的平均。"""
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def fm(x, nd=1):
    """格式化（None → —）。"""
    return "—" if x is None else f"{x:.{nd}f}"


def bounds():
    """逐实例的下界口径表。"""
    rows = []
    for n in sorted(train_names(), key=lambda n: (int(n[1:3]), n)):
        r = jl(os.path.join(T6J, f"{n}__rate.json")); w = jl(os.path.join(T6J, f"{n}__frRw__m4.json")); retry = 0
        if (w is None or w.get("obj") is None) and jl(os.path.join(T6J, f"{n}__frRw__m4__retry1.json")):   # 正式运行失败时取重跑那一份
            w2 = jl(os.path.join(T6J, f"{n}__frRw__m4__retry1.json"))
            if w2.get("obj") is not None: w, retry = w2, 1
        old = jl(os.path.join(T2J, f"{n}__rate__m0.json")); fr = t2_obj(n, "fr"); frR = t2_obj(n, "frR")
        if r is None: continue
        rec = dict(inst=n, m=int(n[1:3]), family=n.split("-")[1], has_L3_2=r["has_L3_2"],
                   R_rate=old["R_rate"], R_rate_int=old["R_rate_int"], R_rate_g=r["R_rate_g"], R_rate_g_int=r["R_rate_g_int"],
                   R_rate_g_oldtools=r["R_rate_g_oldtools"], R_rate_gR=r["R_rate_gR"], R_rate_gR_int=r["R_rate_gR_int"],
                   LP_fr=fr["obj"] if fr else None, LP_frR=frR["obj"] if frR else None,
                   LP_frRw=w["obj"] if w and w.get("obj") is not None else None)
        if w and w.get("obj") is not None:
            rec.update(frRw_viol=max(w["path_viol"], w["tool_viol"], w["wait_viol"], w["M_viol"]), frRw_gap=w["dual_gap"],
                       frRw_index=w["index_mismatch"], frRw_late=w["resets_flow_late"], frRw_sec=w["sec"], frRw_cert_sec=w["cert_sec"],
                       frRw_status=w["lp_status"], frR_sec=frR.get("sec") if frR else None, frRw_retry=retry)
        elif w: rec.update(frRw_status=w.get("lp_status"))
        g = rec["R_rate_g"]
        if rec["LP_fr"]: rec["start_fr"] = 100 * (rec["LP_fr"] - g) / g; rec["start_fr_old"] = 100 * (rec["LP_fr"] - rec["R_rate"]) / rec["R_rate"]
        if rec["LP_frRw"]: rec["start_frRw"] = 100 * (rec["LP_frRw"] - g) / g
        rec["grid_robot"] = 100 * (rec["R_rate_g_oldtools"] - rec["R_rate"]) / rec["R_rate"]
        rec["grid_tool"] = 100 * (g - rec["R_rate_g_oldtools"]) / rec["R_rate"]
        rec["steady_reset_gain"] = 100 * (g - rec["R_rate_gR"]) / g
        cyc = [c for c in r["cycles"] if c["astar"] < c["abar"] - 1e-9]
        rec["types_reset_shorter"] = sorted({c["type"] for c in cyc})
        rows.append(rec)
    return rows


def table_bounds(rows):
    """按 (m, 族) 的启动代价表，以及含 / 不含 L3-2 的分组。"""
    L = ["### 启动代价（按 (m, 族) 平均；% 相对 R_rate_g）\n",
         "| m | 族 | 实例 | R_rate（旧） | R_rate_g | LP_fr | LP_fr^Rw | 旧启动代价 % | 新：LP_fr − R_rate_g % | 新：LP_fr^Rw − R_rate_g % | 机器人网格化 % | 工具占用口径 % | 稳态提前大修收益 % |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    by = defaultdict(list)
    for r in rows: by[(r["m"], r["family"])].append(r)
    for (m, f), rs in sorted(by.items()):
        L.append(f"| {m} | {f} | {len(rs)} | {fm(mean([r['R_rate'] for r in rs]), 2)} | {fm(mean([r['R_rate_g'] for r in rs]), 2)} | "
                 f"{fm(mean([r['LP_fr'] for r in rs]), 2)} | {fm(mean([r['LP_frRw'] for r in rs]), 2)} | {fm(mean([r.get('start_fr_old') for r in rs]))} | "
                 f"{fm(mean([r.get('start_fr') for r in rs]))} | {fm(mean([r.get('start_frRw') for r in rs]))} | {fm(mean([r['grid_robot'] for r in rs]))} | "
                 f"{fm(mean([r['grid_tool'] for r in rs]))} | {fm(mean([r['steady_reset_gain'] for r in rs]), 2)} |")
    L += ["", "### 含 / 不含 L3-2 的实例分开\n",
          "| 分组 | 实例 | 新：LP_fr − R_rate_g % | 新：LP_fr^Rw − R_rate_g % | 启动期提前大修收益 (LP_fr − LP_fr^Rw)/LP_fr % | 稳态提前大修收益 % | 稳态收益 > 0 的实例 |",
          "|---|---|---|---|---|---|---|"]
    for lab, sel in (("含 L3-2", 1), ("不含 L3-2", 0)):
        rs = [r for r in rows if r["has_L3_2"] == sel]
        L.append(f"| {lab} | {len(rs)} | {fm(mean([r.get('start_fr') for r in rs]))} | {fm(mean([r.get('start_frRw') for r in rs]))} | "
                 f"{fm(mean([100 * (r['LP_fr'] - r['LP_frRw']) / r['LP_fr'] for r in rs if r['LP_fr'] and r['LP_frRw']]), 2)} | "
                 f"{fm(mean([r['steady_reset_gain'] for r in rs]), 2)} | {sum(r['steady_reset_gain'] > 1e-6 for r in rs)} |")
    return L


def rlp_rows():
    """规则重跑：fresh_noreset_u（新 / 旧补台）、fresh_rlp_w、es_lpgplan（新 / 旧）。"""
    out = []
    for n in sorted(train_names(m_list=[2, 3, 4, 6]), key=lambda n: (int(n[1:3]), n)):
        new_u = jl(os.path.join(ROOT, "results", "T6", "rlp", f"{n}__fresh_noreset_u.json"))
        old_u = jl(os.path.join(ROOT, "results", "T3b", "jobs", f"{n}__fresh_noreset_u.json"))
        rlpw = jl(os.path.join(ROOT, "results", "T6", "rlp", f"{n}__fresh_rlp_w.json"))
        old_rlpu = jl(os.path.join(ROOT, "results", "T3b", "jobs", f"{n}__fresh_rlp_u.json"))
        new_rlpu = jl(os.path.join(ROOT, "results", "T6", "rlp", f"{n}__fresh_rlp_u.json"))     # 诊断：旧 fresh_rlp_u 换新补台规则
        es_new = [jl(os.path.join(ROOT, "results", "T6", "es", f"{n}__es_lpgplan__s{s}.json")) for s in range(3)]
        es_old = [jl(os.path.join(ROOT, "results", "T3", "jobs", f"{n}__es_lpgplan__s{s}.json")) for s in range(3)]
        rec = dict(inst=n, m=int(n[1:3]), family=n.split("-")[1],
                   noreset_u_new=feas(new_u), noreset_u_old=feas(old_u), rlp_w=feas(rlpw), rlp_u_old=feas(old_rlpu), rlp_u_new=feas(new_rlpu),
                   has_noreset_u_new=int(new_u is not None), has_rlp_w=int(rlpw is not None),
                   noreset_u_new_status=new_u and new_u["status"], rlp_w_status=rlpw and rlpw["status"],
                   rlp_w_resets=rlpw and rlpw.get("resets"), rlp_w_sec=rlpw and rlpw.get("sec"), noreset_u_new_sec=new_u and new_u.get("sec"),
                   es_new=[feas(d) for d in es_new], es_old=[feas(d) for d in es_old],
                   es_new_sec=[d and d.get("sec") for d in es_new], es_old_sec=[d and d.get("sec") for d in es_old])
        out.append(rec)
    return out


def wtl(pairs):
    """(a, b) 列表 → (a 更少, 相同, a 更多, 平均 a − b)，只算两边都有值的。"""
    p = [(a, b) for a, b in pairs if a is not None and b is not None]
    return sum(a < b for a, b in p), sum(a == b for a, b in p), sum(a > b for a, b in p), mean([a - b for a, b in p]), len(p)


def table_rlp(rows):
    """规则对比表。"""
    L = ["### 规则重跑（m ≤ 6；可行解上比较，资源数之差 = 前者 − 后者）\n",
         "| m | 族 | fresh_rlp_w 对 fresh_noreset_u（新补台）：少 / 同 / 多，平均差 | noreset_u 新补台 对 旧补台 | ES lpgplan 新 对 旧（逐种子） | fresh_rlp_w 不可行 | noreset_u 新补台 不可行 |",
         "|---|---|---|---|---|---|---|"]
    by = defaultdict(list)
    for r in rows: by[(r["m"], r["family"])].append(r)
    tot = defaultdict(list)
    for (m, f), rs in sorted(by.items()):
        a = [(r["rlp_w"], r["noreset_u_new"]) for r in rs]; b = [(r["noreset_u_new"], r["noreset_u_old"]) for r in rs]
        c = [(x, y) for r in rs for x, y in zip(r["es_new"], r["es_old"])]
        tot["a"] += a; tot["b"] += b; tot["c"] += c
        fmt = lambda t: f"{t[0]} / {t[1]} / {t[2]}，{fm(t[3], 2)}"
        L.append(f"| {m} | {f} | {fmt(wtl(a))} | {fmt(wtl(b))} | {fmt(wtl(c))} | {sum(r['has_rlp_w'] and r['rlp_w'] is None for r in rs)} | {sum(r['has_noreset_u_new'] and r['noreset_u_new'] is None for r in rs)} |")
    fmt = lambda t: f"{t[0]} / {t[1]} / {t[2]}，{fm(t[3], 2)}"
    L.append(f"| 全部 | | {fmt(wtl(tot['a']))} | {fmt(wtl(tot['b']))} | {fmt(wtl(tot['c']))} | {sum(r['has_rlp_w'] and r['rlp_w'] is None for r in rows)} | {sum(r['has_noreset_u_new'] and r['noreset_u_new'] is None for r in rows)} |")
    L += ["", "### fresh_rlp_w 的改善从哪来（诊断：旧 fresh_rlp_u 只换新补台规则）\n",
          "| 比较 | 少 / 同 / 多，平均差 | 实例 |", "|---|---|---|"]
    for lab, a, b in (("fresh_rlp_u：新补台 对 旧补台（补台规则的作用）", "rlp_u_new", "rlp_u_old"),
                      ("新补台下：fresh_rlp_u 对 fresh_noreset_u", "rlp_u_new", "noreset_u_new"),
                      ("新补台下：fresh_rlp_w 对 fresh_rlp_u（窗口 LP + R_rate_g 外推的作用）", "rlp_w", "rlp_u_new"),
                      ("新补台下：fresh_rlp_w 对 fresh_noreset_u", "rlp_w", "noreset_u_new")):
        t = wtl([(r[a], r[b]) for r in rows])
        L.append(f"| {lab} | {t[0]} / {t[1]} / {t[2]}，{fm(t[3], 2)} | {t[4]} |")
    return L


def main():
    rows = bounds()
    with open(os.path.join(ROOT, "results", "T6_bounds.jsonl.tmp"), "w") as f:
        for r in rows: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(os.path.join(ROOT, "results", "T6_bounds.jsonl.tmp"), os.path.join(ROOT, "results", "T6_bounds.jsonl"))
    rr = rlp_rows()
    with open(os.path.join(ROOT, "results", "T6_rlp.jsonl.tmp"), "w") as f:
        for r in rr: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(os.path.join(ROOT, "results", "T6_rlp.jsonl.tmp"), os.path.join(ROOT, "results", "T6_rlp.jsonl"))
    w = [r for r in rows if r.get("frRw_viol") is not None]
    print(f"LP_fr^Rw：{len(w)} 个有解；最大证书违反 {max(r['frRw_viol'] for r in w):.1e}，最大 |对偶差| {max(abs(r['frRw_gap']) for r in w):.1e}，"
          f"索引不一致 {sum(r['frRw_index'] for r in w)}，窗口外提前大修流 {max(r['frRw_late'] for r in w):.1e}；"
          f"LP_fr^R ≤ LP_fr^Rw ≤ LP_fr 不成立的 {sum(not (r['LP_frR'] - 1e-3 <= r['LP_frRw'] <= r['LP_fr'] + 1e-3) for r in w)} 个\n")
    print("\n".join(table_bounds(rows)))
    print()
    print("\n".join(table_rlp(rr)))


if __name__ == "__main__":
    main()
