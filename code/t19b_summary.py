"""T19 (b) 汇总：两种默认规则加 RDD-FD（偏离窗口 0.5·S_max）。
- 接入构造 + RDD-FD：results/T19/b/<实例>__s<种子>.json（code/t19b_job.py，偏离优先）；
- fresh_rlp_w + RDD-FD：T18 的 W = 0.5·S_max（results/T18/jobs/<实例>__w0.5__s<种子>.json），设置相同（限时、ρ、q_on、(4 + 12)、T_eval、resetopt）；
- 默认规则本身（W = 0）：接入构造取 (a) 的构造台数与 rddfd 的起点定台数（base_R），fresh_rlp_w 取 T12b / T12 的 base；
- 构造用时：接入构造 = (a) 的作业墙钟（8 个 θ × 2 种入场口径，含每个构造的 60 万格判定），与 RDD-FD 限时对照。
可行 = 严格证明或 60 万格复核 0 缺岗；只在可行解上比资源数。
比较：(1) 实例 × 种子平均；(2) 同一种子配对的少 / 同 / 多；(3) 每个实例 3 个种子最好值的少 / 同 / 多。
输出 results/T19b_rdd.jsonl（逐作业）与 markdown 表（打印）。用法：python3 code/t19b_summary.py"""
import os, sys, json
from collections import defaultdict, Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names  # noqa: E402
from t18_summary import job as t18_job, feas  # noqa: E402
from t16_job import TLIM  # noqa: E402


def r_at(rec, t):
    """改进轨迹里 t 秒以前找到的最好台数（轨迹第一项是起点定台数；t 早于起点则返回 None）。"""
    xs = [x[0] for x in rec["traj"] if x[3] <= t]
    return min(xs) if xs else None


def jl(p):
    """读 json；不存在返回 None。"""
    return json.load(open(p)) if os.path.exists(p) else None


def main():
    conn = {r["inst"]: r for r in map(json.loads, open(os.path.join(ROOT, "results", "T19_connector.jsonl"))) if r["kind"] == "steady"}
    rows = []; inst = []
    for n in sorted(train_names([4, 6, 8, 10]), key=lambda n: (int(n[1:3]), n)):
        m = int(n[1:3]); fam = n.split("-")[1]
        a = json.load(open(os.path.join(ROOT, "results", "T19", "steady", f"{n}.json")))
        b0 = t18_job(n, 0.0, None)
        per = dict(conn=[], fresh=[], ded=[])
        for s in range(3):
            c = jl(os.path.join(ROOT, "results", "T19", "b", f"{n}__s{s}.json"))
            f = jl(os.path.join(ROOT, "results", "T18", "jobs", f"{n}__w0.5__s{s}.json"))
            rc = feas(c) if c else None; rf = feas(f) if f else None
            rd = r_at(c, TLIM[m] - a["wall"]) if rc is not None else None
            per["conn"].append(rc); per["fresh"].append(rf); per["ded"].append(rd)
            rows.append(dict(inst=n, m=m, family=fam, seed=s, conn_done=int(c is not None), conn_R=rc, conn_R_ded=rd, conn_base_R=c and c.get("base_R"), conn_proven=c and int(c.get("verify_cert") == 1 and c.get("exact") == 1),
                             conn_wall=c and c.get("wall"), conn_n_active=c and c.get("n_active"), conn_gens=c and c.get("gens"), conn_resets=c and c.get("resets"),
                             check_judge=c and c.get("check_judge"), check_py=c and c.get("check_py"), check_dump=c and (c.get("check_dump_same"), c.get("check_dump_n")),
                             fresh_R=rf, fresh_proven=f and int(f.get("verify_cert") == 1 and f.get("exact") == 1), fresh_n_active=f and f.get("n_active")))
        bc = [x for x in per["conn"] if x is not None]; bf = [x for x in per["fresh"] if x is not None]; bd = [x for x in per["ded"] if x is not None]
        inst.append(dict(inst=n, m=m, family=fam, conn_rule=conn[n]["conn"]["R"], fresh_rule=feas(b0), conn_best=min(bc) if bc else None, conn_best_ded=min(bd) if bd else None, fresh_best=min(bf) if bf else None,
                         n_conn=len(bc), n_fresh=len(bf), build_sec=a["wall"], tlim=TLIM[m], known_best=conn[n]["known_best"]["R"], LB=conn[n]["LB"]))
    tmp = os.path.join(ROOT, "results", "T19b_rdd.jsonl.tmp")
    with open(tmp, "w") as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=False) + "\n")
    os.replace(tmp, tmp[:-4])

    by = defaultdict(list); byi = defaultdict(list)
    for r in rows: by[(r["m"], r["family"])].append(r)
    for r in inst: byi[(r["m"], r["family"])].append(r)

    def cmp(pairs):
        """(接入, fresh) 配对 → '少 / 同 / 多'。"""
        return f"{sum(a < b for a, b in pairs)} / {sum(a == b for a, b in pairs)} / {sum(a > b for a, b in pairs)}"
    L = ["### 两种默认规则加 RDD-FD（偏离窗口 0.5·S_max）：平均资源数（可行解，实例 × 种子；括号内为可行 / 作业数）\n",
         "| m | 族 | 接入构造 + RDD-FD | fresh_rlp_w + RDD-FD | 同种子配对：接入少 / 同 / 多 | 3 个种子最好值：接入少 / 同 / 多 | 严格证明（接入 / fresh） |",
         "|---|---|---|---|---|---|---|"]
    for (m, fam), xs in sorted(by.items()):
        c = [x["conn_R"] for x in xs if x["conn_R"] is not None]; f = [x["fresh_R"] for x in xs if x["fresh_R"] is not None]
        pr = [(x["conn_R"], x["fresh_R"]) for x in xs if x["conn_R"] is not None and x["fresh_R"] is not None]
        pi = [(x["conn_best"], x["fresh_best"]) for x in byi[(m, fam)] if x["conn_best"] is not None and x["fresh_best"] is not None]
        L.append(f"| {m} | {fam} | {np.mean(c):.1f}（{len(c)}/{sum(x['conn_done'] for x in xs)}） | {np.mean(f):.1f}（{len(f)}/{len(xs)}） | {cmp(pr)} | {cmp(pi)} | "
                 f"{sum(x['conn_proven'] or 0 for x in xs)} / {sum(x['fresh_proven'] or 0 for x in xs)} |")
    pr = [(x["conn_R"], x["fresh_R"]) for x in rows if x["conn_R"] is not None and x["fresh_R"] is not None]
    pi = [(x["conn_best"], x["fresh_best"]) for x in inst if x["conn_best"] is not None and x["fresh_best"] is not None]
    L.append(f"\n合计：同种子配对 {len(pr)} 对，接入少 / 同 / 多 = {cmp(pr)}，平均差（接入 − fresh）{np.mean([a - b for a, b in pr]):+.2f}，差的分布 {dict(sorted(Counter(a - b for a, b in pr).items()))}；"
             f"3 个种子最好值 {len(pi)} 个实例，{cmp(pi)}，平均差 {np.mean([a - b for a, b in pi]):+.2f}，分布 {dict(sorted(Counter(a - b for a, b in pi).items()))}。")
    L.append(f"可行：接入构造 + RDD-FD {sum(x['conn_R'] is not None for x in rows)} / {sum(x['conn_done'] for x in rows)}；fresh_rlp_w + RDD-FD {sum(x['fresh_R'] is not None for x in rows)} / {len(rows)}。")
    L += ["", "### 默认规则本身、加 RDD-FD 后与已知最好（按 (m, 族) 平均；最好值取 3 个种子）\n",
          "| m | 族 | 下界 | 接入构造（规则本身） | rddfd 起点定台数 | 接入 + RDD-FD（最好种子） | fresh_rlp_w（规则本身） | fresh + RDD-FD（最好种子） | 已知最好 RDD 解 |", "|---|---|---|---|---|---|---|---|---|"]
    for (m, fam), xs in sorted(byi.items()):
        rs = by[(m, fam)]
        L.append(f"| {m} | {fam} | {np.mean([x['LB'] for x in xs]):.1f} | {np.mean([x['conn_rule'] for x in xs]):.1f} | {np.mean([r['conn_base_R'] for r in rs if r['conn_base_R'] is not None and r['conn_base_R'] >= 0]):.1f} | "
                 f"{np.mean([x['conn_best'] for x in xs if x['conn_best'] is not None]):.1f} | {np.mean([x['fresh_rule'] for x in xs if x['fresh_rule'] is not None]):.1f} | "
                 f"{np.mean([x['fresh_best'] for x in xs if x['fresh_best'] is not None]):.1f} | {np.mean([x['known_best'] for x in xs]):.1f} |")
    L += ["", "### 构造用时（接入构造：8 个 θ × 2 种入场口径，含每个构造的 60 万格判定；秒）\n",
          "| m | 族 | 中位 | 最大 | RDD-FD 限时 | 最大占限时 |", "|---|---|---|---|---|---|"]
    for (m, fam), xs in sorted(byi.items()):
        bsec = [x["build_sec"] for x in xs]
        L.append(f"| {m} | {fam} | {np.median(bsec):.1f} | {max(bsec):.1f} | {xs[0]['tlim']} | {100 * max(bsec) / xs[0]['tlim']:.1f}% |")
    pd = [(x["conn_R_ded"], x["fresh_R"]) for x in rows if x["conn_R_ded"] is not None and x["fresh_R"] is not None]
    pid = [(x["conn_best_ded"], x["fresh_best"]) for x in inst if x["conn_best_ded"] is not None and x["fresh_best"] is not None]
    L.append(f"\n扣构造用时（接入一侧取「限时 − 构造用时」时刻的台数）：同种子配对 {len(pd)} 对，{cmp(pd)}，平均差 {np.mean([a - b for a, b in pd]):+.2f}；"
             f"3 个种子最好值 {len(pid)} 个实例，{cmp(pid)}，平均差 {np.mean([a - b for a, b in pid]):+.2f}；扣后台数变多的作业 {sum(x['conn_R_ded'] != x['conn_R'] for x in rows if x['conn_R_ded'] is not None)} 个。")
    cj = [r["check_judge"] for r in rows if r["check_judge"] is not None]; cp = [r["check_py"] for r in rows if r["check_py"] is not None]
    cd = [r["check_dump"] for r in rows if r["check_dump"]]
    L.append(f"\n核对（接入构造 + RDD-FD 的作业）：最终解重新判定逐字段相同 {sum(cj)} / {len(cj)}；Python 重跑最终解（T_eval 格）与 eval 相同 {sum(cp)} / {len(cp)}；"
             f"抽样缺岗与 eval 相同 {sum(a for a, b in cd)} / {sum(b for a, b in cd)}。完成 {sum(r['conn_done'] for r in rows)} / {len(rows)}。")
    L.append(f"激活偏离数（平均）：接入 {np.mean([r['conn_n_active'] for r in rows if r['conn_n_active'] is not None]):.1f}，fresh {np.mean([r['fresh_n_active'] for r in rows if r['fresh_n_active'] is not None]):.1f}；"
             f"代数（接入，平均）{np.mean([r['conn_gens'] for r in rows if r['conn_gens'] is not None]):,.0f}。")
    L += ["", "### 逐实例（3 个种子最好值）\n", "| 实例 | 下界 | 接入构造 | 接入 + RDD-FD | fresh_rlp_w | fresh + RDD-FD | 已知最好 |", "|---|---|---|---|---|---|---|"]
    for x in inst:
        L.append(f"| {x['inst']} | {x['LB']} | {x['conn_rule']} | {x['conn_best'] if x['conn_best'] is not None else '—'}（{x['n_conn']}/3） | "
                 f"{x['fresh_rule'] if x['fresh_rule'] is not None else '不可行'} | {x['fresh_best'] if x['fresh_best'] is not None else '—'}（{x['n_fresh']}/3） | {x['known_best']} |")
    print("\n".join(L))


if __name__ == "__main__":
    main()
