"""T20-1 下界与时间界的 LP（实验 1）：一个 (测试实例, 种类, H 倍数)。第五批第一次使用测试实例（11-06 用户定）。
种类：
- rate：速率下界 R_rate_g 等（T6 的 t6_job.rate_job：网格口径，LP 与整数，不提前大修 / 允许稳态提前大修）；倍数记 0；
- fr：LP_fr(H)，全新启动、不提前大修（T2 的 t2_job.fresh_job：求解、索引核对、对偶证书、保存对偶与值岗流）；
- frRw：LP_fr^Rw(H)，提前大修窗口固定为 [0, 2·S_max)（T6 的 t6_job.frRw_job；H ≤ 2·S_max 时窗口就是整个时域）。
H = ⌊倍数 · S_max⌋，倍数 ∈ {0.5, 1, 1.5, 2, 3, 4}；单个 LP 限时 1 小时（TASKS T20-1）。
失败（数值失败或超时）时按 T2 / T6 / T17 的做法只换求解设置重跑一次：加参数 retry，行列随机重排（种子 2）+ 关 crossover，
结果写到 <名>__retry1.json，原结果不动；仍失败的在报告里列出。
结果写进 results/T20/lb/jobs/<实例>__<种类>__m<倍数>.json（原子写入，已存在就跳过），对偶写进 results/T20/lb_duals/。
用法：python3 code/t20_lb_job.py <实例> <rate|fr|frRw> <倍数> [retry]"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_inst  # noqa: E402
from freshlb import S_max  # noqa: E402
import t2_job  # noqa: E402
import t6_job  # noqa: E402

OUT = os.path.join(ROOT, "results", "T20", "lb", "jobs")
DUAL = os.path.join(ROOT, "results", "T20", "lb_duals")
TL = 3600


def main():
    """解析参数，调用 T2 / T6 的求解函数，原子写入结果。"""
    inst, kind, mult = sys.argv[1], sys.argv[2], float(sys.argv[3])
    retry = len(sys.argv) > 4 and sys.argv[4] == "retry"
    perm, cross = (2, "off") if retry else (None, "on")
    assert kind in ("rate", "fr", "frRw")
    os.makedirs(OUT, exist_ok=True); os.makedirs(DUAL, exist_ok=True)
    name = f"{inst}__{kind}__m{mult:g}" + ("__retry1" if retry else ""); path = os.path.join(OUT, name + ".json")
    if os.path.exists(path): print("已完成，跳过", name); return
    I = load_inst(inst, "test"); S = S_max(I); H = int(mult * S); W = 2 * S
    rec = dict(inst=inst, m=I.meta["m"], family=I.meta["family"], kind=kind, mult=mult, S=S, H=H, W=min(W, H), tl=TL)
    t0 = time.time(); dp = os.path.join(DUAL, name + ".npz")
    if kind == "rate": rec.update(t6_job.rate_job(I))
    elif kind == "fr": rec.update(t2_job.fresh_job(I, H, False, dp, perm, cross, tl=TL))
    else: rec.update(t6_job.frRw_job(I, H, W, dp, perm, cross, tl=TL))
    if retry: rec.update(tag="retry1", perm=perm, crossover=cross)
    rec["wall"] = round(time.time() - t0, 1)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: v for k, v in rec.items() if k not in ("path_viol_k", "cycles", "y_g", "y_gR", "z_gR")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
