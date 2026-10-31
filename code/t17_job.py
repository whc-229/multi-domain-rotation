"""T17 单个作业：LP 下界随时域 H 的收敛（时间界的松弛一侧）。一个 (训练实例, 种类, H 倍数)。
种类：
- fr：LP_fr(H)，全新启动、不提前大修；求解、索引核对、对偶证书、保存对偶与值岗流都调用 T2 的 t2_job.fresh_job（不改代码）；
- frRw：LP_fr^Rw(H)，提前大修窗口固定为 [0, 2·S_max)；调用 T6 的 t6_job.frRw_job（W = 2·S_max）。
  H < 2·S_max 时没有要删的列，窗口自然就是 [0, H)。
H = ⌊倍数 · S_max⌋（与 T2、T6 的取整相同）。倍数 2、4 不在这里算：fr 复用 T2 的 fr__m2 / fr__m4；
frRw 的 4S 复用 T6 的 frRw__m4，2S 复用 T2 的 frR__m2（H = 2·S_max 时窗口覆盖整个时域，两者是同一个 LP）。
结果写进 results/T17/jobs/<实例>__<种类>__m<倍数>.json（原子写入，已存在就跳过），对偶写进 results/T17_duals/。
用法：python3 code/t17_job.py <实例> <fr|frRw> <倍数>"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train  # noqa: E402
from freshlb import S_max  # noqa: E402
import t2_job  # noqa: E402
import t6_job  # noqa: E402

OUT = os.environ.get("T17_OUT", os.path.join(ROOT, "results", "T17", "jobs"))       # 环境变量只用于开发期试算
DUAL = os.environ.get("T17_DUAL", os.path.join(ROOT, "results", "T17_duals"))


def main():
    """解析参数，调用 T2 / T6 的求解函数，原子写入结果。"""
    inst, kind, mult = sys.argv[1], sys.argv[2], float(sys.argv[3])
    assert kind in ("fr", "frRw") and mult not in (2.0, 4.0)
    os.makedirs(OUT, exist_ok=True); os.makedirs(DUAL, exist_ok=True)
    name = f"{inst}__{kind}__m{mult:g}"; path = os.path.join(OUT, name + ".json")
    if os.path.exists(path): print("已完成，跳过", name); return
    I = load_train(inst); S = S_max(I); H = int(mult * S); W = 2 * S
    rec = dict(inst=inst, m=I.meta["m"], family=I.meta["family"], kind=kind, mult=mult, S=S, H=H, W=min(W, H))
    t0 = time.time()
    dp = os.path.join(DUAL, name + ".npz")
    rec.update(t2_job.fresh_job(I, H, False, dp) if kind == "fr" else t6_job.frRw_job(I, H, W, dp))
    rec["wall"] = round(time.time() - t0, 1)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: v for k, v in rec.items() if k not in ("path_viol_k",)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
