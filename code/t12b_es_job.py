"""T12b 单个作业：严格窗口下重跑 ES lpgplan（与"严格窗口的 fresh_rlp_w"同一场景比较）。
T6 的 ES lpgplan 用 T3 的 fresh_rlp.planf 做提前大修计划，但没有限制执行时刻，重放发现大多数解在 2·S_max 之后还有提前大修，
不属于"只在启动窗口内提前大修"的场景。这里参数与 T6 的 es_lpgplan 完全相同（certs、365 格、10 代、(4 + 12)、60000 / 60 万格），
只加 planW = 2·S_max（rost5 es 的 T12b 选项）。
作业内核对：最终解用 rost5 eval（同一组选项 + planW）重新判定，60 万格的缺岗、判环结果、周期与 es 的输出相同；记下最后一次提前大修时刻。
结果写进 results/T12b/es/<实例>__es_lpgplan_w__s<种子>.json（原子写入，已存在就跳过；运行时加 fcntl 锁）。
用法：python3 code/t12b_es_job.py <实例> <种子>"""
import os, sys, json, time, subprocess, fcntl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")
OUT = os.path.join(ROOT, "results", "T12b", "es")
T3R = os.path.join(ROOT, "results", "T3", "rlp")


def main():
    """跑 es（加 planW），重新判定最终解，原子写入。"""
    name, seed = sys.argv[1], int(sys.argv[2])
    os.makedirs(OUT, exist_ok=True)
    jn = f"{name}__es_lpgplan_w__s{seed}"; path = os.path.join(OUT, jn + ".json")
    if os.path.exists(path): print("已完成，跳过", jn); return
    lk = open(path + ".lock", "w")
    try: fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError: print("另一个进程正在做", jn); sys.exit(4)
    if os.path.exists(path): print("已完成，跳过", jn); return
    I = load_train(name); _, txt = inst_paths(name, "train"); S = S_max(I); W = 2 * S
    ro = [f"lpg={T3R}/{name}.lpg.txt", f"planf={T3R}/{name}__fresh_rlp.planf"]
    args = [R5, "es", txt, "certs", "365", "10", "4", "12", str(seed), "60000", "600000"] + ro + [f"planW={W}"]
    t0 = time.time()
    r = subprocess.run(args, capture_output=True, text=True)
    try: rec = json.loads(r.stdout)
    except Exception: print("rost5 输出无法解析", r.stdout[-300:], r.stderr[-300:]); sys.exit(2)
    rec.update(inst=name, m=I.meta["m"], family=I.meta["family"], kind="es_lpgplan_w", seed=seed, W=W, wall=round(time.time() - t0, 1), cmd=" ".join(args))
    rec["proven"] = int(rec.get("verify_cert") == 1 and rec.get("exact") == 1)
    if rec["R"] >= 0:
        ev = [str(x) for x in rec["Nk"]] + ["|"] + [str(x) for x in rec["Mj"]] + ["|"] + [str(x) for x in rec["w"]] + ["|"] + ro + ["curdef", f"planW={W}"]
        v = json.loads(subprocess.run([R5, "eval", txt, "600000", "1"] + ev, capture_output=True, text=True, check=True).stdout)
        rec["check_eval"] = int(v["gaps"] == rec["verify_gaps"] and v["cert"] == rec["verify_cert"] and (rec["verify_cert"] != 1 or v["period"] == rec["verify_period"]))
        rec["last_reset"] = v["last_reset"]; rec["window_ok"] = int(v["last_reset"] < W)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: rec.get(k) for k in ("inst", "seed", "R", "status", "proven", "last_reset", "W", "window_ok", "check_eval", "wall")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
