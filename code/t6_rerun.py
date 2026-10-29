"""T6 的规则重跑（m ≤ 6 的训练实例），结果写进 results/T6/rlp/ 与 results/T6/es/（原子写入，已存在就跳过）。
- noreset_u：用新补台规则（rost5 的新默认）重跑 T3b 的 fresh_noreset_u，规则文件不变（results/T3/rlp/）；
- rlp_w：新变体 fresh_rlp_w。规则由 make_w 生成：
    * 值岗计划 F、相位计划 FP：LP_fr^Rw(4·S_max) 的值岗流（T6 的 frRw 作业，results/T6_duals/），不缩放（noscale）；
    * 提前大修计划：LP_fr^Rw 中 D < W = 2·S_max 的全部提前大修流（取满启动窗口），按累计流四舍五入成整数次；
    * t ≥ H 后的外推比例：网格口径速率 LP R_rate_g（不提前大修）的班率 y[k,s]；
    * 定台数起点：LP_fr^Rw 的台数向上取整。
- es_lpgplan：用新补台规则重跑 T3 的 ES lpgplan（10 代，种子 0–2），参数与 T3 相同。
缺依赖（T6 的 frRw 或 rate 结果）时退出码 3、不写结果。作业运行时对结果文件加 fcntl 锁，拿不到锁以返回码 4 退出。
- rlp_u（诊断，不在任务书里）：T3b 的 fresh_rlp_u 只换新补台规则重跑，用来拆分 fresh_rlp_w 的改善来自补台规则还是来自窗口 LP 与新外推。
用法：python3 code/t6_rerun.py noreset_u <实例> | rlp_w <实例> | rlp_u <实例> | es_lpgplan <实例> <种子>"""
import os, sys, json, time, fcntl, shutil, subprocess, tempfile, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
import t3_make_rlp as M  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")
OUT = os.path.join(ROOT, "results", "T6")
RULES = os.path.join(OUT, "rules")      # fresh_rlp_w 的规则文件（可由本脚本重新生成，不进仓库）
T6J = os.path.join(OUT, "jobs"); T6D = os.path.join(ROOT, "results", "T6_duals")


def frrw_npz(name):
    """LP_fr^Rw(4S) 的对偶与流量文件；正式运行数值失败、由重跑（__retry1）解出的实例取重跑那一份。"""
    p = os.path.join(T6D, f"{name}__frRw__m4.npz")
    q = os.path.join(T6D, f"{name}__frRw__m4__retry1.npz")
    return p if os.path.exists(p) or not os.path.exists(q) else q


def make_w(name, out):
    """生成 fresh_rlp_w 的规则文件与 planf（见文件头）。返回诊断信息。"""
    I = load_train(name)
    rate = json.load(open(os.path.join(T6J, f"{name}__rate.json")))
    y = {}
    for key, v in rate["y_g"].items():
        k, s = map(int, key.split(",")); y.setdefault(s, {})[k] = v
    d = np.load(frrw_npz(name))
    H, W = int(d["H"]), int(d["W"]); N, Mv, fl = d["N"], d["M"], d["duty_flow"]
    ceil = lambda v: int(math.ceil(v - 1e-6))
    info = M.write_rlp(os.path.join(out, f"{name}__fresh_rlp_w.rlp"), I, H, False, N, [ceil(v) for v in N], [ceil(v) for v in Mv], fl, y, noscale=True)
    info["planf_entries"] = M.write_planf(os.path.join(out, f"{name}__fresh_rlp_w.planf"), H, N, fl, cut=W)
    info.update(H=H, W=W, reset_flow=float(sum(x[5] for x in fl if int(x[4]) == 1)))
    return info


def ensure_w(name):
    """加锁生成规则文件（先写临时目录再改名）；已有完成标记就跳过。"""
    os.makedirs(RULES, exist_ok=True)
    mark = os.path.join(RULES, f"{name}.w.done")
    if os.path.exists(mark): return
    with open(os.path.join(RULES, f"{name}.lock"), "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        if os.path.exists(mark): return
        tmp = tempfile.mkdtemp(dir=RULES, prefix=f".{name}.")
        info = make_w(name, tmp)
        for f in os.listdir(tmp): os.replace(os.path.join(tmp, f), os.path.join(RULES, f))
        shutil.rmtree(tmp)
        with open(mark, "w") as f: f.write(json.dumps(info) + "\n")


def main():
    kind, name = sys.argv[1], sys.argv[2]
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else None
    I = load_train(name); _, txt = inst_paths(name, "train")
    if kind == "noreset_u":
        d = os.path.join(OUT, "rlp"); jn = f"{name}__fresh_noreset_u"
        args = [R5, "rlp", txt, os.path.join(M.OUT, f"{name}__fresh_noreset_u.rlp"), "365", "60000", "600000"]
    elif kind == "rlp_u":      # 诊断：旧的 fresh_rlp_u（LP_fr^R、旧速率外推）只换新补台规则，用来拆分 fresh_rlp_w 的改善来源
        d = os.path.join(OUT, "rlp"); jn = f"{name}__fresh_rlp_u"
        args = [R5, "rlp", txt, os.path.join(M.OUT, f"{name}__fresh_rlp_u.rlp"), "365", "60000", "600000",
                "planf=" + os.path.join(M.OUT, f"{name}__fresh_rlp_u.planf")]
    elif kind == "rlp_w":
        d = os.path.join(OUT, "rlp"); jn = f"{name}__fresh_rlp_w"
        if not (os.path.exists(frrw_npz(name)) and os.path.exists(os.path.join(T6J, f"{name}__rate.json"))):
            print("T6 依赖未齐", name); sys.exit(3)
        ensure_w(name)
        args = [R5, "rlp", txt, os.path.join(RULES, f"{name}__fresh_rlp_w.rlp"), "365", "60000", "600000",
                "planf=" + os.path.join(RULES, f"{name}__fresh_rlp_w.planf")]
    elif kind == "es_lpgplan":
        d = os.path.join(OUT, "es"); jn = f"{name}__es_lpgplan__s{seed}"
        args = [R5, "es", txt, "certs", "365", "10", "4", "12", str(seed), "60000", "600000",
                "lpg=" + os.path.join(M.OUT, f"{name}.lpg.txt"), "planf=" + os.path.join(M.OUT, f"{name}__fresh_rlp.planf")]
    else:
        sys.exit(f"未知种类 {kind}")
    os.makedirs(d, exist_ok=True); path = os.path.join(d, jn + ".json")
    if os.path.exists(path): print("已完成，跳过", jn); return
    lk = open(path + ".lock", "w")      # 两个调度器从作业表两头跑时，同一作业只让一个进程做（单写入者）
    try: fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError: print("另一个进程正在做", jn); sys.exit(4)
    if os.path.exists(path): print("已完成，跳过", jn); return
    t0 = time.time()
    r = subprocess.run(args, capture_output=True, text=True)
    try: rec = json.loads(r.stdout)
    except Exception: print("rost5 输出无法解析", r.stdout[-300:], r.stderr[-300:]); sys.exit(2)
    rec.update(inst=name, m=I.meta["m"], family=I.meta["family"], kind=jn.split("__", 1)[1].split("__s")[0], seed=seed,
               wall=round(time.time() - t0, 1), cmd=" ".join(args))
    rec["proven"] = int(rec.get("verify_cert") == 1 and rec.get("exact") == 1)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: rec.get(k) for k in ("inst", "kind", "seed", "R", "status", "verify_gaps", "wall")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
