"""T16 单个作业：对偶引导偏离位置。一个 (训练实例, β, 种子) 的 RDD-FD（rost5 rddfd），其余设置与 T12 相同：
fresh_rlp_w（严格窗口 planW、后备补足时域 hsfb）+ resetopt，ρ = 0.1、q_on = 0.4、(4 + 12)、T_eval = 6·S_max、互换 0.5、stall 30。
位置分布：β > 0 时，决策 i（站 s、值岗覆盖 [t, t + d)）的权重 w_i = Σ_{u∈[t, t+d)} max(π_{s,u}, 0)，π 取 LP_fr^Rw(4S) 的覆盖对偶
（results/T6_duals/<实例>__frRw__m4[__retry1].npz）；q_i = (1 − β)/N_dev + β·w_i/Σw；槽的位置基因 p 映射到第一个累积 q > p 的决策。
β = 0 时不传对偶文件，与原 RDD 逐位相同。
限时（墙钟）：m = 8、10 按 T12 的收敛曲线定（results/T16/tlim.json，没有时为 2 小时）；m = 4、6 按 T11 的实测定：T11 跑 1000 代在 m = 4 最长 50 秒、m = 6 最长 128 秒，取 300 秒、600 秒。
作业内核对（同 T12）：最终解用 rost5 judge（带同一个位置分布 cdf=）重新判定逐字段相同；抽样 T_eval 缺岗与 rost5 eval 相同。
结果写进 results/T16/jobs/<实例>__b<β>__s<种子>.json（原子写入，已存在就跳过；运行时加 fcntl 锁）。
用法：python3 code/t16_job.py <实例> <β> <种子>"""
import os, sys, json, time, subprocess, fcntl, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
from t11_job import write_geno  # noqa: E402
from t12_job import rule_opts  # noqa: E402
from t6_rerun import frrw_npz  # noqa: E402

R5 = os.environ.get("ROST5_BIN", os.path.join(ROOT, "code", "cpp", "rost5"))
OUT = os.environ.get("T16_OUT", os.path.join(ROOT, "results", "T16", "jobs"))
DW = os.path.join(ROOT, "results", "T16", "dualw")
TLIM = {4: 300, 6: 600, 8: 7200, 10: 7200}
# m = 8、10 的限时按 T12 的收敛曲线定（10-30 晚审核），定下后写进 results/T16/tlim.json（{"8": 秒, "10": 秒}），作业启动时读
_tl = os.path.join(ROOT, "results", "T16", "tlim.json")
if os.path.exists(_tl): TLIM.update({int(k): float(v) for k, v in json.load(open(_tl)).items()})
if os.environ.get("T16_TLIM"): TLIM = {k: float(os.environ["T16_TLIM"]) for k in TLIM}   # 只用于开发期试算


def dual_file(name, W):
    """覆盖对偶 π_{s,t}（t < W，非零的）写成 "s t π"；没有就生成。"""
    os.makedirs(DW, exist_ok=True); p = os.path.join(DW, f"{name}.txt")
    if not os.path.exists(p):
        pi = np.load(frrw_npz(name))["pi"]
        with open(p + f".{os.getpid()}.tmp", "w") as f:
            for s in range(pi.shape[0]):
                for t in range(min(W + 64, pi.shape[1])):     # 值岗可以跨过窗口末端，多写 64 格
                    if pi[s, t] != 0: f.write(f"{s} {t} {float(pi[s, t])!r}\n")
        os.replace(p + f".{os.getpid()}.tmp", p)
    return p


def main():
    name, beta, seed = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
    os.makedirs(OUT, exist_ok=True)
    jn = f"{name}__b{beta:g}__s{seed}"; path = os.path.join(OUT, jn + ".json")
    if os.path.exists(path): print("已完成，跳过", jn); return
    lk = open(path + ".lock", "w")
    try: fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError: print("另一个进程正在做", jn); sys.exit(4)
    if os.path.exists(path): print("已完成，跳过", jn); return
    I = load_train(name); _, txt = inst_paths(name, "train"); S = S_max(I); W = 2 * S; m = I.meta["m"]
    opts = rule_opts(name, S) + ([f"dualw={dual_file(name, W)}", f"beta={beta!r}"] if beta > 0 else [])
    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix=f"t16_{jn}_") as tmpd:
        dump = os.path.join(tmpd, "dump.tsv")
        args = [R5, "rddfd", txt, "365", "60000", "600000", str(6 * S), str(10 ** 9), str(TLIM[m]), "4", "12", str(seed), "0.1", "0.4"] + opts + \
               [f"W={W}", "resetopt", "stall=30", "swap=0.5", "fdsel=best", f"fddump={dump}"]
        r = subprocess.run(args, capture_output=True, text=True)
        try: rec = json.loads(r.stdout)
        except Exception: print("rost5 输出无法解析", r.stdout[-300:], r.stderr[-300:]); sys.exit(2)
        rec.update(inst=name, m=m, family=I.meta["family"], beta=beta, S=S, wall=round(time.time() - t0, 1), cmd=" ".join(args[:-1]))
        # 核对：judge 与 eval 都用同一个位置分布
        cdf = []
        if rec.get("cdf"):
            cf = os.path.join(tmpd, "cdf.txt"); open(cf, "w").write("\n".join(repr(x) for x in rec["cdf"]) + "\n"); cdf = [f"cdf={cf}"]
        ev_opts = [o for o in opts if not o.startswith(("hsfb=", "dualw=", "beta="))]
        if rec["R"] >= 0:
            gp = os.path.join(tmpd, "best.geno"); write_geno(gp, rec["W"], rec["ndev"], rec["qon"], rec["genes"], 1)
            j = json.loads(subprocess.run([R5, "judge", txt, "60000", "600000", ",".join(map(str, rec["Nk"])), ",".join(map(str, rec["Mj"]))]
                                          + ev_opts + cdf + [f"W={rec['W']}", "resetopt", f"rdd={gp}"], capture_output=True, text=True, check=True).stdout)
            rec["check_judge"] = int(j["R"] == rec["R"] and j["feasible"] == 1 and j["status"] == rec["status"] and j["verify_gaps"] == rec["verify_gaps"]
                                     and j["verify_cert"] == rec["verify_cert"] and j["verify_period"] == rec["verify_period"] and j["exact"] == rec["exact"]
                                     and j["resets"] == rec["resets"] and j["last_reset"] == rec["last_reset"])
        n_ok = n_all = 0
        for line in open(dump):
            Nk, Mj, g, genes = line.rstrip("\n").split("\t")
            sp = os.path.join(tmpd, "s.geno"); write_geno(sp, rec["W"], rec["ndev"], rec["qon"], [float(x) for x in genes.split(",")], 1)
            e = json.loads(subprocess.run([R5, "eval", txt, str(rec["Teval"]), "0"] + [str(x) for x in json.loads(Nk)] + ["|"] + [str(x) for x in json.loads(Mj)]
                                          + ["|", "|"] + ev_opts + cdf + [f"rdd={sp}"], capture_output=True, text=True, check=True).stdout)
            n_all += 1; n_ok += int(e["gaps"] == int(g))
        rec.update(check_dump_n=n_all, check_dump_same=n_ok)
    rec["proven"] = int(rec.get("verify_cert") == 1 and rec.get("exact") == 1)
    rec["feasible"] = int(rec.get("R", -1) >= 0 and (rec["proven"] == 1 or rec.get("verify_gaps") == 0))
    rec["window_ok"] = int(rec.get("last_reset", -1) < W)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: rec.get(k) for k in ("inst", "beta", "seed", "base_R", "R", "status", "feasible", "wall", "check_judge", "check_dump_same")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
