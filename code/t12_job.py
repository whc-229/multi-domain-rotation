"""T12 单个作业：大实例（m ∈ {8, 10, 15, 20}）上的"松弛规则 + RDD"，场景 = 全新启动、只在启动窗口内提前大修。
结果写进 results/T12/jobs/<名>.json（原子写入，已存在就跳过；运行时对结果文件加 fcntl 锁，拿不到锁以返回码 4 退出）。

默认规则：fresh_rlp_w（T6 的规则文件与 planf，results/T6/rules/），加两个选项：
  - planW=W：计划提前大修只在 t < W 时执行。原 fresh_rlp_w 的计划按 LP 的值岗开始时刻 D < W 给出，
    但派工推迟时个别提前大修会落到 W 之后几格；严格的场景定义要求都在窗口内，所以 T12 一律加这个选项；
  - hsfb=2H：365 格补足时域下 4 个起点都定不出台数时，用 2H = 8·S_max 再跑一遍（T10 的修复）。
两类作业：
  - base <实例>：默认规则本身定台数 + 60 万格复核（rost5 size），作为"fresh_rlp_w（严格窗口）"对照；
  - rdd <实例> <种子>：RDD-FD（rost5 rddfd），备选表加"本班后提前大修"（resetopt，只在窗口内、相位 ≠ L−1 时可用），
    ρ = 0.1、q_on = 0.4、(4 + 12)、T_eval = 6·S_max、互换概率 0.5、stall = 30（与 T11 的 main 设置相同），预算按墙钟时间（见 TLIM）。
rdd 作业内核对（同 T11）：最终解用 rost5 judge 重新判定，逐字段相同；抽样的 T_eval 缺岗与 rost5 eval 相同。
另外检查每个解的最后一次提前大修时刻 < W（last_reset）。
用法：python3 code/t12_job.py base <实例> | rdd <实例> <种子>"""
import os, sys, json, time, subprocess, fcntl, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
from t11_job import write_geno  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")
OUT = os.environ.get("T12_OUT", os.path.join(ROOT, "results", "T12", "jobs"))   # 环境变量只用于开发期试算
RULES = os.path.join(ROOT, "results", "T6", "rules")
GENS = 10 ** 9                                      # 代数不设上限，预算只看墙钟时间
TLIM = {8: 7200, 10: 7200, 15: 14400, 20: 14400}    # 每个作业的搜索限时（秒），见报告
if os.environ.get("T12_TLIM"): TLIM = {k: float(os.environ["T12_TLIM"]) for k in TLIM}   # 只用于开发期试算


def rule_opts(name, S):
    """默认规则 fresh_rlp_w（严格窗口 + 后备补足时域）的 rost5 选项。"""
    return [f"rlp={RULES}/{name}__fresh_rlp_w.rlp", f"planf={RULES}/{name}__fresh_rlp_w.planf", f"planW={2 * S}", f"hsfb={8 * S}"]


def eval_opts(opts):
    """rost5 eval 第 4 段能接受的选项（去掉只对定台数有用的 hsfb）。"""
    return [o for o in opts if not o.startswith("hsfb=")]


def checks(rec, txt, opts, dump, tmpd):
    """作业内核对：judge 重判最终解；抽样缺岗与 eval 相同。"""
    out = {}
    gp = os.path.join(tmpd, "best.geno"); write_geno(gp, rec["W"], rec["ndev"], rec["qon"], rec["genes"], 1)
    if rec["R"] >= 0:
        j = json.loads(subprocess.run([R5, "judge", txt, "60000", "600000", ",".join(map(str, rec["Nk"])), ",".join(map(str, rec["Mj"]))]
                                      + opts + [f"W={rec['W']}", "resetopt", f"rdd={gp}"], capture_output=True, text=True, check=True).stdout)
        same = (j["R"] == rec["R"] and j["feasible"] == 1 and j["status"] == rec["status"] and j["verify_gaps"] == rec["verify_gaps"]
                and j["verify_cert"] == rec["verify_cert"] and j["verify_period"] == rec["verify_period"] and j["exact"] == rec["exact"]
                and j["resets"] == rec["resets"] and j["last_reset"] == rec["last_reset"])
        out.update(check_judge=int(same), judge=j)
    n_ok = n_all = 0
    for line in open(dump):
        Nk, Mj, g, genes = line.rstrip("\n").split("\t")
        sp = os.path.join(tmpd, "s.geno"); write_geno(sp, rec["W"], rec["ndev"], rec["qon"], [float(x) for x in genes.split(",")], 1)
        args = [R5, "eval", txt, str(rec["Teval"]), "0"] + [str(x) for x in json.loads(Nk)] + ["|"] + [str(x) for x in json.loads(Mj)] + ["|", "|"] \
            + eval_opts(opts) + [f"rdd={sp}"]
        e = json.loads(subprocess.run(args, capture_output=True, text=True, check=True).stdout)
        n_all += 1; n_ok += int(e["gaps"] == int(g))
    out.update(check_dump_n=n_all, check_dump_same=n_ok)
    return out


def main():
    kind, name = sys.argv[1], sys.argv[2]
    seed = int(sys.argv[3]) if kind == "rdd" else None
    os.makedirs(OUT, exist_ok=True)
    jn = f"{name}__base" if kind == "base" else f"{name}__rdd__s{seed}"
    path = os.path.join(OUT, jn + ".json")
    if os.path.exists(path): print("已完成，跳过", jn); return
    lk = open(path + ".lock", "w")
    try: fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError: print("另一个进程正在做", jn); sys.exit(4)
    if os.path.exists(path): print("已完成，跳过", jn); return
    I = load_train(name); _, txt = inst_paths(name, "train"); S = S_max(I); W = 2 * S; m = I.meta["m"]
    opts = rule_opts(name, S)
    t0 = time.time()
    if kind == "base":
        args = [R5, "size", txt, "365", "60000", "600000"] + opts
        rec = json.loads(subprocess.run(args, capture_output=True, text=True, check=True).stdout)
        rec.update(inst=name, m=m, family=I.meta["family"], kind=kind, S=S, W=W, wall=round(time.time() - t0, 1), cmd=" ".join(args))
    else:
        with tempfile.TemporaryDirectory(prefix=f"t12_{jn}_") as tmpd:
            dump = os.path.join(tmpd, "dump.tsv")
            args = [R5, "rddfd", txt, "365", "60000", "600000", str(6 * S), str(GENS), str(TLIM[m]), "4", "12", str(seed), "0.1", "0.4"] + opts + \
                   [f"W={W}", "resetopt", "stall=30", "swap=0.5", "fdsel=best", f"fddump={dump}"]
            r = subprocess.run(args, capture_output=True, text=True)
            try: rec = json.loads(r.stdout)
            except Exception: print("rost5 输出无法解析", r.stdout[-300:], r.stderr[-300:]); sys.exit(2)
            rec.update(inst=name, m=m, family=I.meta["family"], kind=kind, S=S, wall=round(time.time() - t0, 1), cmd=" ".join(args[:-1]))
            t1 = time.time(); rec.update(checks(rec, txt, opts, dump, tmpd)); rec["check_sec"] = round(time.time() - t1, 1)
        rec["base_proven"] = int(rec.get("base_verify_cert") == 1 and rec.get("base_exact") == 1)
    rec["proven"] = int(rec.get("verify_cert") == 1 and rec.get("exact") == 1)
    rec["window_ok"] = int(rec.get("last_reset", -1) < W)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: rec.get(k) for k in ("inst", "kind", "seed", "base_R", "R", "status", "resets", "last_reset", "window_ok", "wall",
                                               "check_judge", "check_dump_same")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
