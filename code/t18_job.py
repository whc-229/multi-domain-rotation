"""T18 单个作业：RDD 启动窗口长度扫描（时间界的编码一侧）。一个 (训练实例, 偏离窗口倍数, 种子) 的 RDD-FD（rost5 rddfd）。
设置与 T12 / T16（β = 0）相同：默认规则 fresh_rlp_w（严格窗口 planW = 2·S_max、后备补足时域 hsfb）+ resetopt，
ρ = 0.1、q_on = 0.4、(4 + 12)、T_eval = 6·S_max、互换 0.5、stall 30；只改偏离窗口 W = ⌊倍数 · S_max⌋：
- 提前大修窗口固定为 [0, 2·S_max)：计划提前大修只在 t < 2·S_max 时执行，"本班后提前大修"备选项也只在 t < 2·S_max 时提供（rost5 的 T18 改动）；
- ρ 不变，所以槽数 K 随 W 内的决策数增长。
W = 0 就是规则本身（T12 base / T12b），W = 2·S_max 与 T16 的 β = 0 相同（复用），这里只跑倍数 0.5、1、3、4。
限时（墙钟）与 T16 相同：读 results/T16/tlim.json（m = 8、10），m = 4、6 为 300、600 秒。
作业内核对（同 T12）：最终解用 rost5 judge 重新判定逐字段相同；抽样 T_eval 缺岗与 rost5 eval 相同。
另记偏离位置：把最终基因型解码，用 Python 参考实现在最终台数下跑到 W，记下每个激活的偏离对应的窗口决策 (编号, 站, 时刻, u)，
并核对参考实现数出的窗口决策数与 rost5 eval 相同。
结果写进 results/T18/jobs/<实例>__w<倍数>__s<种子>.json（原子写入，已存在就跳过；运行时加 fcntl 锁）。
用法：python3 code/t18_job.py <实例> <倍数> <种子>"""
import os, sys, json, time, subprocess, fcntl, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
from t11_job import write_geno  # noqa: E402
from t12_job import rule_opts, RULES  # noqa: E402
from t16_job import TLIM  # noqa: E402   与 T16 相同的限时（含 tlim.json）
import rost5_ref as ref  # noqa: E402

R5 = os.environ.get("ROST5_BIN", os.path.join(ROOT, "code", "cpp", "rost5"))
OUT = os.environ.get("T18_OUT", os.path.join(ROOT, "results", "T18", "jobs"))


def dev_positions(I, name, rec, P2):
    """最终解的激活偏离：[(决策编号, 站, 时刻, u)]，以及参考实现数出的窗口决策数。"""
    D = ref.decode_geno(rec["genes"], rec["W"], rec["ndev"], rec["qon"], True)
    R = ref.parse_rlp(I, f"{RULES}/{name}__fresh_rlp_w.rlp"); PF = ref.parse_planf(I, f"{RULES}/{name}__fresh_rlp_w.planf")
    tr = {}
    r = ref.run(I, rec["Nk"], rec["Mj"], rec["W"], False, rlp=R, planf=PF, rdd=D, plan_until=P2, trace=tr)
    decs = tr.get("decs", [])
    pos = [(c, decs[c][0], decs[c][1], u) for c, u in enumerate(D["dev"]) if u >= 0 and c < len(decs)]
    return pos, r["ndec_w"]


def main():
    """跑 rddfd、作业内核对、记偏离位置，原子写入。"""
    name, mult, seed = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
    assert mult not in (0.0, 2.0)
    os.makedirs(OUT, exist_ok=True)
    jn = f"{name}__w{mult:g}__s{seed}"; path = os.path.join(OUT, jn + ".json")
    if os.path.exists(path): print("已完成，跳过", jn); return
    lk = open(path + ".lock", "w")
    try: fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError: print("另一个进程正在做", jn); sys.exit(4)
    if os.path.exists(path): print("已完成，跳过", jn); return
    I = load_train(name); _, txt = inst_paths(name, "train"); S = S_max(I); P2 = 2 * S; W = int(mult * S); m = I.meta["m"]
    opts = rule_opts(name, S)                    # planW = 2·S_max 固定
    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix=f"t18_{jn}_") as tmpd:
        dump = os.path.join(tmpd, "dump.tsv")
        args = [R5, "rddfd", txt, "365", "60000", "600000", str(6 * S), str(10 ** 9), str(TLIM[m]), "4", "12", str(seed), "0.1", "0.4"] + opts + \
               [f"W={W}", "resetopt", "stall=30", "swap=0.5", "fdsel=best", f"fddump={dump}"]
        r = subprocess.run(args, capture_output=True, text=True)
        try: rec = json.loads(r.stdout)
        except Exception: print("rost5 输出无法解析", r.stdout[-300:], r.stderr[-300:]); sys.exit(2)
        rec.update(inst=name, m=m, family=I.meta["family"], mult=mult, planW=P2, S=S, wall=round(time.time() - t0, 1), cmd=" ".join(args[:-1]))
        ev_opts = [o for o in opts if not o.startswith("hsfb=")]
        if rec["R"] >= 0:
            gp = os.path.join(tmpd, "best.geno"); write_geno(gp, rec["W"], rec["ndev"], rec["qon"], rec["genes"], 1)
            j = json.loads(subprocess.run([R5, "judge", txt, "60000", "600000", ",".join(map(str, rec["Nk"])), ",".join(map(str, rec["Mj"]))]
                                          + ev_opts + [f"W={rec['W']}", "resetopt", f"rdd={gp}"], capture_output=True, text=True, check=True).stdout)
            rec["check_judge"] = int(j["R"] == rec["R"] and j["feasible"] == 1 and j["status"] == rec["status"] and j["verify_gaps"] == rec["verify_gaps"]
                                     and j["verify_cert"] == rec["verify_cert"] and j["verify_period"] == rec["verify_period"] and j["exact"] == rec["exact"]
                                     and j["resets"] == rec["resets"] and j["last_reset"] == rec["last_reset"])
            # 偏离位置（参考实现），并核对窗口决策数与 rost5 eval 相同
            pos, nd_py = dev_positions(I, name, rec, P2)
            e = json.loads(subprocess.run([R5, "eval", txt, str(rec["W"]), "0"] + [str(x) for x in rec["Nk"]] + ["|"] + [str(x) for x in rec["Mj"]]
                                          + ["|", "|"] + ev_opts + [f"rdd={gp}"], capture_output=True, text=True, check=True).stdout)
            rec.update(dev_pos=pos, ndec_final=e["ndec_w"], check_ndec=int(nd_py == e["ndec_w"]))
        n_ok = n_all = 0
        for line in open(dump):
            Nk, Mj, g, genes = line.rstrip("\n").split("\t")
            sp = os.path.join(tmpd, "s.geno"); write_geno(sp, rec["W"], rec["ndev"], rec["qon"], [float(x) for x in genes.split(",")], 1)
            e = json.loads(subprocess.run([R5, "eval", txt, str(rec["Teval"]), "0"] + [str(x) for x in json.loads(Nk)] + ["|"] + [str(x) for x in json.loads(Mj)]
                                          + ["|", "|"] + ev_opts + [f"rdd={sp}"], capture_output=True, text=True, check=True).stdout)
            n_all += 1; n_ok += int(e["gaps"] == int(g))
        rec.update(check_dump_n=n_all, check_dump_same=n_ok)
    rec["proven"] = int(rec.get("verify_cert") == 1 and rec.get("exact") == 1)
    rec["feasible"] = int(rec.get("R", -1) >= 0 and (rec["proven"] == 1 or rec.get("verify_gaps") == 0))
    rec["window_ok"] = int(rec.get("last_reset", -1) < P2)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: rec.get(k) for k in ("inst", "mult", "seed", "base_R", "R", "status", "feasible", "n_active", "wall", "check_judge", "check_dump_same", "check_ndec")},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
