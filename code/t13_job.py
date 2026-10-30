"""T13 / 第五批的单个作业：一个 (方法, 实例, 种子)，限时 TL 秒（墙钟）。所有方法共用同一个外壳，见 code/t13_methods.py。
- 规则类（fifo、pb、gphh）：rost5 size（定台数 + 60 万格复核），与种子无关；
- 搜索类（esdw、eslp、rk、d5）：rost5 fd（表示 w / rk / d5），(4 + 12)，T_eval = 6·S_max，互换 0.5，stall 30，限时 TL；
- rddhand：rost5 rddfd，默认规则 PB，resetopt（只在窗口内），ρ = 0.1；
- rhmip：code/t13_rhmip.py（每步限时 TL/8，总限时 TL）；
- ours：rost5 rddfd，fresh_rlp_w（严格窗口、后备补足时域）+ resetopt，与 T12 相同。
搜索类与 rddhand、ours 在作业内核对（同 T11）：最终解用 rost5 judge 重新判定逐字段相同；抽样的 T_eval 缺岗与 rost5 eval 相同。
结果写进 <输出目录>/<方法>__<实例>__s<种子>.json（原子写入，已存在就跳过；运行时加 fcntl 锁）。
用法：python3 code/t13_job.py <方法> <实例> <种子> <限时秒> [输出目录] [GP 树文件（gphh 用）]"""
import os, sys, json, time, subprocess, fcntl, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import t13_methods as TM  # noqa: E402
from t11_job import write_geno  # noqa: E402

R5 = os.environ.get("ROST5_BIN", os.path.join(ROOT, "code", "cpp", "rost5"))


def eval_args(opts):
    """rost5 eval 的第 3、4 段：权重（w=）放第 3 段，其余选项放第 4 段（去掉只对定台数有用的 facloop、hsfb、N0、M0）。"""
    w = [o for o in opts if o.startswith("w=")]
    rest = [o for o in opts if not o.startswith(("w=", "hsfb=", "N0=", "M0=")) and o != "facloop"]
    return (w[0][2:].split(",") if w else []) + ["|"] + rest


def checks(rec, txt, opts, dump, tmpd, method, W):
    """作业内核对：judge 重判最终解；抽样缺岗与 eval 相同。"""
    out = {}
    if method in ("rddhand", "ours"):
        gp = os.path.join(tmpd, "best.geno"); write_geno(gp, rec["W"], rec["ndev"], rec["qon"], rec["genes"], 1); extra = ["resetopt", f"rdd={gp}"]
        geno_opt = lambda g: (write_geno(os.path.join(tmpd, "s.geno"), rec["W"], rec["ndev"], rec["qon"], g, 1), [f"rdd={os.path.join(tmpd, 's.geno')}"])[1]
    elif method in ("rk", "d5"):
        extra = [f"pd={method}", "pdg=" + ",".join(repr(x) for x in rec["genes"])]
        geno_opt = lambda g: [f"pd={method}", f"pdW={W}", "pdg=" + ",".join(repr(x) for x in g)]
    else:   # w
        extra = ["w=" + ",".join(repr(x) for x in rec["genes"])]
        geno_opt = None
    if rec["R"] >= 0:
        j = json.loads(subprocess.run([R5, "judge", txt, "60000", "600000", ",".join(map(str, rec["Nk"])), ",".join(map(str, rec["Mj"]))]
                                      + opts + [f"W={W}"] + extra, capture_output=True, text=True, check=True).stdout)
        same = (j["R"] == rec["R"] and j["feasible"] == 1 and j["status"] == rec["status"] and j["verify_gaps"] == rec["verify_gaps"]
                and j["verify_cert"] == rec["verify_cert"] and j["verify_period"] == rec["verify_period"] and j["exact"] == rec["exact"]
                and j["resets"] == rec["resets"] and j["last_reset"] == rec["last_reset"])
        out.update(check_judge=int(same), judge=j)
    n_ok = n_all = 0
    for line in open(dump):
        Nk, Mj, g, genes = line.rstrip("\n").split("\t"); genes = [float(x) for x in genes.split(",")]
        if geno_opt is None:   # w：权重放第 3 段
            ea = eval_args([o for o in opts if not o.startswith("w=")]); ea = [repr(x) for x in genes] + ea
        else:
            ea = eval_args(opts) + geno_opt(genes)
        args = [R5, "eval", txt, str(rec["Teval"]), "0"] + [str(x) for x in json.loads(Nk)] + ["|"] + [str(x) for x in json.loads(Mj)] + ["|"] + ea
        e = json.loads(subprocess.run(args, capture_output=True, text=True, check=True).stdout)
        n_all += 1; n_ok += int(e["gaps"] == int(g))
    out.update(check_dump_n=n_all, check_dump_same=n_ok)
    return out


def main():
    method, name, seed, TL = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4])
    out = sys.argv[5] if len(sys.argv) > 5 else os.path.join(ROOT, "results", "T13", "jobs")
    tree = sys.argv[6] if len(sys.argv) > 6 else None
    os.makedirs(out, exist_ok=True)
    jn = f"{method}__{name}__s{seed}"; path = os.path.join(out, jn + ".json")
    if os.path.exists(path): print("已完成，跳过", jn); return
    lk = open(path + ".lock", "w")
    try: fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError: print("另一个进程正在做", jn); sys.exit(4)
    if os.path.exists(path): print("已完成，跳过", jn); return
    I = load_train(name); S = S_max(I); W = 2 * S; _, txt = inst_paths(name, "train")
    opts = TM.opts(method, name, S, tree) if method != "rhmip" else []
    t0 = time.time()
    common = ["365", "60000", "600000"]
    if method in ("fifo", "pb", "gphh"):
        rec = json.loads(subprocess.run([R5, "size", txt] + common + opts, capture_output=True, text=True, check=True).stdout)
    elif method == "rhmip":
        d = os.path.join(out, "rhmip"); os.makedirs(d, exist_ok=True)
        subprocess.run([sys.executable, os.path.join(ROOT, "code", "t13_rhmip.py"), name, str(TL / 8), str(TL), d], capture_output=True, text=True, check=True, env=os.environ)
        rec = json.load(open(os.path.join(d, f"{name}.json")))
    else:
        with tempfile.TemporaryDirectory(prefix=f"t13_{jn}_") as tmpd:
            dump = os.path.join(tmpd, "dump.tsv")
            fdo = [f"W={W}", "stall=30", "swap=0.5", "fdsel=best", f"fddump={dump}"]
            if method in ("rddhand", "ours"):
                args = [R5, "rddfd", txt] + common + [str(6 * S), str(10 ** 9), str(TL), "4", "12", str(seed), "0.1", "0.4"] + opts + fdo + ["resetopt"]
            else:
                args = [R5, "fd", txt] + common + [str(6 * S), str(10 ** 9), str(TL), "4", "12", str(seed), TM.SEARCH[method], "0.1", "0.4"] + opts + fdo
            r = subprocess.run(args, capture_output=True, text=True)
            try: rec = json.loads(r.stdout)
            except Exception: print("rost5 输出无法解析", r.stdout[-300:], r.stderr[-300:]); sys.exit(2)
            rec["Teval"] = rec.get("Teval", 6 * S)
            rec.update(checks(rec, txt, opts, dump, tmpd, method, W)); rec["cmd"] = " ".join(args[:-1])
    rec.update(method=method, inst=name, m=I.meta["m"], family=I.meta["family"], seed=seed, S=S, W=W, TL=TL, wall=round(time.time() - t0, 1))
    rec["proven"] = int(rec.get("verify_cert") == 1 and rec.get("exact") == 1)
    rec["feasible"] = int(rec.get("R", -1) >= 0 and (rec["proven"] == 1 or rec.get("verify_gaps") == 0))
    rec["window_ok"] = int(rec.get("last_reset", -1) < W)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: rec.get(k) for k in ("method", "inst", "seed", "R", "status", "feasible", "window_ok", "wall", "check_judge", "check_dump_same")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
