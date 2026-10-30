"""T11 单个作业：一个 (训练实例, 默认规则, 种子, 设置) 的 RDD-FD（rost5 rddfd）。结果写进 results/T11/jobs/<名>.json（原子写入，已存在就跳过）。
命令：rost5 rddfd <实例> 365 60000 600000 <T_eval> <代数> <限时> 4 12 <种子> 0.1 0.4 <规则选项> W=2·S_max stall=30 swap=<p> fdsel=<best|base>
  - 默认规则、定台数起点、N_dev、K、(μ + λ)、变异与 T8 相同（code/t8_rules.py，不提前大修场景）；
  - 设置（见 SETS）：main = T_eval 6·S_max、互换概率 0.5、按最近可行个体选资源；其余三个各改一项，作敏感性对照。
作业内核对（结果里的 check_*）：
  1. 最终解用 rost5 judge（size 的判定部分：6 万格判环 → 24 万格 → 60 万格复核）重新判定一次，R、可行、状态、复核缺岗、复核判环、周期、exact 逐个相同；
  2. rddfd 抽样记下的（台数, 基因型, T_eval 缺岗）用 rost5 eval（同台数、同基因型、T_eval 格、不判环）重算，缺岗相同。
作业运行时对结果文件加 fcntl 锁；拿不到锁就以返回码 4 退出、不写文件。
用法：python3 code/t11_job.py <实例> <rlp|es> <种子> <设置>"""
import os, sys, json, time, subprocess, fcntl, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import t8_rules  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")
OUT = os.environ.get("T11_OUT", os.path.join(ROOT, "results", "T11", "jobs"))   # 环境变量只用于开发期试算
GENS, TLIM = 1000, 36000
# 设置名 → (T_eval 倍数 × S_max, 互换概率, 选资源用的基因型)
SETS = {"main": (6, 0.5, "best"), "te4": (4, 0.5, "best"), "noswap": (6, 0.0, "best"), "selbase": (6, 0.5, "base")}


def write_geno(path, W, ndev, qon, genes, reset_opt=0):
    """写 rost5 的基因型文件（与 read_geno 的格式相同）。"""
    with open(path, "w") as f:
        f.write(f"W {W}\nndev {ndev}\nqon {qon!r}\nresetopt {reset_opt}\n")
        for a in range(0, len(genes) - 1, 2): f.write(f"slot {genes[a]!r} {genes[a + 1]!r}\n")


def checks(rec, txt, opts, dump, tmpd):
    """作业内的两项核对，返回 dict。"""
    out = {}
    gp = os.path.join(tmpd, "best.geno")
    write_geno(gp, rec["W"], rec["ndev"], rec["qon"], rec["genes"])
    if rec["R"] >= 0:
        j = json.loads(subprocess.run([R5, "judge", txt, "60000", "600000", ",".join(map(str, rec["Nk"])), ",".join(map(str, rec["Mj"]))]
                                      + opts + [f"W={rec['W']}", f"rdd={gp}"], capture_output=True, text=True, check=True).stdout)
        same = (j["R"] == rec["R"] and j["feasible"] == 1 and j["status"] == rec["status"] and j["verify_gaps"] == rec["verify_gaps"]
                and j["verify_cert"] == rec["verify_cert"] and j["verify_period"] == rec["verify_period"] and j["exact"] == rec["exact"])
        out.update(check_judge=int(same), judge=j)
    n_ok = n_all = 0
    for line in open(dump):
        Nk, Mj, g, genes = line.rstrip("\n").split("\t")
        genes = [float(x) for x in genes.split(",")]
        sp = os.path.join(tmpd, "s.geno"); write_geno(sp, rec["W"], rec["ndev"], rec["qon"], genes)
        args = [R5, "eval", txt, str(rec["Teval"]), "0"] + [str(x) for x in json.loads(Nk)] + ["|"] + [str(x) for x in json.loads(Mj)] + ["|"]
        # eval 的权重与选项：线性规则的权重放第 3 段，其余选项（rlp=、curdef、lpg=、planf=）放第 4 段
        w = [o for o in opts if o.startswith("w=")]
        args += (w[0][2:].split(",") if w else []) + ["|"] + [o for o in opts if not o.startswith("w=") and o != "facloop"] + [f"rdd={sp}"]
        e = json.loads(subprocess.run(args, capture_output=True, text=True, check=True).stdout)
        n_all += 1; n_ok += int(e["gaps"] == int(g))
    out.update(check_dump_n=n_all, check_dump_same=n_ok)
    return out


def main():
    name, rule, seed, st = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
    tm, swap, sel = SETS[st]
    os.makedirs(OUT, exist_ok=True)
    jn = f"{name}__{rule}__s{seed}__{st}"; path = os.path.join(OUT, jn + ".json")
    if os.path.exists(path): print("已完成，跳过", jn); return
    lk = open(path + ".lock", "w")
    try: fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError: print("另一个进程正在做", jn); sys.exit(4)
    if os.path.exists(path): print("已完成，跳过", jn); return
    I = load_train(name); _, txt = inst_paths(name, "train"); S = S_max(I); W = 2 * S
    opts = t8_rules.rule_opts(name, rule)
    with tempfile.TemporaryDirectory(prefix=f"t11_{jn}_") as tmpd:
        dump = os.path.join(tmpd, "dump.tsv")
        args = [R5, "rddfd", txt, "365", "60000", "600000", str(tm * S), str(GENS), str(TLIM), "4", "12", str(seed), "0.1", "0.4"] + opts + \
               [f"W={W}", "stall=30", f"swap={swap}", f"fdsel={sel}", f"fddump={dump}"]
        t0 = time.time()
        r = subprocess.run(args, capture_output=True, text=True)
        try: rec = json.loads(r.stdout)
        except Exception: print("rost5 输出无法解析", r.stdout[-300:], r.stderr[-300:]); sys.exit(2)
        rec.update(inst=name, m=I.meta["m"], family=I.meta["family"], rule=rule, setting=st, S=S, wall=round(time.time() - t0, 1), cmd=" ".join(args[:-1]))
        rec["proven"] = int(rec.get("verify_cert") == 1 and rec.get("exact") == 1)
        rec["base_proven"] = int(rec.get("base_verify_cert") == 1 and rec.get("base_exact") == 1)
        t1 = time.time(); rec.update(checks(rec, txt, opts, dump, tmpd)); rec["check_sec"] = round(time.time() - t1, 1)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: rec.get(k) for k in ("inst", "rule", "seed", "setting", "base_R", "R", "status", "wall", "check_judge", "check_dump_same", "check_dump_n")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
