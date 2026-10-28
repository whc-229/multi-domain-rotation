"""T3 单个作业。结果写进 results/T3/jobs/<名>.json（原子写入，已存在就跳过）。
种类：
- R_LP（不搜索）：fresh_rlp（带 LP 提前大修计划）、fresh_noreset、steady，以及诊断变体 fresh_rlp_c、fresh_noreset_c（缩放分母用 ⌈N_LP⌉）；
  命令：rost5 rlp <实例> <规则文件> 365 60000 600000 [planf=…]（定台数参数与 stats_run4.py 的 ES 相同）。
- ES 参照：es_curdef、es_lpgplan，命令与 stats_run4.py / run_lpg.sh 相同：rost5 es <实例> certs 365 10 4 12 <种子> 60000 600000 <选项>。
依赖：T2 的 rate、fr m4、frR m4 结果。缺依赖时退出码 3、不写结果，重跑调度器时会再试。
规则文件由 code/t3_make_rlp.py 生成；同一实例的多个作业用文件锁保证只有一个进程生成。
用法：python3 code/t3_job.py <实例名> <种类> [种子]"""
import os, sys, json, time, fcntl, shutil, subprocess, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
import t3_make_rlp  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")
OUT = os.path.join(ROOT, "results", "T3", "jobs")
RDIR = t3_make_rlp.OUT
RLP_KINDS = ["fresh_rlp", "fresh_noreset", "steady", "fresh_rlp_c", "fresh_noreset_c"]


def deps_ok(name):
    j = t3_make_rlp.T2J; d = t3_make_rlp.T2D
    need = [os.path.join(j, f"{name}__rate__m0.json")] + \
           [t3_make_rlp.dual_path(name, k, 4, d) for k in ("fr", "frR")]
    return all(p is not None and os.path.exists(p) for p in need)


def ensure_rules(name):
    """生成该实例的规则文件（加锁；先写到临时目录再逐个改名，避免半个文件）。"""
    os.makedirs(RDIR, exist_ok=True)
    mark = os.path.join(RDIR, f"{name}.done")
    if os.path.exists(mark): return
    with open(os.path.join(RDIR, f"{name}.lock"), "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        if os.path.exists(mark): return
        tmp = tempfile.mkdtemp(dir=RDIR, prefix=f".{name}.")
        info = t3_make_rlp.make(name, out=tmp)
        for f in os.listdir(tmp): os.replace(os.path.join(tmp, f), os.path.join(RDIR, f))
        shutil.rmtree(tmp)
        with open(mark, "w") as f: f.write(json.dumps(info) + "\n")


def main():
    name, kind = sys.argv[1], sys.argv[2]
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else None
    os.makedirs(OUT, exist_ok=True)
    jn = f"{name}__{kind}" + (f"__s{seed}" if seed is not None else "")
    path = os.path.join(OUT, jn + ".json")
    if os.path.exists(path): print("已完成，跳过", jn); return
    if not deps_ok(name): print("T2 依赖未齐", name); sys.exit(3)
    ensure_rules(name)
    I = load_train(name); _, txt = inst_paths(name, "train")
    if kind in RLP_KINDS:
        args = [R5, "rlp", txt, os.path.join(RDIR, f"{name}__{kind}.rlp"), "365", "60000", "600000"]
        if kind.startswith("fresh_rlp"):
            args.append("planf=" + os.path.join(RDIR, f"{name}__{kind}.planf"))     # fresh_rlp → .planf，fresh_rlp_c → _c.planf
    elif kind == "es_curdef":
        args = [R5, "es", txt, "certs", "365", "10", "4", "12", str(seed), "60000", "600000", "curdef"]
    elif kind == "es_lpgplan":
        args = [R5, "es", txt, "certs", "365", "10", "4", "12", str(seed), "60000", "600000",
                "lpg=" + os.path.join(RDIR, f"{name}.lpg.txt"), "planf=" + os.path.join(RDIR, f"{name}__fresh_rlp.planf")]
    else:
        sys.exit(f"未知种类 {kind}")
    t0 = time.time()
    r = subprocess.run(args, capture_output=True, text=True)
    try: d = json.loads(r.stdout)
    except Exception: print("rost5 输出无法解析", r.stdout[-300:], r.stderr[-300:]); sys.exit(2)
    d.update(inst=name, m=I.meta["m"], family=I.meta["family"], kind=kind, seed=seed, wall=round(time.time() - t0, 1), cmd=" ".join(args))
    d["proven"] = int(d.get("verify_cert") == 1 and d.get("exact") == 1)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(d, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: d.get(k) for k in ("inst", "kind", "seed", "R", "status", "verify_cert", "verify_gaps", "wall")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
