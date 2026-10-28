"""T8 单个作业：一个 (训练实例, 默认规则, 偏离率 ρ, 种子) 的 RDD 进化。结果写进 results/T8/jobs/<名>.json（原子写入，已存在就跳过）。
命令：rost5 rdd <实例> 365 60000 600000 100 4 12 <种子> <ρ> 0.4 <规则选项> W=2·S_max
  - (μ + λ) = (4 + 12)，100 代，初始 16 个个体（个体 0 = 全部槽关闭，即默认规则本身），共 1216 次评估，与 T4b 的 ES 相同；
  - K = ⌈ρ·N_dev / q_on⌉，q_on = 0.4；每个基因以 1/(2K) 的概率重抽；子代没有基因变化时直接沿用父代的评估结果；
  - 适应度：资源数（判不出循环但无缺岗时 +3，与 ES 的 certs 模式相同）；同资源数时比 near_miss（任减 1 个资源后 365 格内缺岗数的最小值）；
  - 场景：全新启动、不提前大修（不加 resetopt）。
默认规则见 code/t8_rules.py（rlp / es / hand）。输出里 base_* 是默认规则本身（同一次运行里算出，与种子无关）。
作业运行时对结果文件加 fcntl 锁；拿不到锁就以返回码 4 退出、不写文件。
用法：python3 code/t8_job.py <实例> <rlp|es|hand> <ρ> <种子>"""
import os, sys, json, time, subprocess, fcntl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import t8_rules  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")
OUT = os.path.join(ROOT, "results", "T8", "jobs")


def main():
    name, rule, rho, seed = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4])
    os.makedirs(OUT, exist_ok=True)
    jn = f"{name}__{rule}__r{rho:g}__s{seed}"; path = os.path.join(OUT, jn + ".json")
    if os.path.exists(path): print("已完成，跳过", jn); return
    lk = open(path + ".lock", "w")      # 两个调度器从作业表两头跑时，同一作业只让一个进程做（单写入者）
    try: fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError: print("另一个进程正在做", jn); sys.exit(4)
    if os.path.exists(path): print("已完成，跳过", jn); return
    I = load_train(name); _, txt = inst_paths(name, "train"); W = 2 * S_max(I)
    opts = t8_rules.rule_opts(name, rule)
    src = t8_rules.best_es(name, "curdef")[1] if rule == "es" else None
    args = [R5, "rdd", txt, "365", "60000", "600000", "100", "4", "12", str(seed), repr(rho), "0.4"] + opts + [f"W={W}"]
    t0 = time.time()
    r = subprocess.run(args, capture_output=True, text=True)
    try: rec = json.loads(r.stdout)
    except Exception: print("rost5 输出无法解析", r.stdout[-300:], r.stderr[-300:]); sys.exit(2)
    rec.update(inst=name, m=I.meta["m"], family=I.meta["family"], rule=rule, es_src=os.path.relpath(src, ROOT) if src else None,
               wall=round(time.time() - t0, 1), cmd=" ".join(args))
    rec["proven"] = int(rec.get("verify_cert") == 1 and rec.get("exact") == 1)
    rec["base_proven"] = int(rec.get("base_verify_cert") == 1 and rec.get("base_exact") == 1)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: rec.get(k) for k in ("inst", "rule", "rho", "seed", "base_R", "R", "status", "ndev", "K", "n_active", "wall")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
