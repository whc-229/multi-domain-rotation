"""T6 作业文件。
- lp：96 个训练实例的 LP_fr^Rw(4·S_max)，按 m 从小到大排（fresh_rlp_w 先要 m ≤ 6 的解）；每个限时 8 小时。
- rate 作业很快（每个 < 1 秒），不进调度器，直接循环跑：python3 code/t6_make_jobs.py rate。
用法：python3 code/t6_make_jobs.py lp > results/T6/jobs_lp.tsv"""
import os, sys, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names  # noqa: E402

J = os.path.join(ROOT, "results", "T6", "jobs")


def by_m():
    """训练实例按 (m, 族, 序号) 排序。"""
    return sorted(train_names(), key=lambda n: (int(n[1:3]), n))


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "lp":
        for n in by_m():
            print(f"{J}/{n}__frRw__m4.json\t{8 * 3600}\tcd {ROOT} && python3 code/t6_job.py {n} frRw 4")
    elif what == "rate":
        for n in by_m():
            subprocess.run([sys.executable, os.path.join(ROOT, "code", "t6_job.py"), n, "rate"], check=True, stdout=subprocess.DEVNULL)
        print("rate 完成")
