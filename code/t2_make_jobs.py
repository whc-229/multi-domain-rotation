"""生成 T2 作业文件：96 个训练实例 × {rate, fr×{2,4}, frR×{2,4}}。按站数从大到小排（长作业先跑，总时长最短）。
用法：python3 code/t2_make_jobs.py > results/T2/jobs.tsv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names, load_train  # noqa: E402

LIM = 6 * 3600 + 1800     # 求解器内部限时 6 小时，外层多留半小时给建模和证书
names = sorted(train_names(), key=lambda n: (-int(n[1:3]), n))
for kind, mult in [("frR", 4), ("fr", 4), ("frR", 2), ("fr", 2), ("rate", 0)]:
    for n in names:
        done = os.path.join(ROOT, "results", "T2", "jobs", f"{n}__{kind}__m{mult:g}.json")
        print(f"{done}\t{LIM}\tcd {ROOT} && python3 code/t2_job.py {n} {kind} {mult}")
