"""生成 T3 作业文件。R_LP 5 种 × 96 个训练实例在前（按站数从大到小）；ES 参照（curdef、lpgplan）× 3 个种子 × m ≤ 10 的训练实例在后，
按站数从小到大（ES 在大实例上很贵，先拿到小规模的完整数据）。
用法：python3 code/t3_make_jobs.py > results/T3/jobs.tsv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names  # noqa: E402
from t3_job import RLP_KINDS  # noqa: E402

J = os.path.join(ROOT, "results", "T3", "jobs")
names = train_names()
for n in sorted(names, key=lambda n: (-int(n[1:3]), n)):
    for kind in RLP_KINDS:
        print(f"{J}/{n}__{kind}.json\t{4 * 3600}\tcd {ROOT} && python3 code/t3_job.py {n} {kind}")
for n in sorted(names, key=lambda n: (int(n[1:3]), n)):
    if int(n[1:3]) > 10: continue
    for seed in range(3):
        for kind in ("es_curdef", "es_lpgplan"):
            print(f"{J}/{n}__{kind}__s{seed}.json\t{8 * 3600}\tcd {ROOT} && python3 code/t3_job.py {n} {kind} {seed}")
