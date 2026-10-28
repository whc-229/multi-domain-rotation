"""T4 作业文件：m ∈ {2, 3} 的 24 个训练实例。
- lb：H = 2·S_max，有无提前大修，限时 1 小时（2 小时试验只到 ⌈LP⌉，H = 3S 未做，见报告）；
- esub：30 代 ES，curdef / lpgplan × 种子 0–2（上界来源）。
用法：python3 code/t4_make_jobs.py > results/T4/jobs.tsv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names  # noqa: E402
J = os.path.join(ROOT, "results", "T4", "jobs")
names = train_names(m_list=[2, 3])
for n in names:
    for cfg in ("lpgplan", "curdef"):
        for sd in range(3):
            print(f"{J}/{n}__esub_{cfg}__s{sd}.json\t{4 * 3600}\tcd {ROOT} && python3 code/t4_job.py esub {n} {cfg} {sd}")
for rs in (0, 1):
    for n in names:
        print(f"{J}/{n}__lb__r{rs}__m2.json\t{3600 + 1800}\tcd {ROOT} && python3 code/t4_job.py lb {n} {rs} 2 3600")
