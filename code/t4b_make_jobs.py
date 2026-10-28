"""T4b 作业文件：m ∈ {2, 3} 的 24 个训练实例 × {curdef, lpgplan} × 种子 0–4，100 代 ES。
按实例规模从小到大排，先拿到 m = 2 的完整数据。
用法：python3 code/t4b_make_jobs.py > results/T4b/jobs.tsv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names  # noqa: E402

J = os.path.join(ROOT, "results", "T4b", "jobs")
for n in train_names(m_list=[2, 3]):
    for cfg in ("curdef", "lpgplan"):
        for sd in range(5):
            print(f"{J}/{n}__es100_{cfg}__s{sd}.json\t{6 * 3600}\tcd {ROOT} && python3 code/t4_job.py es100 {n} {cfg} {sd}")
