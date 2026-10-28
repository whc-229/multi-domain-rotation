"""T3b 作业文件：不缩放的 R_LP 变体（fresh_rlp_u、fresh_noreset_u、steady_u）× 训练实例。
默认 m ≤ 6（TASKS：先跑 m ≤ 6）；加参数 big 生成 m ≥ 8 的部分（_u 明显好于原版时补跑）。
用法：python3 code/t3b_make_jobs.py [big] > results/T3b/jobs[_big].tsv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names  # noqa: E402
from t3_job import U_KINDS  # noqa: E402

J = os.path.join(ROOT, "results", "T3b", "jobs")
big = len(sys.argv) > 1 and sys.argv[1] == "big"
for n in sorted(train_names(), key=lambda n: (int(n[1:3]), n)):
    if (int(n[1:3]) >= 8) != big: continue
    for kind in U_KINDS:
        print(f"{J}/{n}__{kind}.json\t{4 * 3600}\tcd {ROOT} && python3 code/t3_job.py {n} {kind}")
