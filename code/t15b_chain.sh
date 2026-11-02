#!/bin/bash
# T15 两项小检查的排队（11-02 审核：与 T14 并行，排在 T13 冒烟之后）：
# 等主调度器把 T14 的作业全部派出（每个作业要么已完成，要么正在运行），再以 lo 优先级派出 results/T15/jobs_b.tsv
# （2 个 plan 作业 + 120 个长预算 RDD-FD），全机上限 12。两项检查都经由 t15_job.py 启动，正在运行的调度器都数得到。
# 用法：nohup bash code/t15b_chain.sh > results/T15/chain_b.log 2>&1 &
cd /home/hk/多域轮转 || exit 1
PY=/home/hk/miniconda3/envs/ana/bin/python3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
# t14_all_dispatched：T14 作业表里每一行的完成标志已存在，或者它的命令正在运行
until $PY - <<'PYEOF'
import os, sys
sys.path.insert(0, "code")
from sched import procs_args
run = procs_args()
for line in open("results/T14_gphh/jobs.tsv"):
    done, _, cmd = line.rstrip("\n").split("\t")
    if not os.path.exists(done) and cmd.split("&& ")[-1].split(None, 1)[1].strip() not in run: sys.exit(1)
sys.exit(0)
PYEOF
do sleep 120; done
echo "$(date '+%m-%d %H:%M:%S') T14 已全部派出，启动 T15 小检查"
exec $PY code/sched.py results/T15/logs_b lo 12 -- results/T15/jobs_b.tsv
