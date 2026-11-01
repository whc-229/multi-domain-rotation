#!/bin/bash
# T12c 排队（11-01 审核：T16 之后、T17 之前）：等主调度器把 T16 的作业全部派出（每个作业要么已完成，要么正在运行），
# 再以 hihold 优先级派出 T12c 的 24 个作业。hihold 的标记挂到 T12c 全部结束，期间 lo 调度器不派新作业（T17 等它）。
# 用法：nohup bash code/t12c_chain.sh > results/T12c_chain.log 2>&1 &
cd /home/hk/多域轮转 || exit 1
PY=/home/hk/miniconda3/envs/ana/bin/python3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
# t16_all_dispatched：T16 作业表里每一行的完成标志已存在，或者它的命令正在运行
until $PY - <<'EOF'
import os, sys
sys.path.insert(0, "code")
from sched import procs_args
run = procs_args()
for jf in ("results/T16/jobs_m46.tsv", "results/T16/jobs_m810.tsv"):
    for line in open(jf):
        done, _, cmd = line.rstrip("\n").split("\t")
        if not os.path.exists(done) and cmd.split("&& ")[-1].split(None, 1)[1].strip() not in run: sys.exit(1)
sys.exit(0)
EOF
do sleep 120; done
echo "$(date '+%m-%d %H:%M:%S') T16 已全部派出，启动 T12c"
exec $PY code/sched.py results/T12c/logs hihold 12 -- results/T12c/jobs.tsv
