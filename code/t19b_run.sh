#!/bin/bash
# T19 (b) 正式运行：results/T19/jobs_b.tsv 每行 "<实例> <种子> <超时秒>"（m = 10、8 在前），12 个并行；已完成的自动跳过（t19b_job.py 自己判断）。
# 每个作业的输出写 results/T19/logs_b/<实例>__s<种子>.log。用法：bash code/t19b_run.sh
cd "$(dirname "$0")/.." || exit 1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/home/hk/miniconda3/envs/ana/bin/python3
xargs -P 12 -L 1 bash -c 'timeout "$2" '"$PY"' code/t19b_job.py "$0" "$1" > results/T19/logs_b/"$0"__s"$1".log 2>&1; echo "$(date "+%m-%d %H:%M:%S") $0 $1 退出码 $?"' < results/T19/jobs_b.tsv
echo "$(date "+%m-%d %H:%M:%S") 全部结束"
