#!/bin/bash
# T19 正式运行：results/T19_jobs.tsv 每行 "<steady|best> <实例>"，12 个并行；已完成的自动跳过（t19_connector.py 自己判断）。
# 每个作业的输出写 results/T19/logs/<来源>__<实例>.log。用法：bash code/t19_run.sh
cd "$(dirname "$0")/.." || exit 1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/home/hk/miniconda3/envs/ana/bin/python3
xargs -P 12 -L 1 bash -c 'timeout 7200 '"$PY"' code/t19_connector.py "$0" "$1" > results/T19/logs/"$0"__"$1".log 2>&1; echo "$(date +%H:%M:%S) $0 $1 退出码 $?"' < results/T19_jobs.tsv
echo "$(date +%H:%M:%S) 全部结束"
