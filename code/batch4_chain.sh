#!/bin/bash
# 第四批的后续排队（10-31 晚，按"第四批补充二"）：等 T12 调度器（PID 作为第 1 个参数）结束后，
# 1) 按 T12 在 m = 8、10 上的收敛定 T16 / T18 的限时（code/t16_tlim.py → results/T16/tlim.json）；
# 2) 启动两个低优先级调度器：T15（同时最多 4 个）；以及 T16 → T17 → T18 → T13 冒烟 → T14（按文件顺序派出）。
# 用法：nohup bash code/batch4_chain.sh <T12 调度器 PID> > results/batch4_chain.log 2>&1 &
cd /home/hk/多域轮转 || exit 1
PY=/home/hk/miniconda3/envs/ana/bin/python3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
while kill -0 "$1" 2>/dev/null; do sleep 120; done
echo "$(date '+%m-%d %H:%M:%S') T12 调度器已结束"
$PY code/t16_tlim.py || exit 1
nohup $PY code/sched.py results/batch4_logs/T15 lo 12 4 -- results/T15/jobs.tsv > results/batch4_sched_T15.log 2>&1 &
echo "T15 调度器 PID $!"
nohup $PY code/sched.py results/batch4_logs/main lo 12 -- results/T16/jobs_m46.tsv results/T16/jobs_m810.tsv results/T17/jobs.tsv \
    results/T18/jobs_m46.tsv results/T18/jobs_m810.tsv results/T13/smoke_jobs.tsv results/T14_gphh/jobs.tsv > results/batch4_sched_main.log 2>&1 &
echo "主调度器 PID $!"
