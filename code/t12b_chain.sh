#!/bin/bash
# T12b 的严格窗口 ES 重跑：等 T12 调度器（PID 作为第 1 个参数）结束后，以高优先级（第三批）派出 results/T12b/es_jobs.tsv。
# 高优先级调度器从启动到全部派出都挂优先标记，第四批的低优先级调度器在这期间不派新作业。
# 用法：nohup bash code/t12b_chain.sh <T12 调度器 PID> > results/T12b_chain.log 2>&1 &
cd /home/hk/多域轮转 || exit 1
PY=/home/hk/miniconda3/envs/ana/bin/python3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
while kill -0 "$1" 2>/dev/null; do sleep 20; done
echo "$(date '+%m-%d %H:%M:%S') T12 调度器已结束，启动 T12b ES 调度器（hi）"
exec $PY code/sched.py results/T12b/logs hi 12 -- results/T12b/es_jobs.tsv
