#!/bin/bash
# T20 调度第二版（11-06 审核后）：等 T20-1 的 LP 跑完（results/T20/lb/run.log 出现"全部结束"），再派 T20-2 的全部对手（所有 m，新外壳，12 个并行）。
# 第一版 code/t20_chain.sh 的 m ≤ 6 队列没开始就停了（11-06 审核：改外壳后重新生成作业表），作业表 results/T20/opp/jobs_m2-6.tsv 不再用。
# 作业表 jobs_all.tsv 先 m ≤ 6（2160 行）、后 m ≥ 8（2160 行），各自按预算从大到小；m ≥ 8 排后面，开始前还能按"待确认"调整。
# 断了重跑同一条命令即可（已完成的作业跳过；T20-1 已结束时直接派 T20-2）。
cd "$(dirname "$0")/.." || exit 1
until grep -q "全部结束" results/T20/lb/run.log 2>/dev/null; do sleep 120; done
bash code/t20_run.sh results/T20/opp/jobs_all.tsv 12 results/T20/opp/logs results/T20/opp/run_all.log
