#!/bin/bash
# T20 主调度（11-06）：先 T20-1 的 LP（9 个并行，与冒烟补跑、核对错开），再 T20-2 的对手（m ≤ 6，12 个并行）。
# m ≥ 8 的对手等"待确认"里共用外壳的问题定了再派（见 报告/T20_0_冒烟.md）。断了重跑同一条命令即可（已完成的跳过）。
cd "$(dirname "$0")/.." || exit 1
bash code/t20_run.sh results/T20/lb/jobs.tsv 9 results/T20/lb/logs results/T20/lb/run.log
bash code/t20_run.sh results/T20/opp/jobs_m2-6.tsv 12 results/T20/opp/logs results/T20/opp/run_m2-6.log
