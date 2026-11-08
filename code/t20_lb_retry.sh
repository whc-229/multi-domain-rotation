#!/bin/bash
# T20-1 超时 LP 的重跑（11-09 补排）：37 个 LP 在 1 小时限时内没解完（lp_status = "Time limit reached"，没有下界），
# 按 TASKS T20-1 与 code/t20_lb_job.py 已定的做法只换求解设置重跑一次（retry：行列随机重排 + 关 crossover），结果写 <名>__retry1.json，原结果不动。
# 作业表 results/T20/lb/jobs_retry.tsv。等 E2 量测（PID 775923，命令行要对得上）结束、腾出 3 个核再跑，3 个并行，不动总调度。
# 断了重跑同一条命令即可（已完成的跳过）。用法：setsid nohup bash code/t20_lb_retry.sh >> results/T20/lb/retry_launcher.log 2>&1 < /dev/null &
cd "$(dirname "$0")/.." || exit 1
while kill -0 775923 2>/dev/null && grep -q "t20_check_e.py" /proc/775923/cmdline 2>/dev/null; do sleep 120; done
bash code/t20_run.sh results/T20/lb/jobs_retry.tsv 3 results/T20/lb/logs_retry results/T20/lb/run_retry.log
