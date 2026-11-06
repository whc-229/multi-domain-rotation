#!/bin/bash
# T20 无人值守总调度（11-06 晚，TASKS"无人值守安排"）。按顺序跑：
#   1. T20-1 剩余的 LP（若 11-06 起的 T20-1 进程 351119 还在跑，先等它结束，再把同一张表过一遍，补上没跑完的）；
#   2. q2：T20-2 m ≤ 6 的全部对手；
#   3. q3：T20-2 m ≥ 8 的确定性方法（PB、GPHH、RH-MIP）与只用接入构造；
#   4. q4：T20-3（先 m ≤ 6 再 m ≥ 8）；
#   5. q5：T20-2 m ≥ 8 的搜索型对手；
#   6. q6：(e) 装上后，把装上前搜索型方法定不出台数、记无解的作业排到队尾重跑（code/t20_requeue_e.py 现生成；没装 (e) 时为空）。
# 作业表见 code/t20_make_queue.py。已完成的作业由作业脚本自己跳过，断了重跑同一条命令即可续跑。
# 同一时间只能有一个总调度：用 results/T20/master.lock 加锁（子进程继承锁，总调度被杀但作业还在跑时也不会重复启动）。
# 全部跑完写 results/T20/master.done，看门狗（code/t20_watchdog.sh）见到它就不再拉起。
# 用法：setsid nohup bash code/t20_master.sh >> results/T20/master.log 2>&1 < /dev/null &
cd "$(dirname "$0")/.." || exit 1
exec 9> results/T20/master.lock
flock -n 9 || { echo "$(date "+%m-%d %H:%M:%S") 已有总调度在跑，退出"; exit 0; }
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/home/hk/miniconda3/envs/ana/bin/python3
log() { echo "$(date "+%m-%d %H:%M:%S") $*"; }
log "总调度启动 PID $$"
echo $$ > results/T20/master.pid
# 1. T20-1：等旧进程（只认记录的 PID，且命令行要对得上）
while kill -0 351119 2>/dev/null && grep -q "t20_run.sh" /proc/351119/cmdline 2>/dev/null; do sleep 120; done
log "1. T20-1"; bash code/t20_run.sh results/T20/lb/jobs.tsv 9 results/T20/lb/logs results/T20/lb/run.log
log "2. q2"; bash code/t20_run.sh results/T20/queue/q2.tsv 12 results/T20/opp/logs results/T20/queue/run_q2.log
log "3. q3"; bash code/t20_run.sh results/T20/queue/q3.tsv 12 results/T20/queue/logs_q3 results/T20/queue/run_q3.log
log "4. q4"; bash code/t20_run.sh results/T20/queue/q4.tsv 12 results/T20/ours/logs results/T20/queue/run_q4.log
log "5. q5"; bash code/t20_run.sh results/T20/queue/q5.tsv 12 results/T20/opp/logs results/T20/queue/run_q5.log
log "6. q6"; $PY code/t20_requeue_e.py results/T20/queue/q6.tsv
if [ -s results/T20/queue/q6.tsv ]; then bash code/t20_run.sh results/T20/queue/q6.tsv 12 results/T20/queue/logs_q6 results/T20/queue/run_q6.log; fi
# (e) 若在 q6 生成之后才装上，下一轮看门狗拉起时会重新生成 q6；只有 (e) 已装上（或确定不装）且 q6 跑完才算结束
if $PY code/t20_requeue_e.py --final-check; then date > results/T20/master.done; log "全部结束"; else log "q6 之后仍有待重跑的作业（或 (e) 未定），留给看门狗下一轮"; fi
