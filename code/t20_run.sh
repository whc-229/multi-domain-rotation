#!/bin/bash
# T20 调度：作业表每行 "<超时秒> job <阶段> <方法> <实例> <种子>" 或 "<超时秒> lb <实例> <种类> <倍数> [retry]"，P 个并行；
# 已完成的作业自动跳过（作业脚本自己判断），断了重跑同一条命令即可。
# 每个作业的输出写 <日志目录>/<作业名>.log，进度（时刻、作业、退出码）追加到 <进度文件>。
# 用法：bash code/t20_run.sh <作业表> <并行数> <日志目录> <进度文件>
cd "$(dirname "$0")/.." || exit 1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export T20_PY=/home/hk/miniconda3/envs/ana/bin/python3 T20_LOG="$3"
mkdir -p "$3"
xargs -P "$2" -L 1 bash -c '
  to="$0"; kind="$1"; shift 1
  if [ "$kind" = job ]; then name="$2__$3__s$4__$1"; cmd=(code/t20_job.py "$@"); else name="lb__$1__$2__m$3${4:+__$4}"; cmd=(code/t20_lb_job.py "$@"); fi
  timeout "$to" "$T20_PY" "${cmd[@]}" > "$T20_LOG/$name.log" 2>&1; rc=$?
  echo "$(date "+%m-%d %H:%M:%S") $kind $* 退出码 $rc"' < "$1" >> "$4"
echo "$(date "+%m-%d %H:%M:%S") 全部结束" >> "$4"
