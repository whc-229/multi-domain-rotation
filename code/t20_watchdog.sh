#!/bin/bash
# T20 总调度的看门狗（crontab 每 30 分钟跑一次）：总调度没在跑（锁没人拿着）且没写 master.done，就用同一条命令拉起；
# 已完成的作业自动跳过。机器重启后也由它拉起。只做检查和启动，不结束任何进程。
cd /home/hk/多域轮转 || exit 1
[ -f results/T20/master.done ] && exit 0
if flock -n results/T20/master.lock true; then
  echo "$(date "+%m-%d %H:%M:%S") 总调度不在，拉起"
  setsid nohup bash code/t20_master.sh >> results/T20/master.log 2>&1 < /dev/null &
fi
