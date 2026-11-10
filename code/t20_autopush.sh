#!/bin/bash
# T20 每完成一个队列推一次（11-10 审核"之后每完成一个队列推一次"）。由 crontab 每 30 分钟调用，不改正在跑的总调度。
# 判定"队列跑完"：run_<q>.log 里有 t20_run.sh 结束时才写的"全部结束"（不按行数判：总调度被看门狗重启后，已完成的作业会再记一行）；
# 总调度全部结束看 master.done。推送时在提交说明里写上没有结果 JSON 的作业数（应为 0，不为 0 要人工查）。
# 每个事件只推一次：推成功后写标记 results/T20/queue/pushed_<事件>（不进仓库）。
# 只提交 results/T20 下的结果与运行日志（作业 JSON、队列表与进度日志、总调度日志），不碰别的改动。
# 与手动的 git 操作互斥：flock results/T20/git.lock。远端有新提交时先合并（不 stash，不 rebase）；合并冲突就放弃，留给人工处理。
# 用法：bash code/t20_autopush.sh（crontab：*/30 * * * *）
cd "$(dirname "$0")/.." || exit 1
Q=results/T20/queue
log() { echo "$(date "+%m-%d %H:%M:%S") autopush: $*"; }

# 作业表 $1 里没有结果 JSON 的作业数（行格式 "<超时> job <阶段> <方法> <实例> <种子>"）
n_missing() {
    local n=0 to kind st me in se
    while read -r to kind st me in se; do
        [ "$kind" = job ] && [ ! -f "results/T20/$st/jobs/${me}__${in}__s${se}.json" ] && n=$((n + 1))
    done < "$1"
    echo $n
}

# 打印已经发生、还没推过的事件名（q5、q6、done）
events() {
    for q in q5 q6; do
        [ -s "$Q/$q.tsv" ] && [ -f "$Q/run_$q.log" ] && [ ! -f "$Q/pushed_$q" ] || continue
        grep -q "全部结束" "$Q/run_$q.log" && echo "$q"
    done
    [ -f results/T20/master.done ] && [ ! -f "$Q/pushed_done" ] && echo done
}

# 提交并推送一个事件；成功返回 0
push_event() {
    local ev="$1" miss=""
    [ "$ev" != done ] && miss="；没有结果的作业 $(n_missing "$Q/$ev.tsv") 个"
    local p paths=()
    for p in results/T20/opp results/T20/ours results/T20/opp_e results/T20/ours_e results/T20/queue results/T20/master.log results/T20/e_status.json; do
        [ -e "$p" ] && paths+=("$p")   # 不存在的路径会让整条 git add 失败
    done
    git add -- "${paths[@]}" || { log "$ev：git add 失败"; return 1; }
    if git diff --cached --quiet; then log "$ev：没有新结果"; return 0; fi
    git commit -q -m "T20 自动推送：$ev 跑完的结果与运行日志$miss

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" || { log "$ev：commit 失败"; return 1; }
    git fetch -q origin || { log "$ev：fetch 失败"; return 1; }
    if ! git merge -q --no-edit origin/main; then git merge --abort; log "$ev：合并冲突，已放弃合并，留给人工"; return 1; fi
    git push -q origin HEAD:main || { log "$ev：push 失败"; return 1; }
    log "$ev：已推送 $(git rev-parse --short HEAD)$miss"
}

exec 8> results/T20/git.lock
flock -n 8 || { log "有别的 git 操作在进行，下一轮再试"; exit 0; }
for ev in $(events); do
    push_event "$ev" && touch "$Q/pushed_$ev"
done
