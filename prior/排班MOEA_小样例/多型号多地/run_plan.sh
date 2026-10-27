#!/bin/bash
# LP 计划 + 进化派工：对每个实例先由全新启动 LP（带工具、允许提前大修）生成计划，再与 curdef（同一随机数流）配对比较。
for sp in "$@"; do
  [ -f plans/$sp.txt ] || python3 make_plan.py $sp plans/$sp.txt 4 >> plans/make_plan.log
  for sd in 0 1 2; do
    for cfg in plan; do
      if grep -q "\"spec\": \"$sp\", \"cfg\": \"$cfg\", \"seed\": $sd" stats4_plan.jsonl 2>/dev/null; then continue; fi
      if [ $cfg = plan ]; then opt="plan=plans/$sp.txt"; else opt="curdef"; fi
      r=$(cpp/rost4 es inst/$sp.txt certs 365 10 4 12 $sd 60000 600000 $opt)
      python3 -c "import json,sys; d=json.loads(sys.argv[1]); d['spec']='$sp'; d['cfg']='$cfg'; d['seed']=$sd; print(json.dumps(d))" "$r" >> stats4_plan.jsonl
    done
  done
done
