#!/bin/bash
# 消融对照：LP 引导派工 + curdef（与 lpgplan 只差提前大修计划）
for sp in "$@"; do
  [ -f plans/$sp.lpg.txt ] || continue
  for sd in 0 1 2; do
    if grep -q "\"spec\": \"$sp\", \"cfg\": \"lpgc\", \"seed\": $sd" stats4_lpgc.jsonl 2>/dev/null; then continue; fi
    r=$(cpp/rost4 es inst/$sp.txt certs 365 10 4 12 $sd 60000 600000 lpg=plans/$sp.lpg.txt curdef)
    python3 -c "import json,sys; d=json.loads(sys.argv[1]); d['spec']='$sp'; d['cfg']='lpgc'; d['seed']=$sd; print(json.dumps(d))" "$r" >> stats4_lpgc.jsonl
  done
done
