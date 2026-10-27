#!/bin/bash
# 只加 LP 提前大修计划（与 curdef 同维度、同随机数流，最干净的配对）
for sp in "$@"; do
  [ -s plans/$sp.f.txt ] || continue
  for sd in 0 1 2; do
    if grep -q "\"spec\": \"$sp\", \"cfg\": \"planf\", \"seed\": $sd" stats4_planf.jsonl 2>/dev/null; then continue; fi
    r=$(cpp/rost4 es inst/$sp.txt certs 365 10 4 12 $sd 60000 600000 planf=plans/$sp.f.txt)
    python3 -c "import json,sys; d=json.loads(sys.argv[1]); d['spec']='$sp'; d['cfg']='planf'; d['seed']=$sd; print(json.dumps(d))" "$r" >> stats4_planf.jsonl
  done
done
