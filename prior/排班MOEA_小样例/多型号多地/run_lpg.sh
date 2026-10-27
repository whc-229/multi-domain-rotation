#!/bin/bash
# LP 引导：lpg = LP 机型-站点份额作为派工特征 + LP 台数作为定台数起点；lpgplan = 再加按比例的 LP 提前大修计划
for sp in "$@"; do
  if [ ! -f plans/$sp.lpg.txt ]; then python3 make_planf.py $sp plans/$sp.f.txt 4 >> plans/make_planf.log; fi
  for sd in 0 1 2; do
    for cfg in lpg lpgplan; do
      if grep -q "\"spec\": \"$sp\", \"cfg\": \"$cfg\", \"seed\": $sd" stats4_lpg.jsonl 2>/dev/null; then continue; fi
      if [ $cfg = lpg ]; then opt="lpg=plans/$sp.lpg.txt"; else opt="lpg=plans/$sp.lpg.txt planf=plans/$sp.f.txt"; fi
      r=$(cpp/rost4 es inst/$sp.txt certs 365 10 4 12 $sd 60000 600000 $opt)
      python3 -c "import json,sys; d=json.loads(sys.argv[1]); d['spec']='$sp'; d['cfg']='$cfg'; d['seed']=$sd; print(json.dumps(d))" "$r" >> stats4_lpg.jsonl
    done
  done
done
