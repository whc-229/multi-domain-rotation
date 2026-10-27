"""按站点数分组：各方法长期有效数、平均用时；certs 对其他方法的配对资源数差（两者都长期有效时）。
用法：python3 per_m.py 文件..."""
import json, sys
from collections import defaultdict
rows = [json.loads(l) for fn in sys.argv[1:] for l in open(fn)]
rows = [d for d in rows if "error" not in d]
key = lambda d: f"{d['mode']}{d['H'] if d['mode'] == 'finite' else ''}" + ("-gp" if d.get("gp") else "")
valid = lambda d: d["R"] is not None and d["R"] >= 0 and (d["verify_cert"] == 1 or (d["verify_cert"] == -1 and d["verify_gaps"] == 0))
mof = lambda d: int(d["inst"].split("-")[0][1:])
by = defaultdict(lambda: defaultdict(list))
P = defaultdict(dict)
for d in rows:
    by[mof(d)][key(d)].append(d); P[(d["inst"], d["seed"])][key(d)] = d
for m in sorted(by):
    print(f"== m={m}")
    for k in sorted(by[m]):
        L = by[m][k]
        print(f"  {k:12s} 有效 {sum(map(valid, L))}/{len(L)}  证明 {sum(d['verify_cert'] == 1 for d in L)}  "
              f"平均用时 {sum(d['sec'] for d in L)/len(L):.1f}s  R={[d['R'] for d in L]}")
    for B in ["finite365", "finite3000"]:
        dR = []
        for (inst, sd), D in P.items():
            if mof(D[next(iter(D))]) != m or "certs" not in D or B not in D: continue
            if valid(D["certs"]) and valid(D[B]): dR.append(D["certs"]["R"] - D[B]["R"])
        if dR: print(f"  certs − {B}：{len(dR)} 对，差 {dR}，均值 {sum(dR)/len(dR):+.2f}")
    # 每个实例上长期有效解的最好资源数
    best = defaultdict(dict)
    for k, L in by[m].items():
        for d in L:
            if valid(d): best[d["inst"]][k] = min(best[d["inst"]].get(k, 10**9), d["R"])
    for inst in sorted(best): print(f"  {inst} 各方法有效解最好值：{dict(sorted(best[inst].items()))}")
