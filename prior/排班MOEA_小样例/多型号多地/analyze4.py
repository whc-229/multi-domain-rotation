"""提前大修实验汇总：同一 (实例, 种子) 上，基线 certs（rost3 结果文件）与 rost4 各配置的资源数、长期有效性、环内提前大修比例。
用法：python3 analyze4.py <rost4 结果.jsonl> <基线文件...>"""
import json, sys
from collections import defaultdict
from scipy.stats import wilcoxon
valid = lambda d: d.get("R") is not None and d["R"] >= 0 and (d["verify_cert"] == 1 or (d["verify_cert"] == -1 and d["verify_gaps"] == 0))
P = defaultdict(dict)
for l in open(sys.argv[1]):
    d = json.loads(l)
    if "error" in d: continue
    P[(d["inst"], d["seed"])][d["cfg"]] = d
for fn in sys.argv[2:]:
    for l in open(fn):
        d = json.loads(l)
        if "error" in d or d.get("gp") or d["mode"] != "certs": continue
        if (d["inst"], d["seed"]) in P: P[(d["inst"], d["seed"])]["base"] = d
cfgs = sorted({c for D in P.values() for c in D})
print("配置：", cfgs)
for c in cfgs:
    L = [D[c] for D in P.values() if c in D]
    rr = [d["cyc_resets"] / d["cyc_duties"] for d in L if d.get("cyc_duties", -1) > 0]
    print(f"{c:10s} 运行 {len(L):3d}  有效 {sum(map(valid, L)):3d}  证明 {sum(d['verify_cert'] == 1 for d in L):3d}  平均 R {sum(d['R'] for d in L)/len(L):.2f}  "
          f"平均用时 {sum(d['sec'] for d in L)/len(L):.0f}s" + (f"  环内提前大修/值岗 均值 {sum(rr)/len(rr):.3f}，为 0 的 {sum(x == 0 for x in rr)}/{len(rr)}" if rr else ""))
def pair(a, b):
    dR = [D[a]["R"] - D[b]["R"] for D in P.values() if a in D and b in D and valid(D[a]) and valid(D[b])]
    if not dR: return
    nz = [x for x in dR if x != 0]
    p = wilcoxon(nz).pvalue if len(nz) >= 5 else float("nan")
    print(f"{a} − {b}：{len(dR)} 对，均值 {sum(dR)/len(dR):+.2f}，少 {sum(x<0 for x in dR)} 同 {sum(x==0 for x in dR)} 多 {sum(x>0 for x in dR)}，p={p:.3g}")
for a, b in [("curdef", "base"), ("pad", "base"), ("pad", "curdef"), ("reset", "base"), ("reset", "curdef"), ("reset", "pad"), ("hard", "pad"), ("hard", "reset")]:
    pair(a, b)
