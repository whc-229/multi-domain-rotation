"""任意两种方法的配对比较（同实例、同进化种子）：资源数（两者都长期有效时）、用时；Wilcoxon 符号秩检验。
用法：python3 compare.py A B 文件...   （A、B 形如 certs、finite3000、certs-gp）"""
import json, sys
from collections import defaultdict
from scipy.stats import wilcoxon
A, B = sys.argv[1], sys.argv[2]
rows = [json.loads(l) for fn in sys.argv[3:] for l in open(fn)]
rows = [d for d in rows if "error" not in d]
key = lambda d: f"{d['mode']}{d['H'] if d['mode'] == 'finite' else ''}" + ("-gp" if d.get("gp") else "")
valid = lambda d: d["R"] is not None and d["R"] >= 0 and (d["verify_cert"] == 1 or (d["verify_cert"] == -1 and d["verify_gaps"] == 0))
P = defaultdict(dict)
for d in rows: P[(d["inst"], d["seed"])][key(d)] = d
dR, ta, tb, va, vb, n = [], [], [], 0, 0, 0
for k, D in P.items():
    if A in D and B in D:
        a, b = D[A], D[B]; n += 1; va += valid(a); vb += valid(b); ta.append(a["sec"]); tb.append(b["sec"])
        if valid(a) and valid(b): dR.append(a["R"] - b["R"])
nz = [x for x in dR if x != 0]
p = wilcoxon(nz).pvalue if len(nz) >= 5 else float("nan")
print(f"{A} vs {B}：配对 {n}；长期有效 {va}/{n} 对 {vb}/{n}")
if dR:
    print(f"  两者都有效的 {len(dR)} 对：资源数差均值 {sum(dR)/len(dR):+.2f}（{A} 更少 {sum(x<0 for x in dR)}，相同 {sum(x==0 for x in dR)}，更多 {sum(x>0 for x in dR)}），Wilcoxon p={p:.3g}")
print(f"  平均用时 {sum(ta)/len(ta):.1f}s 对 {sum(tb)/len(tb):.1f}s")
