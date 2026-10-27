"""统计实验汇总：每种适应度的长期有效率、资源数（只在有效解上比较）、用时；cert 与 finite3000 的配对检验。
长期有效 = 60 万格复核中被证明永远可行（verify_cert=1），或 60 万格内无缺岗且未判出循环（记为"未证但无缺岗"）。"""
import json, sys, math
from collections import defaultdict
from scipy.stats import wilcoxon
rows = []
for fn in sys.argv[1:]:
    for l in open(fn):
        d = json.loads(l)
        if "error" in d: continue
        rows.append(d)
key = lambda d: f"{d['mode']}{d['H'] if d['mode'] == 'finite' else ''}" + ("-gp" if d.get("gp") else "")
def valid(d): return d["R"] is not None and d["R"] >= 0 and (d["verify_cert"] == 1 or (d["verify_cert"] == -1 and d["verify_gaps"] == 0))
def proven(d): return d["verify_cert"] == 1
by = defaultdict(list)
for d in rows: by[key(d)].append(d)
print("== 各适应度总体 ==")
for k, L in sorted(by.items()):
    n = len(L); nv = sum(valid(d) for d in L); npf = sum(proven(d) for d in L)
    print(f"{k:12s} 运行数 {n:3d}  长期有效 {nv:3d} ({100*nv/n:.0f}%)  其中严格证明 {npf:3d}  平均用时 {sum(d['sec'] for d in L)/n:.1f}s")
# 按站点数
print("== 按站点数：finite365 失效率；cert / finite3000 有效率 ==")
bym = defaultdict(lambda: defaultdict(list))
for d in rows: bym[int(d["inst"].split("-")[0][1:]) if d["inst"].startswith("m") else 1][key(d)].append(d)
for m in sorted(bym):
    parts = []
    for k in sorted(bym[m]):
        L = bym[m][k]; parts.append(f"{k}: 有效 {sum(valid(d) for d in L)}/{len(L)}")
    print(f"m={m}: " + "；".join(parts))
# 配对：同一 (实例, 种子)，cert vs finite3000（只在两者都有效时比较资源数），用时比较
pair = defaultdict(dict)
for d in rows: pair[(d["inst"], d["seed"], d.get("gp", 0))][key(d).replace("-gp", "")] = d
dR, ts, tf, both = [], [], [], 0
for (inst, sd, gp), D in pair.items():
    if "cert" in D and "finite3000" in D:
        a, b = D["cert"], D["finite3000"]
        ts.append(a["sec"]); tf.append(b["sec"])
        if valid(a) and valid(b): dR.append(a["R"] - b["R"]); both += 1
if dR:
    nz = [x for x in dR if x != 0]
    p = wilcoxon(nz).pvalue if len(nz) >= 5 else float("nan")
    print(f"== cert − finite3000 资源数差（两者都有效的 {both} 对）：均值 {sum(dR)/len(dR):+.2f}，"
          f"cert 更少 {sum(x<0 for x in dR)}，相同 {sum(x==0 for x in dR)}，更多 {sum(x>0 for x in dR)}，Wilcoxon p={p:.3g}")
    print(f"== 用时：cert 平均 {sum(ts)/len(ts):.1f}s，finite3000 平均 {sum(tf)/len(tf):.1f}s，比值 {sum(tf)/sum(ts):.2f}")
# finite365 失效时，cert 在同资源数下是否有效
same, more, less, tot = 0, 0, 0, 0
for (inst, sd, gp), D in pair.items():
    if "finite365" in D and "cert" in D and not valid(D["finite365"]) and valid(D["cert"]):
        tot += 1; x = D["cert"]["R"] - D["finite365"]["R"]
        same += (x == 0); more += (x > 0); less += (x < 0)
print(f"== finite365 失效的 {tot} 次运行中，同种子 cert 的资源数：相同 {same}、更多 {more}、更少 {less}")
