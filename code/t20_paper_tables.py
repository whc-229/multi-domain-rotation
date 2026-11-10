"""第五批（T20）论文表格与图。只读 results/T20/，写 manuscript/tables/ 与 manuscript/figures/。
口径（用户 11-05、11-09 定）：
- 下界 LB = ceil(max_H LP_fr^Rw(H))；差距 % = (资源数 - LB) / LB。
- 每个实例先对种子（GPHH 为 5 条规则）在有解的运行上取平均；
- 对比按配对算：每个对手只和本方法在两者都有解的实例上比（平均差距、胜/平/负、Wilcoxon，Holm 校正）；
- 另报求解率；Friedman 平均秩里无解按最差。
- 没跑完的 (方法, m) 不出数字，表里标 pending。
用法：python3 code/t20_paper_tables.py [仓库根目录]
"""
import json, glob, math, os, sys, statistics as st, collections
from scipy.stats import wilcoxon, rankdata

ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
R = lambda *p: os.path.join(ROOT, *p)
MS = [2, 3, 4, 6, 8, 10, 15, 20]
NAME = {'ours': 'C-RDD', 'pb': 'PB', 'gphh': 'GPHH', 'esdw': 'ES', 'rk': 'RK-EA', 'd5': 'DD-EA',
        'rddhand': 'RDD', 'rhmip': 'RH-MIP', 'conn_only': 'w/o search', 'noconn': 'w/o entering',
        'shorteval': 'w/o long-run eval.', 'ours_w0.25': '$0.25$', 'ours_w1': '$1$', 'ours_w2': '$2$'}
RUNS = {'ours': 5, 'pb': 1, 'gphh': 5, 'esdw': 5, 'rk': 5, 'd5': 5, 'rddhand': 5, 'rhmip': 1,
        'conn_only': 1, 'noconn': 3, 'shorteval': 3, 'ours_w0.25': 3, 'ours_w1': 3, 'ours_w2': 3}
OPP = ['pb', 'gphh', 'esdw', 'rk', 'd5', 'rddhand', 'rhmip']

# ---------- 读数据 ----------
lb, lbrow = {}, {}
for l in open(R('results/T20/lb_summary.jsonl')):
    d = json.loads(l); lbrow[d['inst']] = d
    lb[d['inst']] = math.ceil(max(v for v in d['frRw'].values() if v is not None) - 1e-6)
mof = lambda i: int(i[1:3])
runs = collections.defaultdict(lambda: collections.defaultdict(dict))  # method -> inst -> key -> R/None

def put(meth, inst, key, rec):
    ok = rec.get('status') in ('cert', 'uncert')
    runs[meth][inst][key] = rec['R'] if ok else None

for d in ('opp', 'ours'):
    for f in glob.glob(R('results/T20', d, 'jobs', '*.json')):
        meth, inst, key = os.path.basename(f)[:-5].split('__'); put(meth, inst, key, json.load(open(f)))
for f in glob.glob(R('results/T20/opp/rhmip/*.json')):
    rec = json.load(open(f)); put('rhmip', rec['inst'], 's0', rec)
# 装上 (e) 后的重跑：用 _e 的结果替换同名作业（规则见 TASKS 11-09）
for d in ('opp_e', 'ours_e'):
    for f in glob.glob(R('results/T20', d, 'jobs', '*.json')):
        meth, inst, key = os.path.basename(f)[:-5].split('__'); rec = json.load(open(f))
        if rec.get('e_fixed', -1) >= 0: put(meth, inst, key, rec)   # (e) 实际触发才用新结果（11-09 审核）

# q6 会用 (e) 重跑这些 (方法, m) 的一部分作业，跑完前按 pending 处理（q6 结束后删掉对应项）
HOLD = {('d5', 15), ('d5', 20)}

def complete(meth, m):
    if (meth, m) in HOLD: return False
    ins = [i for i in lb if mof(i) == m]
    return all(len(runs[meth].get(i, {})) >= RUNS[meth] for i in ins)

def imean(meth, i):
    v = [x for x in runs[meth].get(i, {}).values() if x is not None]
    return st.mean(v) if v else None

gap = lambda x, i: 100.0 * (x - lb[i]) / lb[i]
GEN = {}   # 生成的表格行，最后写进 06_experiments.tex 的标记之间（Overleaf 只需同步这一个文件）
PEND = r'\textit{pending}'

# ---------- 表：下界与启动代价 ----------
rows = []
for m in MS:
    ins = [i for i in lb if mof(i) == m]
    rate = [lbrow[i]['R_rate_g'] for i in ins]
    fr = [max(lbrow[i]['fr'].values()) for i in ins]
    rw = [max(lbrow[i]['frRw'].values()) for i in ins]
    h0 = st.median(lbrow[i]['frRw_H0'] for i in ins) if 'frRw_H0' in lbrow[ins[0]] else float('nan')
    c_fr = st.mean(100 * (a - b) / b for a, b in zip(fr, rate))
    c_rw = st.mean(100 * (a - b) / b for a, b in zip(rw, rate))
    rows.append(f'{m} & {st.mean(rate):.1f} & {st.mean(lb[i] for i in ins):.1f} & {c_fr:.1f} & {c_rw:.1f} & {h0:g} \\\\')
GEN['tab_bounds'] = ('\n'.join(rows) + '\n')

# ---------- 表：对比（配对） ----------
pv = {}
for o in OPP:
    d = []
    for i in lb:
        if not (complete(o, mof(i)) and complete('ours', mof(i))): continue
        a, b = imean('ours', i), imean(o, i)
        if a is not None and b is not None: d.append(b - a)
    pv[o] = wilcoxon(d).pvalue if len(d) >= 6 and any(x != 0 for x in d) else float('nan')
order = sorted(OPP, key=lambda o: pv[o]); holm = {}; run_max = 0
for k, o in enumerate(order):
    run_max = max(run_max, min(1.0, (len(OPP) - k) * pv[o])); holm[o] = run_max
rows = []
for m in MS:
    ins = [i for i in lb if mof(i) == m]
    if not complete('ours', m):
        rows.append(f'{m} & {PEND} & ' + ' & '.join([PEND] * len(OPP)) + r' \\'); continue
    cells = [f"{st.mean(gap(imean('ours', i), i) for i in ins):.1f}"]
    for o in OPP:
        if not complete(o, m): cells.append(PEND); continue
        both = [i for i in ins if imean(o, i) is not None]
        if not both: cells.append('--'); continue
        dg = st.mean(gap(imean(o, i), i) - gap(imean('ours', i), i) for i in both)
        w = sum(imean('ours', i) < imean(o, i) - 1e-9 for i in both)
        l = sum(imean('ours', i) > imean(o, i) + 1e-9 for i in both)
        cells.append(f'{dg:+.1f} ({w}/{len(both)-w-l}/{l})')
    rows.append(f'{m} & ' + ' & '.join(cells) + r' \\')
rows.append(r'\midrule')
rows.append(r'$p$ (Holm) & & ' + ' & '.join('--' if math.isnan(holm[o]) else (r'$<10^{-4}$' if holm[o] < 1e-4 else f'{holm[o]:.1e}') for o in OPP) + r' \\')
GEN['tab_compare'] = ('\n'.join(rows) + '\n')

# ---------- 表：求解率 ----------
rows = []
for meth in ['ours'] + OPP:
    cells = []
    for m in MS:
        if not complete(meth, m): cells.append('$\\cdot$'); continue
        v = [x for i in lb if mof(i) == m for x in runs[meth][i].values()]
        cells.append(f'{100*sum(x is not None for x in v)/len(v):.0f}')
    rows.append(f'{NAME[meth]} & ' + ' & '.join(cells) + r' \\')
GEN['tab_solve'] = ('\n'.join(rows) + '\n')

# ---------- Friedman 平均秩（无解按最差；只用全部方法都跑完的 m） ----------
meths = ['ours'] + OPP; ranks = collections.defaultdict(list); nrank = 0
for i in lb:
    if not all(complete(x, mof(i)) for x in meths): continue
    vals = [imean(x, i) if imean(x, i) is not None else 1e9 for x in meths]
    for x, r in zip(meths, rankdata(vals)): ranks[x].append(r)
    nrank += 1
fr_line = ' & '.join(f'{st.mean(ranks[x]):.2f}' if ranks[x] else '--' for x in meths)
GEN['tab_rank'] = (
    f'% Friedman 平均秩，实例数 {nrank}\n' + ' & '.join(NAME[x] for x in meths) + r' \\' + '\n' + fr_line + r' \\' + '\n')

# ---------- 表：消融 ----------
rows = []
def abl(meth, label):
    out = []
    for m in MS:
        if meth not in runs or not any(mof(i) == m for i in runs[meth]): continue
        if not (complete(meth, m) and complete('ours', m)): out.append((m, None)); continue
        ins = [i for i in lb if mof(i) == m and imean(meth, i) is not None]
        dg = st.mean(gap(imean(meth, i), i) - gap(imean('ours', i), i) for i in ins)
        w = sum(imean('ours', i) < imean(meth, i) - 1e-9 for i in ins)
        l = sum(imean('ours', i) > imean(meth, i) + 1e-9 for i in ins)
        out.append((m, f'{dg:+.1f} ({w}/{len(ins)-w-l}/{l})'))
    return out
for meth in ['rddhand', 'noconn', 'conn_only', 'shorteval']:
    for m, c in abl(meth, NAME[meth]):
        rows.append(f'{NAME[meth] if meth != "rddhand" else "w/o relaxation rule"} & {m} & {c if c else PEND} \\\\')
    rows.append(r'\midrule')
GEN['tab_ablation'] = ('\n'.join(rows[:-1]) + '\n')

# ---------- 表：偏离窗口长度 ----------
rows = []
for meth in ['ours_w0.25', 'ours', 'ours_w1', 'ours_w2']:
    ins = [i for i in lb if mof(i) in (4, 6)]
    if not all(complete(meth, m) for m in (4, 6)): continue
    g = st.mean(gap(imean(meth, i), i) for i in ins)
    if meth == 'ours': rows.append(f'$0.5$ & {g:.1f} & -- \\\\'); continue
    w = sum(imean(meth, i) < imean('ours', i) - 1e-9 for i in ins); l = sum(imean(meth, i) > imean('ours', i) + 1e-9 for i in ins)
    rows.append(f'{NAME[meth]} & {g:.1f} & {w}/{len(ins)-w-l}/{l} \\\\')
GEN['tab_window'] = ('\n'.join(rows) + '\n')

# ---------- 长期可行评价消融：只看 720 h 选出的解，长期复核缺岗的比例与补救量 ----------
se = [json.load(open(f)) for f in glob.glob(R('results/T20/ours/jobs/shorteval__*.json'))]
short = [r for r in se if r.get('long_feasible') == 0]
add = [r['rescue_long'].get('rescue_add', 0) for r in short if r.get('rescue_long')]
nfail = sum(r.get('status') not in ('cert', 'uncert') for r in se)
GEN['stat_shorteval'] = (
    f'% shorteval 作业 {len(se)} 个；只看 720 h 选出的最终解长期复核缺岗 {len(short)} 个；补救加 {min(add) if add else 0}–{max(add) if add else 0} 个；补救后仍不可行 {nfail} 个\n'
    f'\\newcommand{{\\SEn}}{{{len(se)}}}\\newcommand{{\\SEshort}}{{{len(short)}}}\\newcommand{{\\SEpct}}{{{round(100*len(short)/max(1,len(se)))}}}'
    f'\\newcommand{{\\SEaddmin}}{{{min(add) if add else 0}}}\\newcommand{{\\SEaddmax}}{{{max(add) if add else 0}}}\\newcommand{{\\SEfail}}{{{nfail}}}\n')

# ---------- 图：启动下界随 H 收敛 ----------
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Liberation Serif'], 'font.size': 8, 'mathtext.fontset': 'stix', 'axes.linewidth': 0.6})
H = ['0.5', '1', '1.5', '2', '3', '4']; hx = [float(h) for h in H]
fig, ax = plt.subplots(figsize=(3.5, 2.1))
cols = plt.cm.Blues([0.35 + 0.65 * k / (len(MS) - 1) for k in range(len(MS))])
for c, m in zip(cols, MS):
    ins = [i for i in lb if mof(i) == m]
    y = [100 * st.median(lbrow[i]['frRw'][h] / lbrow[i]['frRw']['4'] for i in ins) for h in H]
    ax.plot(hx, y, color=c, lw=1.2, marker='o', ms=2.5, label=f'$m={m}$')
ax.set_xlabel(r'Horizon $H$ ($\times S_{\max}$)'); ax.set_ylabel(r'$R_0(H)$ / $R_0(4S_{\max})$ (%)')
ax.set_xticks(hx); ax.grid(axis='y', color='#e5e5e5', lw=0.5); [ax.spines[s].set_visible(False) for s in ('top', 'right')]
ax.legend(ncol=2, fontsize=6, frameon=False, loc='lower right', handlelength=1.5)
fig.tight_layout(); fig.savefig(R('manuscript/figures/fig_lpconv.pdf'), metadata={'CreationDate': None}); fig.savefig(R('manuscript/figures/fig_lpconv.png'), dpi=200)

# ---------- 图：区间收紧（每个 m 的平均，以 LB 为 0） ----------
pre = {}
for f in glob.glob(R('results/T20/ours/pre/*.json')):
    p = json.load(open(f)); pre[p['inst']] = p
ms_ok = [m for m in MS if complete('ours', m) and complete('conn_only', m)]
series = collections.OrderedDict([('Rate bound', []), ('Rotation $\\sigma$ (not fresh)', []), ('One entering construction', []),
                                  ('Best entering construction', []), ('C-RDD', [])])
for m in ms_ok:
    ins = [i for i in lb if mof(i) == m and i in pre]
    series['Rate bound'].append(st.mean(gap(math.ceil(lbrow[i]['R_rate_g'] - 1e-6), i) for i in ins))
    series['Rotation $\\sigma$ (not fresh)'].append(st.mean(gap(pre[i]['sigma_R'], i) for i in ins))
    one = []
    for i in ins:
        th = [t for t in pre[i]['thetas'] if t['mode'] == 'first' and t.get('feasible')]
        th = min(th, key=lambda t: t['theta']) if th else None
        if th: one.append(gap(th['R'], i))
    series['One entering construction'].append(st.mean(one))
    series['Best entering construction'].append(st.mean(gap(imean('conn_only', i), i) for i in ins))
    series['C-RDD'].append(st.mean(gap(imean('ours', i), i) for i in ins))
if ms_ok:
    fig, ax = plt.subplots(figsize=(3.5, 2.3)); xs = range(len(ms_ok))
    sty = {'Rate bound': ('#7f7f7f', '--'), 'Rotation $\\sigma$ (not fresh)': ('#9dc3e6', ':'), 'One entering construction': ('#d9822b', '-'),
           'Best entering construction': ('#9dc3e6', '-'), 'C-RDD': ('#1f4e79', '-')}
    ax.axhline(0, color='black', lw=0.6); ax.text(len(ms_ok) - 1, 0.6, 'start-up bound', ha='right', fontsize=6)
    for k, v in series.items():
        c, ls = sty[k]; ax.plot(xs, v, color=c, ls=ls, lw=1.2, marker='o', ms=2.5, label=k)
    ax.set_xticks(list(xs)); ax.set_xticklabels([str(m) for m in ms_ok]); ax.set_xlabel('Number of sites $m$')
    ax.set_ylabel('Above the start-up bound (%)'); ax.grid(axis='y', color='#e5e5e5', lw=0.5)
    [ax.spines[s].set_visible(False) for s in ('top', 'right')]
    ax.legend(fontsize=6, frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.28), handlelength=1.8)
    fig.tight_layout(); fig.savefig(R('manuscript/figures/fig_bracket.pdf'), metadata={'CreationDate': None}); fig.savefig(R('manuscript/figures/fig_bracket.png'), dpi=200)
print('bracket series:', {k: [round(x, 1) for x in v] for k, v in series.items()})
print('complete m (ours):', [m for m in MS if complete('ours', m)], '; Friedman instances:', nrank, '; bracket m:', ms_ok)
print('holm', {o: holm[o] for o in OPP})

# ---------- 把生成的内容写进 06_experiments.tex 的 "% BEGIN 名" 与 "% END 名" 之间 ----------
import re as _re
sec = R('manuscript/sections/06_experiments.tex'); txt = open(sec).read()
for k, v in GEN.items():
    pat = _re.compile(r'([ \t]*)% BEGIN ' + k + r'\n.*?% END ' + k, _re.S)
    if not pat.search(txt): print('缺少标记', k); continue
    txt = pat.sub(lambda mo: f"{mo.group(1)}% BEGIN {k}\n{v.rstrip()}\n{mo.group(1)}% END {k}", txt)
open(sec, 'w').write(txt)
print('已写入', sec, sorted(GEN))
