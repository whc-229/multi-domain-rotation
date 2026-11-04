# 图：新机队接入均匀轮转（第四节例子，tau = (2, 4)，两个岗位，6 台）。
# 数据由构造直接给出，与 code/theory/toy_connector.py 的穷举结果一致。
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams.update({"font.family": "serif", "font.serif": ["Liberation Serif"], "font.size": 8,
                     "mathtext.fontset": "stix", "axes.linewidth": 0.6})
TAU = (2, 4); P = sum(TAU); T = 12
C0, C1, CX, CW = "#1f4e79", "#9dc3e6", "#d9822b", "#bfbfbf"

def duties(b):          # 均匀轮转里偏移为 b 的机器人在 [0, T) 内的值班 (时刻, 相位)
    out = []
    for n in range(-2, 4):
        for i, o in enumerate((0, TAU[0])):
            t = b + o + n * P
            if 0 <= t < T: out.append((t, i))
    return out

fig, axes = plt.subplots(2, 1, figsize=(3.5, 2.75), sharex=True,
                         gridspec_kw={"hspace": 0.45})
def frame(ax, title, labels):
    ax.set_ylim(-0.6, 5.6); ax.invert_yaxis(); ax.set_xlim(-0.2, T)
    ax.set_yticks(range(6)); ax.set_yticklabels(labels)
    ax.set_xticks(range(T)); ax.tick_params(length=2, width=0.5)
    ax.set_title(title, fontsize=8, loc="left", pad=3)
    for x in range(T + 1): ax.axvline(x, color="#eeeeee", lw=0.4, zorder=0)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
def box(ax, t, r, c, hatch=None, ec="none"):
    ax.add_patch(Rectangle((t + 0.05, r - 0.32), 0.9, 0.64, fc=c, ec=ec, lw=0.8, hatch=hatch, zorder=2))

# (a) 均匀轮转
ax = axes[0]; frame(ax, "(a) Even rotation, 6 robots", [str(b) for b in range(6)])
for b in range(6):
    for t, i in duties(b):
        before = t < b
        box(ax, t, b, C0 if i == 0 else C1, ec=CX if before else "none")
# (b) 新机队
ax = axes[1]; frame(ax, "(b) New fleet entering (a)", [f"{b}\u2032" for b in range(6)])
for b in range(6):
    wend = {4: 0, 5: 1}.get(b, b)                  # 4'、5' 只在第一次值班前等待
    if wend > 0: ax.plot([0, wend], [b, b], color=CW, lw=1.0, ls=(0, (2, 1.5)), zorder=1)
    for t, i in duties(b):
        if t >= b: box(ax, t, b, C0 if i == 0 else C1)
for b in (4, 5):                                   # 4'、5'：先补 (a) 中 4、5 号在接入前的值班，再提前大修
    t0 = b - TAU[1]
    box(ax, t0, b, CX)
    ax.annotate("", xy=(b + 0.02, b), xytext=(t0 + 0.97, b),
                arrowprops=dict(arrowstyle="->", color="#555555", lw=0.7), zorder=3)
ax.set_xlabel("Time (duty lengths)", labelpad=1)

from matplotlib.lines import Line2D
h = [Rectangle((0, 0), 1, 1, fc=C0), Rectangle((0, 0), 1, 1, fc=C1),
     Rectangle((0, 0), 1, 1, fc="white", ec=CX, lw=0.8), Rectangle((0, 0), 1, 1, fc=CX),
     Line2D([0], [0], color="#555555", lw=0.7, marker=">", markersize=3, markevery=[1]),
     Line2D([0], [0], color=CW, lw=1.0, ls=(0, (2, 1.5)))]
fig.legend(h, ["phase 0", "phase 1", "before entry", "served first", "early overhaul", "waiting"],
           loc="lower center", ncol=3, fontsize=6, frameon=False, handlelength=1.2, columnspacing=0.8,
           bbox_to_anchor=(0.5, -0.01))
fig.subplots_adjust(left=0.09, right=0.99, top=0.94, bottom=0.24)
import sys; out = sys.argv[1] if len(sys.argv) > 1 else "fig_enter.pdf"
fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"), dpi=200)
