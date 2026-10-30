"""T13：各对手方法在 rost5 上的设置（规则选项、搜索命令）。所有方法共用同一个外壳：
rost5 仿真、同一套定台数（size_rule）与长期可行判定（严格证书或 60 万格复核）、全新启动、只在启动窗口 [0, W) 内提前大修（W = 2·S_max）。
对手只用自己的信息，不用我们的 LP 规则与计划；ES-LP 例外（它本来就用 LP 引导与计划），单列。

方法：
- fifo：先到先派。可用、允许的机器人里就绪最早的（rost5 es 的基础个体 0：w = [1,0,0,0,0,0,1,0]）；工具按保养时长；不提前大修。
- pb：相位平衡（T8 的 hand，rost5 es 的基础个体 1）；不提前大修。
- esdw：逐实例 ES 调派工权重（curdef，8 维），用定台数下降评价（rost5 fd 表示 w）；不提前大修。
- eslp：ES 调权重 + LP 引导（lpz，权重 w13）+ 启动窗口内的 LP_fr^Rw 提前大修计划（planf，严格窗口 planW）；fd 表示 w。
- gphh：GP 树打分（T14 训练出的规则），直接定台数，不逐实例搜索；不提前大修。
- rk：原排班软件的随机键：窗口内每个决策一个基因，取可行候选（按编号）的第 ⌊x·n⌋ 个；窗口后取编号最小的可行候选；fd 表示 rk；不提前大修。
- d5：上一篇的 dec5：候选按出勤少、相位低、编号排序，窗口内按基因取，缺岗时回溯 12 个决策修复；窗口后取第一个；fd 表示 d5；不提前大修。
- rddhand：上一篇的 RDD：默认规则 PB，备选表与 T7 相同（含"本班后提前大修"，只在窗口内），rost5 rddfd；不含 LP 的任何信息。
- ours：我们的方法（T12 的设置）：fresh_rlp_w（严格窗口、后备补足时域）+ RDD-FD + resetopt。
所有搜索型方法（esdw、eslp、rk、d5、rddhand、ours）用同一个 (4 + 12) 框架与定台数下降评价（T_eval = 6·S_max，互换 0.5，stall 30），只有表示不同。"""
import os, sys, math, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train  # noqa: E402

RULES6 = os.path.join(ROOT, "results", "T6", "rules")
RULES13 = os.path.join(ROOT, "results", "T13", "rules")
W_FIFO = [1, 0, 0, 0, 0, 0, 1, 0]
W_PB = [0.1, 0, 10, 1, 0, 0, 1, 0.1]
W_TOOL = [0, 0, 0, 0, 0, 0, 1, 0]      # 只有工具打分（保养时长短的先用），给 rk / d5 / gphh 用
SEARCH = {"esdw": "w", "eslp": "w", "rk": "rk", "d5": "d5"}        # rost5 fd 的表示
METHODS = ["fifo", "pb", "esdw", "eslp", "gphh", "rk", "d5", "rddhand", "rhmip", "ours"]


def wstr(w):
    """权重 → rost5 的 w= 选项。"""
    return "w=" + ",".join(repr(float(x)) for x in w)


def lpg_window(name):
    """ES-LP 用的 LP 引导文件：LP_fr^Rw(4S) 的值岗流在窗口 D < W 内按站的份额（与 T3 的 lpg 同格式，只是换成 LP_fr^Rw）。没有就生成。"""
    os.makedirs(RULES13, exist_ok=True)
    p = os.path.join(RULES13, f"{name}.lpgw.txt")
    if not os.path.exists(p):
        from t6_rerun import frrw_npz
        import t3_make_rlp as M
        d = np.load(frrw_npz(name))
        tmp = p + f".{os.getpid()}.tmp"
        M.write_lpg(tmp, int(d["H"]), d["N"], d["M"], d["duty_flow"])     # 取 D < H/2 = 2·S_max = W 的流
        os.replace(tmp, p)
    return p


def opts(method, name, S, tree_file=None):
    """方法的 rost5 规则选项（不含 W=）。"""
    W = 2 * S
    if method == "fifo": return [wstr(W_FIFO), "facloop"]
    if method in ("pb", "rddhand"): return [wstr(W_PB), "facloop"]
    if method == "esdw": return ["curdef", "facloop"]
    if method == "eslp": return [f"lpg={lpg_window(name)}", f"planf={RULES6}/{name}__fresh_rlp_w.planf", f"planW={W}", "curdef", "facloop"]
    if method == "gphh": return [wstr(W_TOOL), f"tree=@{tree_file}", "facloop"]
    if method in ("rk", "d5"): return [wstr(W_TOOL), "facloop"]
    if method == "ours": return [f"rlp={RULES6}/{name}__fresh_rlp_w.rlp", f"planf={RULES6}/{name}__fresh_rlp_w.planf", f"planW={W}", f"hsfb={8 * S}"]
    raise ValueError(method)
