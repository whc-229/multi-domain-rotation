import math, itertools, sys
from roster_sim import INSTANCES, Instance

def robot_occ_seq(inst, r):
    """单台 r 型机器人第 k 班的'占用时长'序列(从准备开始到再次可用, 含自身值岗)，按维修超周期给出"""
    p = inst.robots[r]; T = inst.tools
    unin_min = p["tool_num"] * min(T[t]["uninstall"] for t in p["tools"])
    base = p["prep"] + p["max_install"] + p["go"] + p["duty"] + p["post_consume"] + p["back"] + unin_min
    left = [m[0] for m in p["maint"]]
    seq, states = [], {}
    k = 0
    while True:
        key = tuple(round(l, 6) for l in left)
        if key in states:
            mu = states[key]; break
        states[key] = k
        left = [l - p["consume"] for l in left]
        rec = p["charge"]; lvl = None
        for lv in range(len(p["maint"]) - 1, -1, -1):
            if left[lv] < p["consume"] - 1e-9:
                rec = p["maint"][lv][1]
                for l2 in range(lv + 1): left[l2] = p["maint"][l2][0]
                lvl = lv; break
        seq.append(base + rec); k += 1
    return seq, mu   # seq[mu:] 为周期部分(超周期), seq[:mu] 为暂态

def occ_at(seq, mu, i):
    if i < len(seq): return seq[i]
    L = len(seq) - mu
    return seq[mu + (i - mu) % L]

def tool_occ(inst, r):
    p = inst.robots[r]; T = inst.tools
    unin_min = p["tool_num"] * min(T[t]["uninstall"] for t in p["tools"])
    base = p["max_install"] + p["go"] + p["duty"] + p["post_consume"] + p["back"] + unin_min
    return {t: base + T[t]["maintain"] for t in p["tools"]}

def max_shifts_slots(seq, mu, d, S):
    """单台机器人在 S 个等长班次槽中最多能值多少班（状态只依赖已值班次数 → 贪心最早即最优）"""
    pos, k = 0, 0
    while pos < S:
        pos += math.ceil(occ_at(seq, mu, k) / d - 1e-9); k += 1
    return k

def max_shifts_cont(seq, mu, H):
    pos, k = 0.0, 0
    while pos < H - 1e-9:
        pos += occ_at(seq, mu, k); k += 1
    return k

def report(name):
    inst = INSTANCES[name]()
    H = inst.horizon
    print(f"\n=== {name}  H={H}")
    rows = []
    for r, p in enumerate(inst.robots):
        seq, mu = robot_occ_seq(inst, r)
        cyc = seq[mu:]; d = p["duty"]
        rho = d * len(cyc) / sum(cyc)                       # 连续稳态值岗占比
        abar = sum(math.ceil(o / d - 1e-9) for o in cyc) / len(cyc)  # 槽化后每班平均占用槽数
        to = tool_occ(inst, r); tbest = min(to.values())
        tslots = min(math.ceil(v / d - 1e-9) for v in to.values())
        S = math.ceil(H / d)
        kH = max_shifts_slots(seq, mu, d, S)
        N_H = math.ceil(S / kH)
        M_inf = p["tool_num"] * tslots
        # 有限H 工具下界: 每把工具在 S 槽内最多用 ceil(S/tslots) 次
        M_H = math.ceil(S * p["tool_num"] / math.ceil(S / tslots))
        print(f" 类型{r}: d={d}, 超周期班数={len(cyc)}(暂态{mu}), 占用序列={[round(o,1) for o in cyc]}")
        print(f"   连续 1/rho={1/rho:.2f}  槽化 abar={abar:.2f} -> N_inf>={math.ceil(abar-1e-9)};"
              f" 工具占用≈{tbest:.1f}h={tslots}槽 -> M_inf>={M_inf};  有限H: N_H>={N_H}, M_H>={M_H}, 合计>={N_H+M_H} (∞: {math.ceil(abar-1e-9)+M_inf})")
        rows.append((r, N_H, M_H))
    return inst, rows

if __name__ == "__main__":
    for n in ["L1-720", "L2-720", "L3-720", "L1-1440", "L3-1440"]:
        report(n)
