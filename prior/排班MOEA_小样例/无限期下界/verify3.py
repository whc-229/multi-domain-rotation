import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import math, re, inspect, numpy as np, sys
import roster_sim as rs
from roster_sim import INSTANCES, Instance, evaluate
from lb import robot_occ_seq
sys.setrecursionlimit(100000)

# ---------- 1) 精确有限状态图: 找可达环, 并记录每槽选择的相位
def plan(inst, r, N):
    seq, mu = robot_occ_seq(inst, r); d = inst.robots[r]["duty"]
    g = [math.ceil(o/d - 1e-9) for o in seq]; L = len(g)
    nph = lambda ph: ph+1 if ph+1 < L else mu
    def succ(s):
        out = {}
        for (w, ph) in set(s):
            if w == 0:
                lst = list(s); lst.remove((w, ph))
                new = [(max(w2-1, 0), p2) for (w2, p2) in lst] + [(g[ph]-1, nph(ph))]
                out.setdefault(tuple(sorted(new)), ph)
        return list(out.items())
    start = tuple(sorted([(0, 0)]*N)); color = {start: 1}
    stack = [(start, iter(succ(start)))]; path = [start]; ch = []
    while stack:
        s, it = stack[-1]
        for t, ph in it:
            c = color.get(t, 0)
            if c == 1:
                k = path.index(t); ch.append(ph); return ch, k, len(path) - k, mu, L
            if c == 0:
                color[t] = 1; stack.append((t, iter(succ(t)))); path.append(t); ch.append(ph); break
        else:
            color[s] = 2; stack.pop(); path.pop()
            if ch: ch.pop()
    return None

# ---------- 2) evaluate 的可挂钩副本
src = inspect.getsource(rs.evaluate)
src = src.replace("def evaluate(", "def evaluate_hook(POL_R, POL_T, LOG, ")
src = re.sub(r"sel = int\(x\[k\] \* \(len\(avail\) \+ nT\)\)\n\s*sel = min\(max\(sel, 0\), len\(avail\) \+ nT - 1\)",
             "sel = POL_R(k, avail, robots); LOG.append([sel, len(avail) + nT])", src)
src = src.replace("off = int(x[k + S] * len(cand)) % len(cand)",
                  "off = POL_T(k, cand, tools, R[rb['type']]['tool_num']); LOG[-1] += [off, len(cand)]")
ns = dict(rs.__dict__); exec(src, ns); evaluate_hook = ns["evaluate_hook"]

def run(name, H, r, N):
    base = INSTANCES[name](); inst = Instance(f"{name}@{H}", H, base.level, tool_dec=TOOL_DEC)
    pl = plan(inst, r, N); ch, k0, per, mu, L = pl
    def phase(n): return n if n < L else mu + (n - mu) % (L - mu)
    def pol_r(k, avail, robots):
        want = ch[k] if k < len(ch) else ch[k0 + (k - k0) % per]
        for i, j in enumerate(avail):
            if robots[j]["type"] == r and phase(robots[j]["duties"]) == want: return i
        if want == 0: return len(avail) + r        # 新开一台
        raise RuntimeError(f"plan mismatch k={k} want={want} avail_ph={[(phase(robots[j]['duties'])) for j in avail]} n={len(robots)} allph={[(phase(rb['duties']),round(rb['avail'],1)) for rb in robots]}")
    TT = inst.tools
    def pol_t(k, cand, tools, need):
        best = None
        for off in range(len(cand)):
            new = 0; got = 0; un = 0
            for c in range(len(cand)):
                kind, v = cand[(off + c) % len(cand)]
                if got == need: break
                if kind == "old":
                    got += 1; un += TT[tools[v]["type"]]["uninstall"]
                else:
                    n_new = 1 if inst.tool_dec == "orig" else need - got
                    got += n_new; new += n_new; un += n_new * TT[v]["uninstall"]
            key = (un, new)
            if best is None or key < best[0]: best = (key, off)
        return best[1]
    LOG = []
    x = np.zeros(inst.n_var)
    evaluate_hook(pol_r, pol_t, LOG, x, inst)
    for k, (sel, nr, off, nc) in enumerate(LOG):
        x[k] = (sel + 0.5) / nr; x[k + inst.S] = (off + 0.5) / nc
    out = evaluate(x, inst)       # 用原始 evaluate 独立复核
    return out, (k0, per)

import sys
TOOL_DEC = sys.argv[1] if len(sys.argv) > 1 else 'group'
for name, N in [("L1-720", 7), ("L1-1440", 7), ("L2-720", 5), ("L3-720", 6), ("L3-1440", 6)]:
    for H in [INSTANCES[name]().horizon, 20000]:
        try:
            out, (k0, per) = run(name, H, 0, N)
        except RuntimeError as e:
            print(f"{name:8s} H={H:6.0f}: 解码器无法复现该周期计划"); continue
        print(f"{name:8s} H={H:6.0f}: 前缀{k0}班+周期{per}班  原evaluate可行={out['feasible']}  机器人={out.get('n_robot')} 工具={out.get('n_tool')} f={np.round(out['f'],2)}")
