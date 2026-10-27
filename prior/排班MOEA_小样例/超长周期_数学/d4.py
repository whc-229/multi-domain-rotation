import os, sys; _h=os.path.dirname(os.path.abspath(__file__)); sys.path[:0]=[os.path.join(_h,'..','启发式验证'), os.path.join(_h,'..')]
"""D4：允许提前维修。每班结束时，机器人可选择：充电（若未被强制维修），或任一级维修 ℓ（≥ 被强制的级别），
维修 ℓ 重置第 0..ℓ 级剩余工时。占用 = 基础占用 + 恢复时长，按值岗时长取整成槽。全新启动，同型，K 个点。"""
import math, itertools, sys, json
from roster_sim import INSTANCES

def params(name, r):
    inst = INSTANCES[name](); p = inst.robots[r]; T = inst.tools
    unin = p["tool_num"] * min(T[t]["uninstall"] for t in p["tools"])
    base = p["prep"] + p["max_install"] + p["go"] + p["duty"] + p["post_consume"] + p["back"] + unin
    return dict(d=p["duty"], base=base, consume=p["consume"], charge=p["charge"], maint=p["maint"])

def options(P, left, early):
    """返回 [(槽间隔, 新剩余工时)]"""
    c = P["consume"]; l2 = [round(x - c, 6) for x in left]; M = P["maint"]
    forced = None
    for lv in range(len(M) - 1, -1, -1):
        if l2[lv] < c - 1e-9: forced = lv; break
    outs = []
    def do(lv):
        nl = list(l2)
        for k in range(lv + 1): nl[k] = M[k][0]
        return (math.ceil((P["base"] + M[lv][1]) / P["d"] - 1e-9), tuple(nl))
    if forced is None:
        outs.append((math.ceil((P["base"] + P["charge"]) / P["d"] - 1e-9), tuple(l2)))
        if early: outs += [do(lv) for lv in range(len(M))]
    else:
        outs.append(do(forced))
        if early: outs += [do(lv) for lv in range(forced + 1, len(M))]
    return list(dict.fromkeys(outs))

def feasible(P, N, K, early, lim=400_000):
    full = tuple(m[0] for m in P["maint"])
    s0 = tuple(sorted([(0, full)] * N))
    def succ(st):
        av = [i for i, (w, _) in enumerate(st) if w == 0]
        if len(av) < K: return []
        cls = {}
        for i in av: cls.setdefault(st[i], []).append(i)
        keys = list(cls); out = set()
        # 选 K 台（按类别多重集），每台再选恢复方式
        def rec(q, left_n, pick):
            if left_n == 0:
                chosen = pick
                opt_lists = [options(P, st[i][1], early) for i in chosen]
                for combo in itertools.product(*opt_lists):
                    new = []
                    cs = set(chosen)
                    for i, (w, l) in enumerate(st):
                        if i not in cs: new.append((max(w - 1, 0), l))
                    for (gap, nl) in combo: new.append((gap - 1, nl))
                    out.add(tuple(sorted(new)))
                return
            if q == len(keys): return
            for c in range(min(left_n, len(cls[keys[q]])), -1, -1):
                rec(q + 1, left_n - c, pick + cls[keys[q]][:c])
        rec(0, K, [])
        return list(out)
    color = {s0: 1}; nodes = 0; stack = [(s0, iter(succ(s0)))]
    while stack:
        st, it = stack[-1]
        for t in it:
            nodes += 1
            if nodes > lim: return None
            c = color.get(t, 0)
            if c == 1: return True
            if c == 0: color[t] = 1; stack.append((t, iter(succ(t)))); break
        else: color[st] = 2; stack.pop()
    return False

def nstar(P, K, early, lo):
    for N in range(lo, lo + 8):
        r = feasible(P, N, K, early)
        if r is None: return ("?", N)
        if r: return N
if __name__ == "__main__":
    out = []
    for name in ["L1-720", "L2-720", "L3-720"]:
        for r in [0, 1, 2]:
            P = params(name, r)
            for K in [1, 2]:
                lo = max(K + 1, 2)
                a = nstar(P, K, False, lo); b = nstar(P, K, True, lo)
                line = dict(inst=name, type=r, K=K, 强制维修=a, 允许提前维修=b)
                out.append(line); print(json.dumps(line, ensure_ascii=False), flush=True)
    json.dump(out, open("d4_results.json", "w"), ensure_ascii=False, indent=1)
