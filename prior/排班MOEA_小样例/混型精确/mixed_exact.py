import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
"""混型（值岗时长不同）+ 工具的精确可行性搜索：事件时间状态（相对当前值岗开始时刻），时间单位 0.1h。
与 roster_sim.evaluate 的规则一致：
  机器人可用 ⇔ avail ≤ start − go − max_install − prep；工具可用 ⇔ avail ≤ start − go − max_install（按所装机器人型号）
  卸载时长 = 所带工具卸载时长之和；工具在 p4_end + 保养 后可用；维修按各级剩余工时判断。
给定配置（各型机器人数、各型工具数），判断：
  无限期：从全新状态出发能否进入环；有限期：能否覆盖 H 小时。"""
import itertools, math, sys
from roster_sim import INSTANCES

def T(x): return int(round(x * 10))

class Model:
    def __init__(self, inst):
        self.inst = inst
        self.R = []
        for r in inst.robots:
            self.R.append(dict(d=T(r["duty"]), lead=T(r["go"] + r["max_install"] + r["prep"]),
                               tlead=T(r["go"] + r["max_install"]), post_back=T(r["post_consume"] + r["back"]),
                               charge=T(r["charge"]), consume=T(r["consume"]), m=r["tool_num"], tools=r["tools"],
                               maint=[(T(a), T(b)) for a, b in r["maint"]]))
        self.Tl = [dict(un=T(t["uninstall"]), mt=T(t["maintain"])) for t in inst.tools]
        self.tclamp = max(r["tlead"] for r in self.R)

    def start(self, nr, nt):
        robots = tuple(sorted((r, 0, tuple(a for a, _ in self.R[r]["maint"])) for r in range(len(nr)) for _ in range(nr[r])))
        tools = tuple(sorted((t, -self.tclamp) for t in range(len(nt)) for _ in range(nt[t])))
        return (robots, tools)

    def succ(self, st):
        robots, tools = st
        out = []
        seen_r = set()
        for i, (r, w, left) in enumerate(robots):
            if w != 0 or (r, left) in seen_r: continue
            seen_r.add((r, left)); P = self.R[r]
            # 可用工具类别（类型, 相对可用时刻）
            avail = [k for k, (t, a) in enumerate(tools) if t in P["tools"] and a <= -P["tlead"]]
            cls = {}
            for k in avail: cls.setdefault(tools[k], []).append(k)
            keys = list(cls)
            combos = set()
            def rec(q, left_n, pick):
                if left_n == 0: combos.add(tuple(pick)); return
                if q == len(keys): return
                for c in range(min(left_n, len(cls[keys[q]])), -1, -1):
                    rec(q + 1, left_n - c, pick + cls[keys[q]][:c])
            rec(0, P["m"], [])
            for combo in combos:
                out.append(self.apply(st, i, combo))
        return out

    def apply(self, st, i, combo):
        robots, tools = st
        r, w, left = robots[i]; P = self.R[r]; d = P["d"]
        un = sum(self.Tl[tools[k][0]]["un"] for k in combo)
        p4 = d + P["post_back"] + un                      # 相对本班开始
        left2 = [l - P["consume"] for l in left]; rec = P["charge"]
        for lv in range(len(P["maint"]) - 1, -1, -1):
            if left2[lv] < P["consume"]:
                rec = P["maint"][lv][1]
                for l2 in range(lv + 1): left2[l2] = P["maint"][l2][0]
                break
        avail_r = p4 + rec
        newR = []
        for k, (r2, w2, l2) in enumerate(robots):
            if k == i:
                # 相对下一班开始(=+d)：需要 avail ≤ start' − lead
                newR.append((r2, max(0, avail_r - d + self.R[r2]["lead"]), tuple(left2)))
            else:
                newR.append((r2, max(0, w2 - d), l2))
        newT = []
        cs = set(combo)
        for k, (t, a) in enumerate(tools):
            a2 = (p4 + self.Tl[t]["mt"]) if k in cs else a
            newT.append((t, max(-self.tclamp, a2 - d)))
        return (tuple(sorted(newR)), tuple(sorted(newT))), d

def infinite_feasible(M, nr, nt, node_limit=200000):
    s0 = M.start(nr, nt); color = {s0: 1}; nodes = 0
    stack = [(s0, iter(M.succ(s0)))]
    while stack:
        st, it = stack[-1]
        for (t, d) in it:
            nodes += 1
            if nodes > node_limit: return None, nodes
            c = color.get(t, 0)
            if c == 1: return True, nodes
            if c == 0:
                color[t] = 1; stack.append((t, iter(M.succ(t)))); break
        else:
            color[st] = 2; stack.pop()
    return False, nodes

def finite_feasible(M, nr, nt, H, node_limit=200000):
    Ht = T(H); s0 = M.start(nr, nt); dead = set(); nodes = 0
    stack = [(s0, 0, iter(M.succ(s0)))]
    while stack:
        st, el, it = stack[-1]
        for (t, d) in it:
            nodes += 1
            if nodes > node_limit: return None, nodes
            e2 = el + d
            if e2 >= Ht: return True, nodes
            if (t, e2) in dead: continue
            stack.append((t, e2, iter(M.succ(t)))); break
        else:
            dead.add((st, el)); stack.pop()
    return False, nodes
