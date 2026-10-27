"""机器人—工具排班仿真（Python 版），仿照：
  - 算法库/目标函数/排班/problem.cpp  的 robot_arrangement / robot_arrangement_moea
  - 目标函数智能配置/graduate_code_judge_double/Evaluate.py  的逐班次解码逻辑

场景：值班地点必须始终且仅有 1 台机器人值岗，连续覆盖 horizon 小时。
每个班次派 1 台机器人，携带固定数量的适配工具。
  机器人流程：准备节点(P100…) → P0 安装工具(各工具安装时长串行相加) → P1 去程 → P2 值岗
              → 收尾节点(P203…) → P3 返程 → P4 卸载 → A4 扣减各级维修剩余可用时间
              → J3/J2/J1 多级维修判断（最严重的先判）→ 维修或 P5 充电 → 可用
  工具流程：P0 安装 → 随机器人去程/值岗/返程 → P4 卸载 → 保养 → 可用

编码（随机键，长度 2×最大班次数）：
  前半段 x[i]：第 i 班次选机器人，int(x·(可用机器人数+机器人类型数))，落在可用列表内则复用，否则新开该类型
  后半段 x[i+S]：第 i 班次工具选择偏移，候选 = 可用的已有适配工具 + 各适配类型"新开"选项，从偏移处循环取满所需数量
目标（与 robot_arrangement_moea 一致）：
  f1 = 机器人数 + 工具数（min_num）
  f2 = 0.3·std(各机器人值岗次数) + 0.4·std(各机器人值岗总时长) + 0.3·std(相邻"首次维修开始时间"间隔)（max_balance）
"""
import math
from dataclasses import dataclass, field

import numpy as np

# ------------------------------------------------------------------ 实例定义
TOOLS = [  # (数量上限, 安装, 卸载, 保养)  —— view.md 工具参数总表
    dict(max_num=10, install=1, uninstall=1, maintain=8),
    dict(max_num=10, install=2, uninstall=2, maintain=6),
    dict(max_num=10, install=1, uninstall=1, maintain=4),
    dict(max_num=10, install=2, uninstall=2, maintain=3),
]

# view.md 机器人参数总表 + Evaluate.py 的三级维修参数
ROBOTS = [
    dict(max_num=10, tools=[1, 2, 3], tool_num=2, go=1, back=1, charge=3, duty=10,
         prep_consume=0.3 + 0.4 + 0.2, prep_free=0.5 + 0.1, post_consume=0.2 + 0.1,
         maint1=[(36, 72)],                                   # 一级：(连续可用时长, 维修时长)
         maint3=[(25, 20), (40, 30), (70, 70)]),              # 一/二/三级
    dict(max_num=10, tools=[0, 1], tool_num=3, go=2, back=2, charge=4, duty=6,
         prep_consume=0.3 + 0.2 + 0.1, prep_free=0.6 + 0.1, post_consume=0.1 + 0.1,
         maint1=[(50, 100)],
         maint3=[(35, 10), (55, 40), (80, 100)]),
    dict(max_num=10, tools=[2, 3], tool_num=2, go=2, back=2, charge=3, duty=10,
         prep_consume=0.3 + 0.3 + 0.1, prep_free=0.8 + 0.1, post_consume=0.3 + 0.1,
         maint1=[(60, 120)],
         maint3=[(45, 30), (80, 60), (100, 120)]),
]


@dataclass
class Instance:
    name: str
    horizon: float
    level: int                 # 维修级数 1/2/3
    stock_scale: float = 1.0   # 资源上限缩放（<1 更紧）
    decoder: str = "orig"      # orig: 原解码；theta: 阈值解码；rank: 阈值 + 可用机器人按出勤次数升序
    theta: float = 0.85        # 阈值解码中"复用"区间的宽度
    tool_dec: str = "orig"     # orig: 原工具解码（遇到"新开"选项每类只开 1 把）；group: 遇到"新开 t 型"时一次开够所缺数量
    robots: list = field(default_factory=list)
    tools: list = field(default_factory=list)

    def __post_init__(self):
        self.robots, self.tools = [], []
        for r in ROBOTS:
            r = dict(r)
            r["maint"] = r["maint1"] if self.level == 1 else r["maint3"][: self.level]
            r["max_num"] = max(1, int(round(r["max_num"] * self.stock_scale)))
            # 每班次消耗的连续可用工时 = 准备中消耗的部分 + 去程 + 值岗 + 收尾 + 返程
            r["consume"] = r["prep_consume"] + r["go"] + r["duty"] + r["post_consume"] + r["back"]
            r["prep"] = r["prep_consume"] + r["prep_free"]
            r["max_install"] = max(TOOLS[t]["install"] for t in r["tools"]) * r["tool_num"]
            self.robots.append(r)
        for t in TOOLS:
            t = dict(t)
            t["max_num"] = max(1, int(round(t["max_num"] * self.stock_scale)))
            self.tools.append(t)
        self.min_duty = min(r["duty"] for r in self.robots)
        self.S = math.ceil(self.horizon / self.min_duty) + 1           # 最大班次数
        self.n_var = 2 * self.S
        # 第一班值岗开始时刻：保证任何类型都来得及准备（避免负时间）
        self.T0 = max(r["prep"] + r["max_install"] + r["go"] for r in self.robots)


INSTANCES = {
    "L1-720":   lambda: Instance("L1-720", 720, 1),
    "L2-720":   lambda: Instance("L2-720", 720, 2),
    "L3-720":   lambda: Instance("L3-720", 720, 3),
    "L3-720紧": lambda: Instance("L3-720紧", 720, 3, stock_scale=0.6),
    "L1-1440":  lambda: Instance("L1-1440", 1440, 1),
    "L3-1440":  lambda: Instance("L3-1440", 1440, 3),
}

# ------------------------------------------------------------------ 仿真
INFEAS = 30000.0


def evaluate(x, inst, detail=False):
    """返回 dict: f(2,), feasible, progress(已覆盖比例), fail_shift, shift_robot(每班次机器人编号), ..."""
    R, T, S = inst.robots, inst.tools, inst.S
    res_r = [r["max_num"] for r in R]
    res_t = [t["max_num"] for t in T]
    robots = []   # 每台: type, avail, left[级], duties, duty_time, first_repair
    tools = []    # 每件: type, avail
    duty_start = inst.T0
    end_time = inst.T0 + inst.horizon
    shift_robot, shift_tools = [], []
    n_avail_log, n_oldtool_log, n_cand_log, new_shift = [], [], [], []
    k = 0
    fail = None
    while duty_start < end_time - 1e-9:
        if k >= S:
            fail = k; break
        # ---- 2.1 可用机器人（时间来得及 + 各级剩余可用时长足够）
        avail = []
        for j, rb in enumerate(robots):
            p = R[rb["type"]]
            latest_prep_start = duty_start - p["go"] - p["max_install"] - p["prep"]
            if rb["avail"] <= latest_prep_start + 1e-9 and min(rb["left"]) >= p["consume"] - 1e-9:
                avail.append(j)
        # ---- 2.2 解码机器人
        nT = len(R)
        n_avail_log.append(len(avail))
        if inst.decoder == "orig":
            sel = int(x[k] * (len(avail) + nT))
            sel = min(max(sel, 0), len(avail) + nT - 1)
        else:
            if inst.decoder == "rank":   # 最闲的排前面（出勤次数少、其次剩余工时多）
                avail.sort(key=lambda j: (robots[j]["duties"], -min(robots[j]["left"])))
            th = inst.theta
            if avail and x[k] < th:
                sel = min(int(x[k] / th * len(avail)), len(avail) - 1)
            else:
                u = (x[k] - th) / (1 - th) if avail else x[k]
                sel = len(avail) + min(int(max(u, 0) * nT), nT - 1)
        if sel < len(avail):
            j = avail[sel]
        else:
            ty = sel - len(avail)
            for _ in range(nT):
                if res_r[ty] > 0:
                    break
                ty = (ty + 1) % nT
            else:
                fail = k; break
            res_r[ty] -= 1
            robots.append(dict(type=ty, avail=-1e9, left=[m[0] for m in R[ty]["maint"]],
                               duties=0, duty_time=0.0, first_repair=None))
            j = len(robots) - 1
            new_shift.append(k)
        rb = robots[j]; p = R[rb["type"]]
        # ---- 2.3 解码工具：候选 = 时间来得及的已有适配工具 + 各适配类型的"新开"选项
        install_start_bound = duty_start - p["go"] - p["max_install"]
        cand = [("old", t) for t, tl in enumerate(tools)
                if tl["type"] in p["tools"] and tl["avail"] <= install_start_bound + 1e-9]
        n_oldtool_log.append(len(cand))
        cand += [("new", ty) for ty in p["tools"]]
        n_cand_log.append(len(cand))
        off = int(x[k + S] * len(cand)) % len(cand)
        chosen = []
        for c in range(len(cand)):
            kind, v = cand[(off + c) % len(cand)]
            if len(chosen) == p["tool_num"]:
                break
            if kind == "old":
                chosen.append(v)
            else:
                # orig: 每个"新开"选项只开 1 把 → 一次需要开 2 把时不可能同型（会被迫混入卸载更慢的类型）
                # group: 一个"新开 v 型"选项开够当前所缺数量（受库存限制）
                n_new = 1 if inst.tool_dec == "orig" else p["tool_num"] - len(chosen)
                for _ in range(n_new):
                    if res_t[v] <= 0:
                        break
                    res_t[v] -= 1; tools.append(dict(type=v, avail=-1e9)); chosen.append(len(tools) - 1)
        while len(chosen) < p["tool_num"]:          # 一轮不够：继续从有库存的适配类型新开
            ok = False
            for ty in p["tools"]:
                if res_t[ty] > 0 and len(chosen) < p["tool_num"]:
                    res_t[ty] -= 1; tools.append(dict(type=ty, avail=-1e9)); chosen.append(len(tools) - 1); ok = True
            if not ok:
                break
        if len(chosen) < p["tool_num"]:
            fail = k; break
        # ---- 2.4 执行班次（时间从值岗开始时刻反推）
        inst_sum = sum(T[tools[t]["type"]]["install"] for t in chosen)
        unin_sum = sum(T[tools[t]["type"]]["uninstall"] for t in chosen)
        duty_end = duty_start + p["duty"]
        p4_end = duty_end + p["post_consume"] + p["back"] + unin_sum
        for t in chosen:
            tools[t]["avail"] = p4_end + T[tools[t]["type"]]["maintain"]
        rb["duties"] += 1; rb["duty_time"] += p["duty"]
        # A4 扣减 → 多级维修判断（最严重级先判）
        rb["left"] = [l - p["consume"] for l in rb["left"]]
        t_next = p4_end + p["charge"]
        for lv in range(len(p["maint"]) - 1, -1, -1):
            if rb["left"][lv] < p["consume"] - 1e-9:
                t_next = p4_end + p["maint"][lv][1]
                for l2 in range(lv + 1):
                    rb["left"][l2] = p["maint"][l2][0]
                if rb["first_repair"] is None:
                    rb["first_repair"] = p4_end
                break
        rb["avail"] = t_next
        shift_robot.append(j); shift_tools.append(chosen)
        duty_start = duty_end
        k += 1

    progress = min(1.0, (duty_start - inst.T0) / inst.horizon)
    out = dict(feasible=fail is None, progress=progress, fail_shift=fail,
               shift_robot=shift_robot, n_shift=k, n_avail=n_avail_log,
               n_oldtool=n_oldtool_log, n_cand=n_cand_log, new_shift=new_shift)
    if fail is not None:
        out["f"] = np.array([INFEAS, INFEAS])
        return out
    f1 = len(robots) + len(tools)
    cnt = np.array([r["duties"] for r in robots], float)
    dur = np.array([r["duty_time"] for r in robots], float)
    rep = np.sort([r["first_repair"] for r in robots if r["first_repair"] is not None])
    sd = lambda v: float(np.std(v)) if len(v) > 1 else 0.0
    f2 = 0.3 * sd(cnt) + 0.4 * sd(dur) + 0.3 * sd(np.diff(rep) if len(rep) > 1 else [])
    out["f"] = np.array([f1, f2])
    out["n_robot"], out["n_tool"] = len(robots), len(tools)
    out["robot_duties"] = cnt
    return out


if __name__ == "__main__":
    import time
    rng = np.random.default_rng(0)
    for nm, mk in INSTANCES.items():
        inst = mk()
        t = time.time(); fs = []; feas = 0
        for _ in range(300):
            o = evaluate(rng.random(inst.n_var), inst)
            feas += o["feasible"]
            if o["feasible"]:
                fs.append(o["f"])
        fs = np.array(fs) if fs else np.zeros((0, 2))
        print(f"{nm:9s} 变量数={inst.n_var:4d} 随机解可行率={feas/300:.2f} "
              f"f1范围={fs[:,0].min() if len(fs) else '-'}~{fs[:,0].max() if len(fs) else '-'} "
              f"单次评价={1000*(time.time()-t)/300:.2f}ms")
