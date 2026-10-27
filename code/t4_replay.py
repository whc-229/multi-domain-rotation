"""T4：上界解的独立重放。
1) 把 fresh_tail_mdrr 的周期整数解按周期 P 延拓到 Tmax（下标 τ ≥ T+P 的变量取 τ' = T + (τ − T) mod P 处的值）；
2) 分解成显式排班：值岗按出发时刻 v 升序，给每班分配一台"同型、同相位、已就绪"的机器人（同一节点上的机器人可互换，所以贪心总能成功，
   失败说明解不满足流守恒）；工具按占用起点升序，给每个组合分配空闲的具体工具；
3) 交给 rost5 replay：用仿真器自己的状态机逐条检查（机器人就绪、相位、提前大修合法、工具空闲、每个岗位到点都有人接班），
   报告缺岗数和违规数。0 缺岗、0 违规即说明该无限延拓排班在 [0, Tmax) 内可行（60 万格复核）。
用法：python3 code/t4_replay.py <ub 结果 json> [Tmax=600000]"""
import os, sys, json, math, heapq, tempfile, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, U, ToolInst, load_train, inst_paths  # noqa: E402
from fresh_tail_mdrr import combos  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")


def expand(items, T, P, idx_pos, tmax):
    """items: [[…, 下标, …, 数量]]；返回按下标展开到 tmax 的 (下标, 记录) 生成器。下标 < T+P 的直接用，之后按周期 P 重复 [T, T+P)。"""
    base = [it for it in items if it[idx_pos] < T + P]
    pat = [it for it in items if T <= it[idx_pos] < T + P]
    for it in base: yield it[idx_pos], it
    c = 1
    while True:
        any_ = False
        for it in pat:
            tau = it[idx_pos] + c * P
            if tau < tmax: any_ = True; yield tau, it
        if not any_: break
        c += 1


def build_schedule(I, sol, tmax):
    TI = ToolInst(I)
    T, P = sol["T"], sol["P"]
    nk = len(I.T)
    v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]
    # 机器人
    heaps = {(k, i): [] for k in range(nk) for i in range(I.L[k])}
    rid = 0
    for k in range(nk):
        for _ in range(sol["N"][k]): heapq.heappush(heaps[(k, 0)], (v0[k], rid)); rid += 1
    duties = sorted((tau, it[0], it[1], it[3], it[4], it[5]) for tau, it in expand(sol["duties"], T, P, 2, tmax))  # (v, k, i, s, r, 数量)
    assign = {}     # (k, s, D) → [(机器人, r)]
    for (v, k, i, s, r, c) in duties:
        for _ in range(c):
            hp = heaps[(k, i)]
            if not hp or hp[0][0] > v: raise RuntimeError(f"机器人分解失败 v={v} k={k} i={i}")
            ready, rob = heapq.heappop(hp)
            L = I.L[k]
            nx = I.step(k, L - 1 if r else i, s)[1]; ni = 0 if r else (i + 1) % L
            heapq.heappush(heaps[(k, ni)], (v + nx, rob))
            assign.setdefault((k, s, v + I.lead(s)), []).append([rob, r, None])
    # 工具
    th = {j: [] for j in range(TI.nt)}; uid = {j: 0 for j in range(TI.nt)}
    for j in range(TI.nt):
        for _ in range(sol["M"][j]): heapq.heappush(th[j], (-1000000, uid[j])); uid[j] += 1
    cb = {k: combos(TI, k) for k in range(nk)}
    zs = sorted((max(0, tau - TI.pre(it[0], it[1])), tau, it[0], it[1], it[3], it[4]) for tau, it in expand(sol["z"], T, P, 2, tmax))
    for (a, D, k, s, ci, c) in zs:
        cnt, unin = cb[k][ci]
        for _ in range(c):
            picks = []
            for j, n in cnt:
                for _ in range(n):
                    if not th[j] or th[j][0][0] > a: raise RuntimeError(f"工具分解失败 a={a} j={j}")
                    f, u = heapq.heappop(th[j]); picks.append((j, u))
            for (j, u) in picks:
                heapq.heappush(th[j], (D + I.d[k] + TI.post(k, s, j, unin), u))
            lst = assign.get((k, s, D), [])
            slot = next((x for x in lst if x[2] is None), None)
            if slot is None: raise RuntimeError(f"组合多于值岗 {(k, s, D)}")
            slot[2] = picks
    lines = []
    for (k, s, D), lst in assign.items():
        if D >= tmax: continue
        for rob, r, picks in lst:
            if picks is None: raise RuntimeError(f"值岗没有工具组合 {(k, s, D)}")
            lines.append((D, s, rob, r, picks))
    lines.sort(key=lambda x: (x[0], x[1], x[2]))
    return lines


def replay(name, sol, tmax=600000):
    I = load_train(name); _, txt = inst_paths(name, "train")
    lines = build_schedule(I, sol, tmax + 2000)
    with tempfile.NamedTemporaryFile("w", suffix=".sched", delete=False) as f:
        f.write("N " + " ".join(map(str, sol["N"])) + "\n" + "M " + " ".join(map(str, sol["M"])) + "\n")
        for (D, s, rob, r, picks) in lines:
            f.write(f"{D} {s} {rob} {r} {len(picks)} " + " ".join(f"{j} {u}" for j, u in picks) + "\n")
        path = f.name
    try:
        out = subprocess.run([R5, "replay", txt, path, str(tmax)], capture_output=True, text=True)
        return json.loads(out.stdout)
    finally:
        os.remove(path)


if __name__ == "__main__":
    sol = json.loads(open(sys.argv[1]).read())
    tmax = int(sys.argv[2]) if len(sys.argv) > 2 else 600000
    print(json.dumps(replay(sol["inst"], sol, tmax)))
