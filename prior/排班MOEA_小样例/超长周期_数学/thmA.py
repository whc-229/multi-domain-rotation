import os, sys; _h=os.path.dirname(os.path.abspath(__file__)); sys.path[:0]=[_h, os.path.join(_h,"..","启发式验证"), os.path.join(_h,"..")]
"""定理 A 数值核对：K 为 L 的倍数时，稳态（任意相位、等待为 0 的初态）下 N = K·ā = K·S/L 台可行。
构造：S 台机器人，第 t 台从时刻 t 开始按间隔序列 g 循环工作 ⇒ 每格恰好 L 台。直接验证构造 + 用精确搜索再确认。"""
import itertools, math
from slotmodel import Robot, canon
import kpost

def construct_check(g, m):
    S, L = sum(g), len(g); A = [0]
    for x in g[:-1]: A.append(A[-1] + x)
    cover = [0] * S
    for rep in range(m):
        for t in range(S):
            for a in A: cover[(t + a) % S] += 1
    return all(c == m * L for c in cover)

def steady_exact(g, N, K, lim=300000):
    R = [Robot(g, 0, 0) for _ in range(N)]; L = len(g)
    # 构造给出的相位分布：直接用构造的初态（第 t 台在时刻 0 的状态）
    S = sum(g); A = [0]
    for x in g[:-1]: A.append(A[-1] + x)
    st = []
    for rep in range(N // S):
        for t in range(S):
            # 时刻 0 时第 t 台：找它下一次工作在 A+t+SZ 中 ≥0 的最小值，相位为对应下标
            nxt = min(((t + a) % S, i) for i, a in enumerate(A))
            st.append((nxt[0], nxt[1]))
    from conj2 import exact_from
    return exact_from(R, tuple(st), K, lim)

for g in [[3,10],[2,9],[4,4,5],[2,2,5],[3,5,2,6],[4,4,5,4,9],[3,3,15]]:
    S, L = sum(g), len(g)
    for m in [1, 2]:
        K = m * L; N = m * S
        print(g, f"K={K}, N=K·ā={N}", "构造覆盖每格恰好K台:", construct_check(g, m), "| 精确搜索从构造初态找到环:", steady_exact(g, N, K) if N <= 26 else "跳过")
