"""K=1 失败例：如果允许初始相位任意（机器人不是同时全新），N=⌈ā⌉ 能否无限期覆盖？区分'稳态不可行'和'启动同步造成的不可行'"""
import math, random, itertools
from slotmodel import Robot, canon
import kpost
random.seed(1)
cases = []
for trial in range(250):
    L = random.randint(1, 4)
    g = [random.randint(2, 6) for _ in range(L - 1)] + [random.randint(4, 16)]
    cases.append(g)
def exact_from(R, s0, K, lim):
    color = {canon(R, s0): 1}; nodes = 0
    stack = [(s0, iter(kpost.succK(R, s0, K)))]
    while stack:
        st, it = stack[-1]
        for t in it:
            nodes += 1
            if nodes > lim: return None
            c = canon(R, t); col = color.get(c, 0)
            if col == 1: return True
            if col == 0: color[c] = 1; stack.append((t, iter(kpost.succK(R, t, K)))); break
        else: color[canon(R, st)] = 2; stack.pop()
    return False
stats = {"全新可行": 0, "全新不可行_任意相位可行": 0, "任意相位也不可行": 0, "未决": 0}
ex_desync, ex_hard = [], []
for g in cases:
    L = len(g); a = sum(g) / L; N = math.ceil(a - 1e-9)
    R = [Robot(g, 0, 0) for _ in range(N)]
    ok, _ = kpost.exact(R, 1, 150000)
    if ok is None: stats["未决"] += 1; continue
    if ok: stats["全新可行"] += 1; continue
    found = False; und = False
    for phs in itertools.combinations_with_replacement(range(L), N):
        s0 = tuple((0, p) for p in phs)
        r = exact_from(R, s0, 1, 3000)
        if r is None: und = True
        if r: found = True; break
    if found: stats["全新不可行_任意相位可行"] += 1; ex_desync.append((g, round(a, 2), N))
    elif und: stats["未决"] += 1
    else: stats["任意相位也不可行"] += 1; ex_hard.append((g, round(a, 2), N, round(N - a, 2)))
print(stats)
print("启动同步造成的例子:", ex_desync[:15])
print("稳态也不可行的例子（g, ā, N=⌈ā⌉, 余量）:")
for e in sorted(ex_hard, key=lambda e: e[3])[:60]: print(e)
