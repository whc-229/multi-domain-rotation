import os, sys; _h=os.path.dirname(os.path.abspath(__file__)); sys.path[:0]=[_h, os.path.join(_h,"..","启发式验证"), os.path.join(_h,"..")]
"""K=1 失败例：如果允许初始相位任意（机器人不是同时全新），N=⌈ā⌉ 能否无限期覆盖？区分'稳态不可行'和'启动同步造成的不可行'"""
import math, random, itertools
from slotmodel import Robot, canon
import kpost
random.seed(1)
cases = []
for trial in range(0):
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
if __name__=='__main__':
    pass
