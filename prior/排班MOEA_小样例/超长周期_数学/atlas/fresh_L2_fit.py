"""L=2 全新启动流体超额 e(a,G)=α−S/2 的公式拟合核对：猜想 e = a·m/(2S−d)，m=min(r,2a−r)，r=S mod 2a。
对每个 (a,r) 用多个 q 值反解 d，看 d 是否只依赖 (a,r)。同时检验 H=2S 与 H=4S 的值一致。"""
from fresh_lp3 import fresh_lp
from fractions import Fraction as F
import json
out = []
for a in range(1, 8):
    for G in range(a, 10 * a + 2):
        S = a + G; r = S % (2 * a); m = min(r, 2 * a - r)
        v2 = F(fresh_lp([a, G], 2 * S)).limit_denominator(4000)
        v4 = F(fresh_lp([a, G], 4 * S)).limit_denominator(4000)
        e = v4 - F(S, 2)
        d = (2 * S - F(a * m) / e) if e else None
        out.append(dict(a=a, G=G, S=S, r=r, e=str(e), d=str(d), same2S4S=(v2 == v4)))
        print(a, G, S, r, e, d, v2 == v4, flush=True)
json.dump(out, open("fresh_L2_fit.json", "w"), indent=0)
