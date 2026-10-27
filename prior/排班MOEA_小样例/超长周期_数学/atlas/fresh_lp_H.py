from fresh_lp3 import fresh_lp
from fractions import Fraction
for g in [[3,10],[2,9],[4,4,5],[3,5,2,6],[2,2,5],[3,3,15],[4,6,5,4],[2,6,3],[1,4],[2,4]]:
    S=sum(g); row=[]
    for m in [1,1.5,2,3,4,6,8]:
        H=int(m*S); v=fresh_lp(g,H); row.append(f"{m}S:{Fraction(v).limit_denominator(500)}")
    print(g, round(S/len(g),4), row, flush=True)
