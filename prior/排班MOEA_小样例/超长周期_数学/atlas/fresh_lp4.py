import json
from fractions import Fraction
from fresh_lp3 import fresh_lp
rows = []
for L in [3, 4]:
    for c in [1, 2, 3]:
        for G in range(c, 3 * c * L + 3):
            g = [c] * (L - 1) + [G]; S = sum(g); a = S / L
            v = fresh_lp(g, 6 * S)
            rows.append((L, c, G, S, round(v - a, 4), S % (L * c)))
for r in rows: print(r)
