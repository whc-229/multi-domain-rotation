"""零余量 K=1：A 平移铺满 Z_S ⇔（L 为素数幂时）Coven–Meyerowitz 判据。
CM：S_A = {素数幂 s : Φ_s | A(x)}；T1：Π_{s∈S_A} Φ_s(1) = |A|；T2：S_A 中不同素数的幂之积 m 也满足 Φ_m | A。
铺满 Z_S 还要求 S_A 中的素数幂都整除 S（见笔记的引理）。"""
import json, cmath, math, itertools
def pp(n):  # 素数幂 -> (p,k) or None
    for p in range(2, n + 1):
        if n % p == 0:
            k = 0; m = n
            while m % p == 0: m //= p; k += 1
            return (p, k) if m == 1 else None
    return None
def phi_div(A, d):
    z = cmath.exp(2j * math.pi / d)
    return abs(sum(z ** a for a in A)) < 1e-7
def cm_tiles_ZS(A, S):
    L = len(A); SA = [s for s in range(2, S + 1) if pp(s) and phi_div(A, s)]
    if any(S % s for s in SA): return False, SA
    T1 = math.prod(pp(s)[0] for s in SA) == L
    T2 = True
    for r in range(2, 6):
        for comb in itertools.combinations(SA, r):
            if len(set(pp(s)[0] for s in comb)) == r and not phi_div(A, math.prod(comb)): T2 = False
    return T1 and T2, SA
def is_pp(n): return n == 1 or pp(n) is not None
rows = []
for src in ["atlas_random_80.json"]:
    for item in json.load(open(src)):
        g = item["g"]; S = sum(g); L = len(g)
        if S % L: continue
        A = [sum(g[:i]) for i in range(L)]
        r1 = item["rows"][0]
        tiles = r1.get("tiling"); gap0 = (r1["N_wait"] == r1["LB"])
        cm, SA = cm_tiles_ZS(A, S)
        rows.append((g, L, is_pp(L), gap0, tiles, cm, SA))
agree = sum(1 for r in rows if r[3] == r[5]); print("零余量 K=1 实例", len(rows), "CM 判据与网络 ILP 一致", agree)
for r in rows:
    if r[3] != r[5]: print("不一致:", r)
print("L 为素数幂的实例:", sum(1 for r in rows if r[2]), " 非素数幂:", sum(1 for r in rows if not r[2]))
for r in rows[:12]: print(r)
