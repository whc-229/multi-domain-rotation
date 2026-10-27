"""复核 C2 的疑似反例 [5,2,2,11,3]、K=1（W=2）：周期 P 放大到 6S，看是否存在 5 台的周期排班。"""
import time, json
from net import solve
g = [5, 2, 2, 11, 3]; K = 1; S = sum(g)
for P in range(5, 6 * S + 1, 5):
    t = time.time(); v, x = solve(g, K, P, time_limit=120)
    print(json.dumps(dict(P=P, N=None if v is None else round(v, 3), sec=round(time.time() - t, 1))), flush=True)
    if v is not None and round(v) <= 5: print("找到 5 台的周期排班，C2 未被否定"); break
