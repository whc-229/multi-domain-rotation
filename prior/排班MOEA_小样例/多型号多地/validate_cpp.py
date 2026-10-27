"""C++ 与 Python 的逐项对齐：同一实例、同一配比、同一权重，比较缺岗数（固定时域）与证书/周期（判环模式）。"""
import subprocess, json, random
import numpy as np
from export_inst import build
from tools_model import ToolInst
from certsim import run, start_point
random.seed(1); rng = np.random.default_rng(1)
tot = 0; bad = 0
for spec in ["L1", "L3", "2s1", "3s1"]:
    I = build(spec); TI = ToolInst(I); N0, M0 = start_point(I, TI)
    for trial in range(12):
        Nk = [n + random.randint(0, 2) for n in N0]; Mj = [m + random.randint(0, 2) for m in M0]
        w = rng.normal(0, 1, 8).tolist() if trial > 1 else [0.1, 0, 10, 1, 0, 0, 1, 0.1]
        for Tmax, cert in [(3000, False), (60000, True)]:
            py = run(I, Nk, Mj, w, Tmax, TI, certify=cert)
            args = ["cpp/rost", "eval", f"inst/{spec}.txt", str(Tmax), str(int(cert))] + [str(x) for x in Nk] + ["|"] + [str(x) for x in Mj] + ["|"] + [repr(x) for x in w]
            cc = json.loads(subprocess.run(args, capture_output=True, text=True).stdout)
            pc = {None: -1, False: 0, True: 1}[py["certified"]]
            if not cert: ok = (py["gaps"] == cc["gaps"])
            else: ok = (pc == cc["cert"]) and (pc != 1 or py["period"] == cc["period"]) and (pc == -1 and py["gaps"] == cc["gaps"] or pc != -1)
            tot += 1; bad += (not ok)
            if not ok: print("不一致", spec, trial, Tmax, cert, py["gaps"], py["certified"], py["period"], cc)
print("比对", tot, "组，不一致", bad)
