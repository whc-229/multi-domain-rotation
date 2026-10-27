"""T0：把 prior/ 里 validate_cpp2.py 的 96 组比对（同一实例、同一配比、同一权重，比较固定时域缺岗数与判环结果）
指向任意一个已编译的 C++ 仿真器（rost2 / rost3 / rost4），不改 prior/ 下的文件。
用法：python3 code/t0_validate.py <rost 可执行文件> [...]"""
import os, sys, subprocess, json, random
import numpy as np

PRIOR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "prior", "排班MOEA_小样例", "多型号多地")
sys.path.insert(0, PRIOR)
from export_inst import build            # noqa: E402
from tools_model import ToolInst         # noqa: E402
from certsim import run, start_point     # noqa: E402


def check(exe):
    """与 validate_cpp2.py 相同的随机数流和判据；返回 (组数, 不一致数)。"""
    random.seed(1); rng = np.random.default_rng(1)
    tot = bad = 0
    for spec in ["L1", "L3", "2s1", "3s1"]:
        I = build(spec); TI = ToolInst(I); N0, M0 = start_point(I, TI)
        for trial in range(12):
            Nk = [n + random.randint(0, 2) for n in N0]; Mj = [m + random.randint(0, 2) for m in M0]
            w = rng.normal(0, 1, 8).tolist() if trial > 1 else [0.1, 0, 10, 1, 0, 0, 1, 0.1]
            for Tmax, cert in [(3000, False), (60000, True)]:
                py = run(I, Nk, Mj, w, Tmax, TI, certify=cert)
                args = [exe, "eval", os.path.join(PRIOR, "inst", f"{spec}.txt"), str(Tmax), str(int(cert))] + \
                       [str(x) for x in Nk] + ["|"] + [str(x) for x in Mj] + ["|"] + [repr(x) for x in w]
                cc = json.loads(subprocess.run(args, capture_output=True, text=True).stdout)
                pc = {None: -1, False: 0, True: 1}[py["certified"]]
                if not cert:
                    ok = py["gaps"] == cc["gaps"]
                else:
                    ok = (pc == cc["cert"]) and (pc == -1 or py["period"] == cc["period"]) and (py["gaps"] == cc["gaps"]) \
                         and (pc != 1 or cc.get("exact", 1) == 1)
                tot += 1; bad += (not ok)
                if not ok: print("不一致", exe, spec, trial, Tmax, cert, py["gaps"], py["certified"], py["period"], cc)
    return tot, bad


if __name__ == "__main__":
    for exe in sys.argv[1:]:
        t, b = check(exe)
        print(json.dumps(dict(exe=os.path.basename(exe), groups=t, mismatch=b), ensure_ascii=False), flush=True)
