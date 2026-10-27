"""T1 核对。
(a) 逐项核对：每个实例的 C++ 读入结果（code/cpp/inst_dump）与由 JSON 重新算出的 Python 值逐项比对
    （步长表、pre、lead、K、工具参数、起点台数等）。这是数据完整性检查，对全部 256 个文件做，只报告不一致个数。
(b) 仿真核对（与 validate_cpp*.py 同法）：只在 K1 族的训练实例上，用 prior 的 Python 仿真 certsim.run 与 C++ rost4 eval
    比对固定时域缺岗数和判环结果。rost4 / certsim 每站只有一个空闲时刻，只支持 K_s = 1，K12 族的仿真核对放到 T3（rost5）。
用法：python3 code/t1_check.py"""
import os, sys, json, math, glob, random, subprocess
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, PRIOR, INST_DIR, load_json, ToolInst, U, train_names, inst_paths  # noqa: E402
from tools_model import tool_lb   # noqa: E402
from certsim import run           # noqa: E402

DUMP = os.path.join(ROOT, "code", "cpp", "inst_dump")
ROST4 = os.path.join(PRIOR, "cpp", "rost4")


def expected(I):
    """由 Python 模型算出 C++ 应读到的全部字段。"""
    TI = ToolInst(I)
    _, xN, xM = tool_lb(TI, integer=True)
    types = []
    for k in range(len(I.T)):
        R = TI.R[k]
        types.append(dict(d=I.d[k], L=I.L[k], v0=math.ceil(I.T[k]["off_h"] / U - 1e-9), flex=len(I.elig[k]), elig=I.elig[k],
                          q=R["q"], C=R["C"], post=float(R["post"]), back=float(R["back"]),
                          step=[[list(I.step(k, i, s)) for s in I.elig[k]] for i in range(I.L[k])],
                          pre=[TI.pre(k, s) for s in I.elig[k]]))
    return dict(name=I.name, nk=len(I.T), ns=len(I.sites), nt=TI.nt, types=types,
                sites=[[float(s["delta"]), I.lead(j), s["K"]] for j, s in enumerate(I.sites)],
                tools=[[float(t["install"]), float(t["uninstall"]), float(t["maintain"])] for t in TI.TT],
                N0=[int(round(v)) for v in xN], M0=[int(round(v)) for v in xM],
                nelig=[sum(1 for k in range(len(I.T)) if s in I.elig[k]) for s in range(len(I.sites))])


def diff(a, b, path=""):
    """递归比较，返回不一致项列表。"""
    out = []
    if isinstance(a, dict):
        for key in a:
            out += diff(a[key], b.get(key), path + "." + key)
    elif isinstance(a, list):
        if not isinstance(b, list) or len(a) != len(b): return [path]
        for i, (x, y) in enumerate(zip(a, b)): out += diff(x, y, f"{path}[{i}]")
    elif a != b:
        out.append(f"{path}: {a} != {b}")
    return out


def check_fields():
    """(a) 全部实例逐项核对。"""
    n = bad = 0
    for split in ["train", "test"]:
        for jp in sorted(glob.glob(os.path.join(INST_DIR, split, "*.json"))):
            I = load_json(jp); tp = jp[:-5] + ".txt"
            cc = json.loads(subprocess.run([DUMP, tp], capture_output=True, text=True, check=True).stdout)
            d = diff(expected(I), cc); n += 1
            if d:
                bad += 1
                if split == "train": print("不一致", os.path.basename(jp), d[:3])
    return n, bad


def check_sim(trials=3):
    """(b) K1 训练实例上 certsim.run 与 rost4 eval 比对。"""
    random.seed(1); rng = np.random.default_rng(1)
    tot = bad = 0
    for name in train_names(fam_list=["K1"]):
        jp, tp = inst_paths(name, "train")
        I = load_json(jp); TI = ToolInst(I)
        _, xN, xM = tool_lb(TI, integer=True); N0 = [int(round(v)) for v in xN]; M0 = [int(round(v)) for v in xM]
        Tc = 60000 if len(I.sites) <= 4 else 20000
        for trial in range(trials):
            Nk = [n + random.randint(0, 3) for n in N0]; Mj = [m + random.randint(0, 3) for m in M0]
            w = rng.normal(0, 1, 8).tolist() if trial > 0 else [0.1, 0, 10, 1, 0, 0, 1, 0.1]
            for Tmax, cert in [(3000, False), (Tc, True)]:
                py = run(I, Nk, Mj, w, Tmax, TI, certify=cert)
                args = [ROST4, "eval", tp, str(Tmax), str(int(cert))] + [str(x) for x in Nk] + ["|"] + [str(x) for x in Mj] + ["|"] + [repr(x) for x in w]
                cc = json.loads(subprocess.run(args, capture_output=True, text=True).stdout)
                pc = {None: -1, False: 0, True: 1}[py["certified"]]
                if not cert: ok = py["gaps"] == cc["gaps"]
                else: ok = (pc == cc["cert"]) and (pc == -1 or py["period"] == cc["period"]) and py["gaps"] == cc["gaps"] and (pc != 1 or cc["exact"] == 1)
                tot += 1; bad += (not ok)
                if not ok: print("仿真不一致", name, trial, Tmax, cert, py["gaps"], py["certified"], py["period"], cc)
    return tot, bad


if __name__ == "__main__":
    n, b = check_fields(); print(json.dumps(dict(check="fields", files=n, mismatch=b)), flush=True)
    t, b2 = check_sim(); print(json.dumps(dict(check="sim_K1_train", groups=t, mismatch=b2)), flush=True)
