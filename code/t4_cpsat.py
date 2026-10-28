"""T4 上界的另一种求解：同一个周期尾整数模型（code/fresh_tail_mdrr.py 的 build()），用 OR-Tools CP-SAT 求解。
原因：HiGHS MIP 在 2 站实例上 2 小时找不到任何整数可行解（报告/T4_小实例精确.md）。CP-SAT 在等式约束很多的调度模型上找可行解通常更强。
运行环境：独立虚拟环境 .venv_cpsat（OR-Tools 依赖 protobuf，不装进共用的 conda 环境）。单线程（num_workers = 1）。
用法：.venv_cpsat/bin/python code/t4_cpsat.py <实例> <reset 0/1> <T> <P> <限时秒> <输出 json> [变量上界=80]"""
import os, sys, json, time
import numpy as np
from ortools.sat.python import cp_model   # 必须先于 highspy 导入：两者各带一份 HiGHS，后导入的会解析到错误的符号
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import load_train  # noqa: E402
from fresh_tail_mdrr import build  # noqa: E402


def solve(I, T, P, reset, tl, ub=80):
    t0 = time.time()
    M = build(I, T, P, reset)
    A = M["A"].tocsr(); n = A.shape[1]
    md = cp_model.CpModel()
    x = [md.NewIntVar(0, ub, f"x{j}") for j in range(n)]
    for r in range(A.shape[0]):
        a, b = A.indptr[r], A.indptr[r + 1]
        md.Add(sum(int(round(A.data[q])) * x[A.indices[q]] for q in range(a, b)) == int(round(M["rhs"][r])))
    md.Minimize(sum(int(round(M["cost"][j])) * x[j] for j in range(n) if M["cost"][j] != 0))
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = float(tl); s.parameters.num_workers = 1
    tb = time.time() - t0
    st = s.Solve(md)
    out = dict(T=T, P=P, reset=int(reset), Hx=M["Hx"], G=M["G"], ncol=n, nrow=A.shape[0], status=s.StatusName(st), solver="cpsat",
               build_sec=round(tb, 1), sec=round(time.time() - t0, 1), bound=s.BestObjectiveBound())
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        xv = np.array([s.Value(v) for v in x])
        out["obj"] = float(s.ObjectiveValue())
        out["N"] = [int(v) for v in xv[M["iN"]:M["iM"]]]; out["M"] = [int(v) for v in xv[M["iM"]:]]
        out["duties"] = [list(M["meta"][j][1:]) + [int(xv[j])] for j in range(n) if xv[j] > 0 and M["meta"][j][0] == "d"]
        out["z"] = [list(M["meta"][j][1:]) + [int(xv[j])] for j in range(n) if xv[j] > 0 and M["meta"][j][0] == "z"]
    else:
        out["obj"] = None
    return out


if __name__ == "__main__":
    name, rs, T, P, tl, outp = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), sys.argv[6]
    ub = int(sys.argv[7]) if len(sys.argv) > 7 else 80
    if os.path.exists(outp): print("已完成，跳过", outp); sys.exit(0)
    I = load_train(name)
    r = solve(I, T, P, bool(rs), tl, ub); r["inst"] = name; r["m"] = I.meta["m"]; r["family"] = I.meta["family"]; r["kind"] = "ub"
    with open(outp + ".tmp", "w") as f: f.write(json.dumps(r) + "\n")
    os.replace(outp + ".tmp", outp)
    print(json.dumps({k: v for k, v in r.items() if k not in ("duties", "z")}))
