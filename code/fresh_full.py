"""整体 LP：直接用 prior 的 freshlb.fresh_lb 建模（不抄代码），只把求解换成 highspy。

做法：fresh_lb 在函数内部 `from scipy.optimize import linprog`，所以调用期间临时把 scipy.optimize.linprog 换成本文件的
_solve，截下它建好的 (c, A_eq, b_eq)，用 highspy 求解后按 scipy 的格式返回。fresh_lb 的其余逻辑（取 N、M、提前大修流量等）不变。
原因（T2 报告有记录）：本机 scipy 1.15.2 的 highs-ipm 在 2 站实例上 5 分钟未解完，而 highspy 可以选单纯形 / 内点并设置时限。"""
import os, sys, time, types
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mdrr  # noqa: F401,E402  （设置 prior 路径）
import scipy.optimize  # noqa: E402
import highspy  # noqa: E402
import freshlb  # noqa: E402

_cfg = dict(solver="ipm", time_limit=3600.0, keep=None)


def highs_from_eq(c, A, b, solver, time_limit):
    """min c·x，A x = b，x ≥ 0。返回 (Highs 对象, 状态字符串)。"""
    A = A.tocsc()
    h = highspy.Highs(); h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    h.setOptionValue("solver", solver)
    h.setOptionValue("time_limit", float(time_limit))
    lp = highspy.HighsLp()
    n, m = A.shape[1], A.shape[0]
    lp.num_col_ = n; lp.num_row_ = m
    lp.col_cost_ = np.asarray(c, float); lp.col_lower_ = np.zeros(n); lp.col_upper_ = np.full(n, highspy.kHighsInf)
    lp.row_lower_ = np.asarray(b, float); lp.row_upper_ = np.asarray(b, float)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A.indptr; lp.a_matrix_.index_ = A.indices; lp.a_matrix_.value_ = A.data
    lp.a_matrix_.num_col_ = n; lp.a_matrix_.num_row_ = m
    h.passModel(lp)
    h.run()
    return h, h.modelStatusToString(h.getModelStatus())


def _solve(c, A_eq=None, b_eq=None, bounds=None, method=None, options=None):
    """替换 scipy.optimize.linprog：只支持 fresh_lb 用到的形式（等式、x ≥ 0）。"""
    h, st = highs_from_eq(c, A_eq, b_eq, _cfg["solver"], _cfg["time_limit"])
    ok = st == "Optimal"
    sol = h.getSolution()
    x = np.array(sol.col_value) if ok else None
    if _cfg["keep"] is not None:
        _cfg["keep"].update(c=c, A=A_eq, b=b_eq, status=st, row_dual=np.array(sol.row_dual) if ok else None, x=x)
    return types.SimpleNamespace(x=x, fun=h.getInfo().objective_function_value if ok else None, status=0 if ok else st)


def fresh_full(I, H, reset=False, solver="ipm", time_limit=3600, keep=None, want_plan=False):
    """调用 freshlb.fresh_lb（tools=True，LP），用 highspy 求解。keep 为 dict 时存下矩阵、对偶与解。"""
    _cfg.update(solver=solver, time_limit=time_limit, keep=keep)
    old = scipy.optimize.linprog
    scipy.optimize.linprog = _solve
    try:
        t0 = time.time()
        r = freshlb.fresh_lb(I, H, tools=True, reset=reset, integer=False, time_limit=time_limit, want_plan=want_plan)
        r["sec"] = round(time.time() - t0, 1)
    finally:
        scipy.optimize.linprog = old
    return r


if __name__ == "__main__":
    # 速度测试：python3 code/fresh_full.py <prior spec 或 训练实例名> <H 倍数> <reset 0/1> <solver>
    from export_inst import build
    spec, mult, rs, solver = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    I = mdrr.load_train(spec) if spec.startswith("m") else build(spec)
    H = int(mult * freshlb.S_max(I))
    r = fresh_full(I, H, reset=bool(rs), solver=solver)
    print(dict(spec=spec, H=H, reset=rs, solver=solver, obj=r["obj"], sec=r["sec"], n=r["n"]), flush=True)
