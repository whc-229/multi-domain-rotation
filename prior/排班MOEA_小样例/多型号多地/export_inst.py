"""把实例导出为 C++ 读取的纯文本（所有取整好的步长表都由 Python 预先算好，保证两边一致）。"""
import math, sys
from gen import make, SIX
from mtms import real_types, Inst, U
from tools_model import ToolInst, tool_lb

def build(spec):
    T = real_types()
    if spec.startswith("L"):
        names = [f"{spec}-{r}" for r in range(3)]
        I = Inst({n: T[n] for n in names}, [dict(E=names, delta=0, K=1)]); I.name = spec + "-single"
    else:
        m, sd = spec.split("s"); I = make(int(m), SIX, int(sd))
    return I

def export(I, path, tools=True):
    TI = ToolInst(I) if tools else None
    nk, ns = len(I.T), len(I.sites)
    out = []
    out.append(f"{I.name}")
    out.append(f"{nk} {ns} {TI.nt if tools else 0}")
    v0 = [math.ceil(t["off_h"] / U - 1e-9) for t in I.T]
    for k in range(nk):
        flex = len(I.elig[k])
        line = [I.d[k], I.L[k], v0[k], flex, len(I.elig[k])] + I.elig[k]
        if tools:
            R = TI.R[k]; line += [R["q"], len(R["C"])] + R["C"]
        else:
            line += [0, 0]
        out.append(" ".join(map(str, line)))
        if tools:
            R = TI.R[k]; out.append(" ".join(repr(float(x)) for x in [R["post"], R["back"]]))
        else:
            out.append("0.0 0.0")
        # 步长表：每个相位 i、每个允许站点 s：lead nx；以及 pre
        for i in range(I.L[k]):
            out.append(" ".join(f"{I.step(k, i, s)[0]} {I.step(k, i, s)[1]}" for s in I.elig[k]))
        out.append(" ".join(str(TI.pre(k, s) if tools else 0) for s in I.elig[k]))
    for s in range(ns):
        out.append(f"{repr(float(I.sites[s]['delta']))} {I.lead(s)} {I.sites[s]['K']}")
    if tools:
        for j in range(TI.nt):
            tt = TI.TT[j]; out.append(" ".join(repr(float(x)) for x in [tt["install"], tt["uninstall"], tt["maintain"]]))
        _, xN, xM = tool_lb(TI, integer=True)
        out.append(" ".join(str(int(round(v))) for v in xN)); out.append(" ".join(str(int(round(v))) for v in xM))
    else:
        from mtms import fluid_lb
        _, x = fluid_lb(I, integer=True)
        out.append(" ".join(str(int(round(v))) for v in x)); out.append("")
    open(path, "w").write("\n".join(out) + "\n")
    return I, TI

if __name__ == "__main__":
    for spec in sys.argv[1].split(","):
        I = build(spec); export(I, f"inst/{spec}.txt", tools=(len(sys.argv) < 3 or sys.argv[2] != "notools"))
        print("exported", spec, I.name)
