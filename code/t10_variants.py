"""T10 诊断实验：m = 20 上 fresh_rlp_w 失败的 3 个实例，换掉一个技术设置重跑定台数，看哪一个能让规则成功。
结果写进 results/T10/variants/<实例>__<变体>.json（原子写入，已存在就跳过）。变体：
- clip2 / clip4：lag 截断 C（站计划 CQ_{k,s} 与相位计划 CP_k）放大到 2 倍 / 4 倍（原来是 ±40 格对应的班数）；
- hs1520：定台数补足阶段的时域从 365 格加长到 1520 格（= 2H，覆盖窗口结束和 LP 时域结束）；
- why：证书阶段用旧补台规则（直接加归因信号最大的资源，certadd=why）；
- cut570：LP 计划只用到 H' = 570 = H − S_max，之后就按速率 LP 外推（改的是规则逻辑，只作观察，不作为修复）。
其余与 T6 的 fresh_rlp_w 相同（同一份 LP 解、同一份速率 LP、同一个 planf）。
用法：python3 code/t10_variants.py <实例> <变体>"""
import os, sys, json, time, subprocess, math, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from t6_rerun import frrw_npz, T6J  # noqa: E402
import t3_make_rlp as M  # noqa: E402

R5 = os.path.join(ROOT, "code", "cpp", "rost5")
SRC = os.path.join(ROOT, "results", "T6", "rules")
OUT = os.path.join(ROOT, "results", "T10", "variants")
RDIR = os.path.join(ROOT, "results", "T10", "rules")


def scale_clip(src, dst, f):
    """复制规则文件，把 p 行的截断 CQ 与 CP 行都乘 f。"""
    L = []
    for line in open(src):
        a = line.split()
        if a and a[0] == "p": a[4] = str(int(a[4]) * f); line = " ".join(a) + "\n"
        elif a and a[0] == "CP": line = "CP " + " ".join(str(int(x) * f) for x in a[1:]) + "\n"
        L.append(line)
    open(dst, "w").write("".join(L))


def make_cut(name, dst, Hc):
    """LP 计划只用到 Hc：值岗流只保留 D < Hc，规则文件的 H 写 Hc，其余（速率外推、台数起点）与 fresh_rlp_w 相同。"""
    I = load_train(name)
    rate = json.load(open(os.path.join(T6J, f"{name}__rate.json")))
    y = {}
    for key, v in rate["y_g"].items():
        k, s = map(int, key.split(",")); y.setdefault(s, {})[k] = v
    d = np.load(frrw_npz(name)); N, Mv, fl = d["N"], d["M"], d["duty_flow"]
    fl = np.array([x for x in fl if x[3] < Hc])
    ceil = lambda v: int(math.ceil(v - 1e-6))
    M.write_rlp(dst, I, Hc, False, N, [ceil(v) for v in N], [ceil(v) for v in Mv], fl, y, noscale=True)


def main():
    name, var = sys.argv[1], sys.argv[2]
    os.makedirs(OUT, exist_ok=True); os.makedirs(RDIR, exist_ok=True)
    path = os.path.join(OUT, f"{name}__{var}.json")
    if os.path.exists(path): print("已完成，跳过", name, var); return
    _, txt = inst_paths(name, "train")
    rlp = os.path.join(SRC, f"{name}__fresh_rlp_w.rlp"); planf = os.path.join(SRC, f"{name}__fresh_rlp_w.planf")
    Hs, extra = "365", []
    if var in ("clip2", "clip4"):
        dst = os.path.join(RDIR, f"{name}__{var}.rlp"); scale_clip(rlp, dst, int(var[-1])); rlp = dst
    elif var == "hs1520": Hs = "1520"
    elif var == "why": extra = ["certadd=why"]
    elif var == "cut570":
        dst = os.path.join(RDIR, f"{name}__{var}.rlp"); make_cut(name, dst, 570); rlp = dst
    else: raise ValueError(var)
    args = [R5, "rlp", txt, rlp, Hs, "60000", "600000", "planf=" + planf] + extra
    t0 = time.time()
    r = subprocess.run(args, capture_output=True, text=True, check=True)
    rec = json.loads(r.stdout); rec.update(variant=var, wall=round(time.time() - t0, 1), cmd=" ".join(args))
    rec["proven"] = int(rec.get("verify_cert") == 1 and rec.get("exact") == 1)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: rec.get(k) for k in ("inst", "variant", "R", "status", "verify_gaps", "proven", "wall")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
