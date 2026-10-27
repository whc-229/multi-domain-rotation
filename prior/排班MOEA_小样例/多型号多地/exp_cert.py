"""证书适应度试点：三种适应度 × 若干带工具实例。
  finite365  ：可行 ⇔ 720h（365 格）内无缺岗（原问题设定）
  finite3000 ：可行 ⇔ 3000 格（约 8 个月）内无缺岗
  cert       ：先按 365 格筛，再要求 2 万格内无缺岗，并尽量给出状态重复证书（未证罚 0.5）
最终解统一用 20 万格的状态重复检测复核：certified=True 为严格证明永远可行，False 为出现缺岗。"""
import json, sys, time
from gen import make, SIX
from mtms import real_types, Inst
from tools_model import ToolInst
from certsim import es, run

def inst_list(spec):
    T = real_types(); out = []
    for item in spec.split(","):
        if item.startswith("L"):
            names = [f"{item}-{r}" for r in range(3)]
            out.append((item + "-单站", Inst({n: T[n] for n in names}, [dict(E=names, delta=0, K=1)])))
        else:
            m, sd = item.split("s"); I = make(int(m), SIX, int(sd)); out.append((I.name, I))
    return out

gens = int(sys.argv[2]) if len(sys.argv) > 2 else 6
MODES = [("finite", 365), ("finite", 3000), ("cert", 365)]
if len(sys.argv) > 3: MODES = [(m, int(h)) for m, h in (x.split(":") for x in sys.argv[3].split(","))]
for name, I in inst_list(sys.argv[1]):
    TI = ToolInst(I)
    for mode, H in MODES:
        t0 = time.time(); f, r, w, hist = es(I, mode, TI, H=H, gens=gens)
        v = run(I, r["Nk"], r["Mj"], w, 200000, TI, certify=True)
        rec = dict(inst=name, fitness=f"{mode}{H if mode == 'finite' else ''}", R=r["R"], Nk=r["Nk"], Mj=r["Mj"], status=r["status"],
                   verify_gaps=v["gaps"], verify_certified=v["certified"], enter=v["enter"], period=v["period"],
                   sec=round(time.time() - t0), hist=hist, w=[round(x, 3) for x in w])
        print(json.dumps(rec, ensure_ascii=False), flush=True)
