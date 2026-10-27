"""T1：MDRR 实例生成与训练 / 测试划分（划分在生成时固定，写死在本文件里）。

以 prior/.../gen.py 的 make() 为基础：
- 机型取 SIX（L1-0..2、L3-0..2）；每站随机限定 2–3 种；额外单程路程 δ_s ∈ {0, 1, 2, 4} h；共享工具池（ToolInst）。
- 站数 m ∈ {2, 3, 4, 6, 8, 10, 15, 20}；族 K1（K_s = 1）、K12（K_s 在 {1, 2} 中均匀抽取）。
- 每个 (m, 族) 16 个实例：序号 0–5 为训练，6–15 为测试。
- 随机种子 = md5("MDRR-v1|m|族|序号") 的前 8 位十六进制，只由 (m, 族, 序号) 决定。
抽样顺序（逐站）：先抽允许机型（与 make() 相同的 rng.sample / rng.randint），再抽 δ_s，K12 族再抽 K_s。
输出 bench/instances/{train,test}/<名>.json 与 .txt（C++ 文本，沿用 export_inst.export 的格式，每站一行 "δ lead K"）、
bench/split.json。用法：python3 code/gen_mdrr.py"""
import os, sys, json, random, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import Inst, real_types, tool_params, robot_tool_params, BENCH, inst_paths  # noqa: E402
from export_inst import export   # noqa: E402  （prior，C++ 文本格式）
from gen import SIX              # noqa: E402

M_LIST = [2, 3, 4, 6, 8, 10, 15, 20]
FAMILIES = ["K1", "K12"]
N_TRAIN, N_TEST = 6, 10
DELTAS = (0, 1, 2, 4)


def split_of(idx):
    """固定划分：序号 < 6 为训练，其余为测试。"""
    return "train" if idx < N_TRAIN else "test"


def seed_of(m, fam, idx):
    """由 (m, 族, 序号) 决定的种子。"""
    return int(hashlib.md5(f"MDRR-v1|{m}|{fam}|{idx}".encode()).hexdigest()[:8], 16)


def make_mdrr(m, fam, idx):
    """生成一个实例；返回 (Inst, 元数据)。"""
    seed = seed_of(m, fam, idx)
    rng = random.Random(seed); T = real_types()
    sites = []
    for _ in range(m):
        E = rng.sample(SIX, rng.randint(2, 3))
        delta = rng.choice(DELTAS)
        K = 1 if fam == "K1" else rng.choice((1, 2))
        sites.append(dict(E=sorted(E), delta=delta, K=K))
    used = sorted({t for s in sites for t in s["E"]})
    name = f"m{m:02d}-{fam}-{idx:02d}"
    I = Inst({t: T[t] for t in used}, sites, cost={t: 1.0 for t in used}, name=name)
    meta = dict(name=name, m=m, family=fam, idx=idx, seed=seed, split=split_of(idx))
    return I, meta


def to_json(I, meta):
    """实例的完整 JSON（含工具参数，便于独立复现）。键排序、固定缩进，保证字节级可复现。"""
    RP = robot_tool_params()
    d = dict(meta)
    d.update(types={n: t for n, t in zip(I.tn, I.T)}, sites=I.sites, cost=I.cost,
             tools=tool_params(), robot_tools={n: RP[n] for n in I.tn}, grid_h=2.0)
    return json.dumps(d, ensure_ascii=False, sort_keys=True, indent=1) + "\n"


def main():
    entries = []
    for m in M_LIST:
        for fam in FAMILIES:
            for idx in range(N_TRAIN + N_TEST):
                I, meta = make_mdrr(m, fam, idx)
                jp, tp = inst_paths(meta["name"], meta["split"])
                os.makedirs(os.path.dirname(jp), exist_ok=True)
                with open(jp, "w") as f: f.write(to_json(I, meta))
                export(I, tp, tools=True)
                entries.append(meta)
    with open(os.path.join(BENCH, "split.json"), "w") as f:
        f.write(json.dumps(dict(version="MDRR-v1", rule="idx<6 train, else test", instances=entries),
                           ensure_ascii=False, indent=1) + "\n")
    print("生成", len(entries), "个实例；训练", sum(e["split"] == "train" for e in entries),
          "测试", sum(e["split"] == "test" for e in entries))


if __name__ == "__main__":
    main()
