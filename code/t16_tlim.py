"""定 T16（以及沿用同一限时的 T18）在 m = 8、10 上的墙钟限时（10-30 晚审核的规则）：
看 T12 在 m = 8、10 上的 RDD-FD 运行（限时 2 小时）：每个运行的"最终解" = 结束时的最好可行资源数，"找到时刻" = traj 里第一次达到它的秒数。
某个 m 上多数运行（> 一半）在 1 小时内找到最终解，就用 1 小时（3600 秒），否则用 2 小时（7200 秒）。
没有可行解的运行算"1 小时内没找到"。结果写进 results/T16/tlim.json（{"8": 秒, "10": 秒}），统计写进 results/T16/tlim_basis.json；
tlim.json 已存在时不覆盖，只打印。
用法：python3 code/t16_tlim.py"""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402


def found_at(traj):
    """traj = [[R, 代, 评价次数, 秒, 状态], …]：返回 (最终最好可行 R, 第一次达到它的秒数)；没有可行解返回 (None, None)。"""
    best = None
    for R, g, ev, t, st in traj:
        if R >= 0 and (best is None or R < best): best = R
    if best is None: return None, None
    return best, min(t for R, g, ev, t, st in traj if R == best)


def main():
    out = os.path.join(ROOT, "results", "T16", "tlim.json")
    basis = {}; tl = {}
    for m in (8, 10):
        rows = []
        for p in sorted(glob.glob(os.path.join(ROOT, "results", "T12", "jobs", f"m{m:02d}-*__rdd__s*.json"))):
            d = json.load(open(p)); R, t = found_at(d["traj"])
            rows.append(dict(job=os.path.basename(p)[:-5], R=R, t_found=t, wall=d.get("wall")))
        n = len(rows); k = sum(1 for r in rows if r["t_found"] is not None and r["t_found"] <= 3600)
        tl[str(m)] = 3600 if n > 0 and k * 2 > n else 7200
        basis[str(m)] = dict(n=n, within_1h=k, tlim=tl[str(m)], runs=rows)
        print(f"m = {m}：{n} 个运行，{k} 个在 1 小时内找到最终解 → 限时 {tl[str(m)]} 秒")
    tmp = os.path.join(ROOT, "results", "T16", "tlim_basis.json.tmp")
    with open(tmp, "w") as f: json.dump(basis, f, ensure_ascii=False, indent=1)
    os.replace(tmp, tmp[:-4])
    if os.path.exists(out): print("tlim.json 已存在，不覆盖：", open(out).read().strip()); return
    with open(out + ".tmp", "w") as f: json.dump(tl, f)
    os.replace(out + ".tmp", out)
    print("写入", out, tl)


if __name__ == "__main__":
    main()
