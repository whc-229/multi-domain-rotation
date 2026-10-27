"""对比实验：变体 × 实例 × 种子
用法: python experiment.py 变体组 种子数 代数 进程数 实例(逗号分隔|all) 输出.pkl
变体组: v1 = 第一轮（分级惩罚/仿真变异/分层精英，原解码）
        v2 = 第二轮（解码器对比）
"""
import itertools, pickle, sys, time
from multiprocessing import Pool
import roster_sim as R
import moea

GROUPS = {
    "v1": {
        "基线(死罚)":        dict(),
        "+分级惩罚":         dict(graded=True),
        "+分级+仿真变异":    dict(graded=True, sim=True),
        "+分级+分层精英":    dict(graded=True, level=True),
        "完整(A+B)":         dict(graded=True, sim=True, level=True),
    },
    "v2": {
        "原解码(死罚)":           dict(decoder="orig"),
        "原解码+分级":            dict(decoder="orig", graded=True),
        "原解码+分级+偏向初始化":  dict(decoder="orig", graded=True, seed_greedy=True),
        "阈值解码+分级":          dict(decoder="theta", graded=True),
        "排序解码+分级":          dict(decoder="rank", graded=True),
        "排序解码+分级+偏向初始化": dict(decoder="rank", graded=True, seed_greedy=True),
        "排序解码+分级+仿真变异":  dict(decoder="rank", graded=True, sim=True),
    },
    "v3": {   # 正式对比
        "原解码(死罚)":           dict(decoder="orig"),
        "原解码+分级":            dict(decoder="orig", graded=True),
        "阈值解码+分级":          dict(decoder="theta", graded=True),
        "排序解码+分级":          dict(decoder="rank", graded=True),
        "排序解码+分级+仿真变异":  dict(decoder="rank", graded=True, sim=True),
    },
    "v5": {   # 工具解码修正（允许一次新开多把同型工具）前后对比，均为 v4 中最好的完整方案
        "完整方案+原工具解码":   dict(decoder="rank", graded=True, sim=True, seed_rule=True, tool_dec="orig"),
        "完整方案+分组工具解码": dict(decoder="rank", graded=True, sim=True, seed_rule=True, tool_dec="group"),
    },
    "v4": {   # 补充：规则播种（与 v3 合并分析）
        "原解码+分级+规则播种":         dict(decoder="orig", graded=True, seed_rule=True),
        "排序解码+分级+仿真变异+规则播种": dict(decoder="rank", graded=True, sim=True, seed_rule=True),
    },
}


def job(a):
    group, inst_name, var, seed, gens = a
    kw = dict(GROUPS[group][var])
    inst = R.INSTANCES[inst_name]()
    inst.decoder = kw.pop("decoder", "orig")
    inst.tool_dec = kw.pop("tool_dec", "orig")
    t = time.time()
    r = moea.run(inst, gens=gens, seed=seed, **kw)
    r.update(inst=inst_name, var=var, seed=seed, secs=time.time() - t)
    return r


if __name__ == "__main__":
    group, seeds, gens, procs = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    insts = list(R.INSTANCES) if sys.argv[5] == "all" else sys.argv[5].split(",")
    out = sys.argv[6]
    jobs = list(itertools.product([group], insts, GROUPS[group], range(seeds), [gens]))
    t = time.time()
    with Pool(procs) as p:
        res = p.map(job, jobs, chunksize=1)
    pickle.dump(res, open(out, "wb"))
    print(f"{len(res)} runs, {time.time() - t:.0f}s")
