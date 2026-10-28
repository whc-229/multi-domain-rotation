"""T7 / T8 的三种默认规则（不提前大修场景），返回 rost5 的规则选项列表。
1. rlp：R_LP 的 fresh_noreset_u（results/T3/rlp/<实例>__fresh_noreset_u.rlp）；
2. es：每个实例上 ES 的最好权重（两阶段做法）。默认取 curdef 配置（不提前大修场景下 ES 的饱和值就来自它）；
   lpgplan 的权重是配着提前大修计划调出来的，去掉计划后在 16 个试算实例里 14 个不如 curdef（见 T7 报告），只作对照；
3. hand：手写规则 = rost5 es 的基础个体 1（相位平衡 PB：10·flex − 相位缺额 + 0.1·等待），定台数用放大起点循环（facloop）。
用法：python3 code/t8_rules.py <实例> <rlp|es|hand>  → 打印规则选项"""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402

RDIR = os.path.join(ROOT, "results", "T3", "rlp")
HAND = {"fifo": [1, 0, 0, 0, 0, 0, 1, 0], "pb": [0.1, 0, 10, 1, 0, 0, 1, 0.1]}   # rost5 es 的基础个体 0、1


def feasible(d):
    """严格证明或 60 万格复核 0 缺岗。"""
    return d.get("R") is not None and d["R"] >= 0 and (d.get("proven") == 1 or d.get("verify_gaps") == 0)


def best_es(name, cfg="lpgplan"):
    """该实例 ES（cfg 配置）所有已有结果里资源最少的可行解（同资源取文件名排序靠前的）。
    来源：T4b 的 100 代（m ≤ 3）、T4 的 30 代、T3 的 10 代，以及 T6 用新补台规则重跑的 10 代。返回 (记录, 文件)。"""
    pats = [f"results/T4b/jobs/{name}__es100_{cfg}__s*.json", f"results/T4/jobs/{name}__esub_{cfg}__s*.json",
            f"results/T3/jobs/{name}__es_{cfg}__s*.json", f"results/T6/es/{name}__es_{cfg}__s*.json"]
    best = None
    for p in pats:
        for f in sorted(glob.glob(os.path.join(ROOT, p))):
            d = json.load(open(f))
            if feasible(d) and (best is None or d["R"] < best[0]["R"]): best = (d, f)
    return best


def rule_opts(name, rule, es_cfg="curdef", hand="pb"):
    """规则选项（传给 rost5 size / rdd）。"""
    if rule == "rlp": return [f"rlp={RDIR}/{name}__fresh_noreset_u.rlp"]
    if rule == "es":
        d, _ = best_es(name, es_cfg)
        opts = ["w=" + ",".join(repr(float(x)) for x in d["w"])]
        if es_cfg == "lpgplan": opts += [f"lpg={RDIR}/{name}.lpg.txt"]
        return opts + ["curdef"]
    if rule == "hand": return ["w=" + ",".join(repr(float(x)) for x in HAND[hand]), "facloop"]
    raise ValueError(rule)


if __name__ == "__main__":
    print(" ".join(rule_opts(sys.argv[1], sys.argv[2])))
