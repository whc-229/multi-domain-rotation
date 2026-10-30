"""T14 GPHH 训练：遗传规划进化派工打分规则（GP 树），只用训练实例。
训练集：m ∈ {2, 3, 4, 6, 8, 10} × 族 {K1, K12} 各取 2 个实例（idx 00、01），共 24 个；m ≥ 15 不放进训练集。
规则：机器人打分 = GP 树（特征 wait、nextocc、flex、deficit（curdef 口径）、duty、phase，运算 + − × 保护除 max min，与 rost5 的 GP 相同），
  打分小的先派；工具按保养时长（w6 = 1）；不提前大修。定台数：rost5 size（放大起点循环），与 T13 的其他规则相同。
适应度（越小越好）：训练集上 资源数 / ⌈LP_fr^Rw(4S)⌉ 的平均；
  - 训练时的判定用 6 万格判环（Tcert = Tver = 60000），判出循环（cert）按资源数算，判不出循环但无缺岗（uncert）按 资源数 + 3 算（与 ES 的 certs 模式相同），
    定不出台数或超时记 PEN = 3.0（即 3 倍下界）；
  - 竞速：实例按规模从小到大评价，部分和已超过当前种群最差个体的总和时停止，记为无穷大（不可能进入下一代）。
进化：(μ + λ) = (4 + 12)；初始个体 0 = FIFO 树 wait，1 = PB 近似树 (10·flex − deficit)，其余随机生长（深度 ≤ 4）；
  子代 0.5 概率子树交叉、0.5 概率子树变异（新子树深度 ≤ 3），深度上限 7（超了保留父代）；与 rost5 es 的 gp 模式相同的算子。
预算：代数上限 GENS 或墙钟 BUDGET 秒，先到为止；每个实例一次定台数限时 TL 秒。
输出 results/T14_gphh/run<种子>/：log.jsonl（每代最好适应度、最好树、评价次数、用时）、best.tree（最好树，全精度前缀式）、evals.jsonl（每次评价）。
中断后重跑同一条命令：从 log.jsonl 的最后一代种群续跑。
用法：python3 code/t14_gphh.py <种子> [代数上限] [墙钟预算秒] [单实例限时秒]"""
import os, sys, json, time, random, subprocess, math, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, inst_paths, train_names  # noqa: E402

R5 = os.environ.get("ROST5_BIN", os.path.join(ROOT, "code", "cpp", "rost5"))
OUT = os.path.join(ROOT, "results", "T14_gphh")
FEATS = ["wait", "nextocc", "flex", "deficit", "duty", "phase"]; OPS = ["+", "-", "*", "/", "max", "min"]
PEN = 3.0; MU, LAM = 4, 12


def train_set():
    """训练集：每个 (m, 族) 的 idx 00、01；按规模排序。"""
    names = [n for n in train_names([2, 3, 4, 6, 8, 10]) if n.endswith("-00") or n.endswith("-01")]
    return sorted(names, key=lambda n: (int(n[1:3]), n))


def lbs():
    """⌈LP_fr^Rw(4S)⌉。"""
    return {r["inst"]: math.ceil(r["LP_frRw"] - 1e-6) for r in map(json.loads, open(os.path.join(ROOT, "results", "T6_bounds.jsonl")))}


# ---------- 树：嵌套元组 (op, a, b) 或叶子（特征名 / 浮点常数） ----------
def grow(rng, d, maxd):
    """随机生长（与 rost5 GP::grow 同样的规则：到深度上限或深度 > 1 时以 0.3 概率停；叶子 0.8 概率是特征，否则 U(−2, 2) 常数）。"""
    if d >= maxd or (d > 1 and rng.random() < 0.3):
        return rng.choice(FEATS) if rng.random() < 0.8 else rng.uniform(-2, 2)
    return (rng.choice(OPS), grow(rng, d + 1, maxd), grow(rng, d + 1, maxd))


def ser(t):
    """前缀式全精度文本（rost5 tree= 读这个格式）。"""
    if isinstance(t, tuple): return f"({t[0]} {ser(t[1])} {ser(t[2])})"
    return t if isinstance(t, str) else repr(float(t))


def parse(src):
    """前缀式文本 → 树。"""
    tok = src.replace("(", " ( ").replace(")", " ) ").split(); p = [0]
    def rec():
        a = tok[p[0]]; p[0] += 1
        if a == "(":
            op = tok[p[0]]; p[0] += 1; l = rec(); r = rec(); p[0] += 1; return (op, l, r)
        return a if a in FEATS else float(a)
    return rec()


def nodes(t, path=()):
    """所有节点的路径。"""
    out = [path]
    if isinstance(t, tuple): out += nodes(t[1], path + (1,)) + nodes(t[2], path + (2,))
    return out


def get(t, p):
    """取路径 p 处的子树。"""
    for x in p: t = t[x]
    return t


def put(t, p, sub):
    """把路径 p 处换成 sub，返回新树。"""
    if not p: return sub
    l = list(t); l[p[0]] = put(t[p[0]], p[1:], sub); return tuple(l)


def depth(t):
    """深度（叶子为 1）。"""
    return 1 + max(depth(t[1]), depth(t[2])) if isinstance(t, tuple) else 1


def mutate(rng, t):
    """子树变异：随机节点换成深度 ≤ 3 的随机子树；超过深度 7 时保留原树。"""
    out = put(t, rng.choice(nodes(t)), grow(rng, 0, 3))
    return out if depth(out) <= 7 else t


def crossover(rng, a, b):
    """子树交叉：a 的随机节点换成 b 的随机子树；超过深度 7 时保留 a。"""
    out = put(a, rng.choice(nodes(a)), get(b, rng.choice(nodes(b))))
    return out if depth(out) <= 7 else a


# ---------- 评价 ----------
def size_one(tree_s, name, Tc, tl, tmpd):
    """在一个实例上用规则定台数（训练口径），返回 (资源数或 None, 状态, 用时)。"""
    _, txt = inst_paths(name, "train")
    tf = os.path.join(tmpd, "t.tree"); open(tf, "w").write(tree_s)
    t0 = time.time()
    try:
        r = subprocess.run([R5, "size", txt, "365", str(Tc), str(Tc), "w=0,0,0,0,0,0,1,0", f"tree=@{tf}", "curdef", "facloop"],
                           capture_output=True, text=True, timeout=tl)
        d = json.loads(r.stdout)
    except subprocess.TimeoutExpired:
        return None, "timeout", time.time() - t0
    if d["R"] < 0: return None, "fail", time.time() - t0
    return d["R"], d["status"], time.time() - t0


def fitness(tree, insts, LB, cap_sum, tl, tmpd, log):
    """竞速评价：返回 (平均比值或 inf, 逐实例记录)。"""
    s = ser(tree); tot = 0.0; per = []
    for n in insts:
        R, st, sec = size_one(s, n, 60000, tl, tmpd)
        ratio = PEN if R is None else (R + (3 if st == "uncert" else 0)) / LB[n]
        tot += ratio; per.append((n, R, st, round(sec, 1)))
        if tot > cap_sum: log.write(json.dumps(dict(tree=s, raced_out=True, per=per)) + "\n"); log.flush(); return math.inf, per
    log.write(json.dumps(dict(tree=s, fit=tot / len(insts), per=per)) + "\n"); log.flush()
    return tot / len(insts), per


def main():
    seed = int(sys.argv[1]); gens = int(sys.argv[2]) if len(sys.argv) > 2 else 10 ** 6
    budget = float(sys.argv[3]) if len(sys.argv) > 3 else 86400.0; tl = float(sys.argv[4]) if len(sys.argv) > 4 else 300.0
    d = os.path.join(OUT, f"run{seed}"); os.makedirs(d, exist_ok=True)
    insts = train_set(); LB = lbs(); rng = random.Random(seed)
    tmpd = tempfile.mkdtemp(prefix=f"gphh{seed}_")
    lp = os.path.join(d, "log.jsonl"); ev = open(os.path.join(d, "evals.jsonl"), "a")
    t_used = 0.0; g0 = 0; pop = None
    if os.path.exists(lp):     # 续跑：取最后一代的种群与随机数状态
        last = [json.loads(l) for l in open(lp)][-1]
        pop = [(parse(x["tree"]), x["fit"]) for x in last["pop"]]; g0 = last["gen"] + 1; t_used = last["sec"]
        rng.setstate(eval(last["rng"]))   # noqa: S307  只读自己写的文件
    t0 = time.time() - t_used
    if pop is None:
        init = ["wait", ("-", ("*", 10.0, "flex"), "deficit")] + [grow(rng, 0, 4) for _ in range(MU + LAM - 2)]
        pop = []
        for t in init:
            f, _ = fitness(t, insts, LB, math.inf, tl, tmpd, ev); pop.append((t, f))
        pop.sort(key=lambda x: x[1]); pop = pop[:MU]
        with open(lp, "a") as f:
            f.write(json.dumps(dict(gen=0, best=pop[0][1], best_tree=ser(pop[0][0]), pop=[dict(tree=ser(t), fit=v) for t, v in pop],
                                    sec=round(time.time() - t0, 1), rng=repr(rng.getstate()))) + "\n")
        g0 = 1
    g = g0
    while g < gens and time.time() - t0 < budget:
        kids = []
        cap = sum(1 for _ in insts) * pop[-1][1] if math.isfinite(pop[-1][1]) else math.inf
        for _ in range(LAM):
            a = pop[rng.randrange(MU)][0]
            t = crossover(rng, a, pop[rng.randrange(MU)][0]) if rng.random() < 0.5 else mutate(rng, a)
            f, _ = fitness(t, insts, LB, cap, tl, tmpd, ev); kids.append((t, f))
        pop = sorted(kids + pop, key=lambda x: x[1])[:MU]     # 同适应度时子代在前
        with open(lp, "a") as f:
            f.write(json.dumps(dict(gen=g, best=pop[0][1], best_tree=ser(pop[0][0]), pop=[dict(tree=ser(t), fit=v) for t, v in pop],
                                    sec=round(time.time() - t0, 1), rng=repr(rng.getstate()))) + "\n")
        g += 1
    with open(os.path.join(d, "best.tree.tmp"), "w") as f: f.write(ser(pop[0][0]) + "\n")
    os.replace(os.path.join(d, "best.tree.tmp"), os.path.join(d, "best.tree"))
    print(json.dumps(dict(seed=seed, gens=g, best=pop[0][1], tree=ser(pop[0][0]), sec=round(time.time() - t0, 1))))


if __name__ == "__main__":
    main()
