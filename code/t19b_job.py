"""T19 (b) 单个作业：接入构造作默认规则 + RDD-FD（rost5 rddfd），一个 (训练实例, 种子)。只用训练实例。
默认规则 = (a) 里该实例资源最少的接入构造（results/T19_connector.jsonl 的 steady 行：θ、入场口径）：
  稳态松弛规则 σ（H = 0、全新启动）+ 计划文件 planx（[0, T_plan) 按计划派，加 planxlate、记账相位列）+ planW = max(2·S_max, T_plan)。
  规则文件的起点台数 N0 / M0 改成构造台数 N_σ + ΔN、M_σ（写到 results/T19/rules/<实例>__conn.rlp），定台数从这里出发。
RDD-FD 与 T18 的 W = 0.5·S_max 相同：偏离窗口 W = ⌊0.5·S_max⌋，ρ = 0.1、q_on = 0.4、(4 + 12)、T_eval = 6·S_max、互换 0.5、stall 30、resetopt，
限时同 T18（results/T16/tlim.json；m = 4、6 为 300、600 秒）。偏离优先（rost5 的 devfirst，11-05 审核）：
  偏离窗口内先按计划派，再按偏离改派；改派到别的机器人时这一班的计划条目作废、计划的本班后提前大修取消，之后对不上的条目跳过、由规则补。
作业内核对：
  - 最终解用 rost5 judge 重新判定，逐字段相同；
  - 抽样的 T_eval 缺岗与 rost5 eval 相同（同 T11 / T18）；
  - 最终解用 Python 参考实现在最终台数下跑 T_eval 格：缺岗、提前大修次数、窗口决策数与 rost5 eval 相同。
结果写进 results/T19/b/<实例>__s<种子>.json（原子写入，已存在就跳过；运行时加 fcntl 锁）。
用法：python3 code/t19b_job.py <实例> <种子>"""
import os, sys, json, time, subprocess, fcntl, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
from t11_job import write_geno  # noqa: E402
from t16_job import TLIM  # noqa: E402
import rost5_ref as ref  # noqa: E402

R5 = os.environ.get("ROST5_BIN", os.path.join(ROOT, "code", "cpp", "rost5"))     # 环境变量只用于开发期试算
OUT = os.environ.get("T19B_OUT", os.path.join(ROOT, "results", "T19", "b"))
MULT = 0.5


def conn_spec(name):
    """(a) 里该实例最好的接入构造：返回 dict(I, txt, S, rule（台数改成构造台数的规则文件）, plan, planW, Nk, Mj, opts（rost5 规则选项）, py（Python 参数）)。"""
    row = [r for r in map(json.loads, open(os.path.join(ROOT, "results", "T19_connector.jsonl"))) if r["inst"] == name and r["kind"] == "steady"][0]
    c = row["conn"]
    I = load_train(name); S = S_max(I); _, txt = inst_paths(name, "train")
    plan = os.path.join(ROOT, "results", "T19", "plans", f"steady__{name}__th{c['theta']}__{c['mode']}.planx")
    src = os.path.join(ROOT, "results", "T19", "rules", f"{name}__steady__h0fresh.rlp")
    rule = os.path.join(ROOT, "results", "T19", "rules", f"{name}__conn.rlp")
    Nk = [row["N_sigma"][k] + c["dN"][k] for k in range(len(row["N_sigma"]))]; Mj = list(row["M_sigma"])
    if not os.path.exists(rule):
        lines = []
        for ln in open(src):
            a = ln.split()
            if a and a[0] == "N0": ln = "N0 " + " ".join(map(str, Nk)) + "\n"
            elif a and a[0] == "M0": ln = "M0 " + " ".join(map(str, Mj)) + "\n"
            lines.append(ln)
        tmp = rule + f".{os.getpid()}.tmp"; open(tmp, "w").writelines(lines); os.replace(tmp, rule)
    opts = [f"rlp={rule}", f"planx={plan}", "planxlate", f"planW={c['planW']}", "devfirst"]
    X = ref.parse_planx(plan); X["late"] = True
    py = dict(rlp=ref.parse_rlp(I, rule), planx=X, plan_until=c["planW"], devfirst=True)
    return dict(I=I, txt=txt, S=S, row=row, rule=rule, plan=plan, planW=c["planW"], Nk=Nk, Mj=Mj, opts=opts, py=py)


def c_json(args):
    """运行 rost5，返回 JSON。"""
    return json.loads(subprocess.run([R5] + [str(a) for a in args], capture_output=True, text=True, check=True).stdout)


def main():
    """跑 rddfd、作业内核对，原子写入。"""
    name, seed = sys.argv[1], int(sys.argv[2])
    os.makedirs(OUT, exist_ok=True)
    jn = f"{name}__s{seed}"; path = os.path.join(OUT, jn + ".json")
    if os.path.exists(path): print("已完成，跳过", jn); return
    lk = open(path + ".lock", "w")
    try: fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError: print("另一个进程正在做", jn); sys.exit(4)
    if os.path.exists(path): print("已完成，跳过", jn); return
    sp = conn_spec(name); I, txt, S = sp["I"], sp["txt"], sp["S"]; m = I.meta["m"]; W = int(MULT * S)
    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix=f"t19b_{jn}_") as tmpd:
        dump = os.path.join(tmpd, "dump.tsv")
        args = [R5, "rddfd", txt, "365", "60000", "600000", str(6 * S), str(10 ** 9), str(TLIM[m]), "4", "12", str(seed), "0.1", "0.4"] + sp["opts"] + \
               [f"W={W}", "resetopt", "stall=30", "swap=0.5", "fdsel=best", f"fddump={dump}"]
        r = subprocess.run(args, capture_output=True, text=True)
        try: rec = json.loads(r.stdout)
        except Exception: print("rost5 输出无法解析", r.stdout[-300:], r.stderr[-300:]); sys.exit(2)
        rec.update(inst=name, m=m, family=I.meta["family"], mult=MULT, S=S, planW=sp["planW"], conn=sp["row"]["conn"], conn_R=sp["row"]["conn"]["R"],
                   wall=round(time.time() - t0, 1), cmd=" ".join(args[:-1]))
        if rec["R"] >= 0:
            gp = os.path.join(tmpd, "best.geno"); write_geno(gp, rec["W"], rec["ndev"], rec["qon"], rec["genes"], 1)
            j = c_json(["judge", txt, "60000", "600000", ",".join(map(str, rec["Nk"])), ",".join(map(str, rec["Mj"]))] + sp["opts"] + [f"W={rec['W']}", "resetopt", f"rdd={gp}"])
            rec["check_judge"] = int(j["R"] == rec["R"] and j["feasible"] == 1 and j["status"] == rec["status"] and j["verify_gaps"] == rec["verify_gaps"]
                                     and j["verify_cert"] == rec["verify_cert"] and j["verify_period"] == rec["verify_period"] and j["exact"] == rec["exact"]
                                     and j["resets"] == rec["resets"] and j["last_reset"] == rec["last_reset"])
            # 参考实现重跑最终解（T_eval 格）
            e = c_json(["eval", txt, rec["Teval"], "0"] + rec["Nk"] + ["|"] + rec["Mj"] + ["|", "|"] + sp["opts"] + [f"rdd={gp}"])
            D = ref.decode_geno(rec["genes"], rec["W"], rec["ndev"], rec["qon"], True)
            p = ref.run(I, rec["Nk"], rec["Mj"], rec["Teval"], False, rdd=D, **sp["py"])
            rec["check_py"] = int(p["gaps"] == e["gaps"] and p["resets"] == e["resets"] and p["ndec_w"] == e["ndec_w"])
            rec["check_py_detail"] = dict(cpp=[e["gaps"], e["resets"], e["ndec_w"]], py=[p["gaps"], p["resets"], p["ndec_w"]])
        n_ok = n_all = 0
        for line in open(dump):
            Nk, Mj, g, genes = line.rstrip("\n").split("\t")
            gq = os.path.join(tmpd, "s.geno"); write_geno(gq, rec["W"], rec["ndev"], rec["qon"], [float(x) for x in genes.split(",")], 1)
            e = c_json(["eval", txt, rec["Teval"], "0"] + json.loads(Nk) + ["|"] + json.loads(Mj) + ["|", "|"] + sp["opts"] + [f"rdd={gq}"])
            n_all += 1; n_ok += int(e["gaps"] == int(g))
        rec.update(check_dump_n=n_all, check_dump_same=n_ok)
    rec["proven"] = int(rec.get("verify_cert") == 1 and rec.get("exact") == 1)
    rec["feasible"] = int(rec.get("R", -1) >= 0 and (rec["proven"] == 1 or rec.get("verify_gaps") == 0))
    rec["window_ok"] = int(rec.get("last_reset", -1) < sp["planW"])
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: rec.get(k) for k in ("inst", "seed", "base_R", "conn_R", "R", "status", "feasible", "n_active", "wall", "check_judge", "check_dump_same", "check_py")},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
