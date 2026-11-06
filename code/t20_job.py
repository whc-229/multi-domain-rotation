"""T20 第五批的单个作业：一个 (阶段, 方法, 测试实例, 种子)。第五批第一次使用测试实例（11-06 用户定），参数都已由训练实例定好。
预算（单核墙钟，方法自己的预处理都算在内）：m = 2、3、4：300 秒；m = 6：600 秒；m = 8、10：1800 秒；m = 15、20：3600 秒。
所有方法共用同一个外壳（code/t13_methods.py）加 11-05 用户定的补救（rost5 的 rescue：定台数起点再加 200%、300%；起点复核不通过时
按第一次缺岗逐个加资源再判定，最多 10 次）。提前大修窗口一律 [0, 2·S_max)。
11-06 审核后：补足阶段至少补 60 次，之后只要 365 格缺岗在最近 5 次补台里刷新过最低值就继续（rost5 rescue 时）；
size 命令另加 hardtl（总限时）让定台数受作业预算约束：pb / gphh 取预算，conn_only 取预算 − 预处理，RH-MIP 的两次定台数各取预算。

方法：
- 对手（T20-2）：
  - pb：rost5 size（确定性，种子只取 0）；
  - gphh：T14 的 5 条规则，"种子" = 规则编号 0–4（results/T14_gphh/run<种子>/best.tree），rost5 size；
  - esdw、rk、d5：rost5 fd（表示 w / rk / d5），(4 + 12)，T_eval = 6·S_max，互换 0.5，stall 30；
    d5 的起点定台数限时 sizetl = 预算 / 12（每档；技术选择，见 报告/T20_0_冒烟.md）；
  - rddhand：rost5 rddfd，默认规则 PB，偏离窗口 2·S_max，resetopt，ρ = 0.1（同 T13）；
  - rhmip：code/t13_rhmip.py，第 0 步限时 = 预算 / 2，其余 7 步各 (预算 / 2) / 7，总限时 = 预算。
- 本方法与消融（T20-3；预处理见 code/t20_ours.py，按实例算一次、结果缓存，预处理用时从每个作业的预算里扣）：
  - ours：接入构造（严格窗口）+ 偏离优先的 RDD-FD，偏离窗口 ⌊0.5·S_max⌋；
  - ours_w<倍数>：窗口扫描，偏离窗口 ⌊倍数·S_max⌋（0.25、1、2），其余同 ours；
  - conn_only：只用接入构造的规则（rost5 size，确定性）；
  - noconn：σ 规则（H = 0 形式）直接从全新启动跑，不构造启动期，其余同 ours；
  - shorteval：去掉长期可行评价。rddfd 加 finite，定台数、判定、适应度都只看 360 格（720 h）；
    最后按长期口径（rost5 judge 6 万 / 60 万格，不加 finite）判定最终解，不可行时按共用外壳补救（judge 加 rescue），记补救后的资源数。
作业内核对（搜索型方法，同 T19 (b)）：最终解 rost5 judge 重判逐字段相同；抽样的 T_eval 缺岗与 rost5 eval 相同；
本方法一类另用 Python 参考实现重跑最终解（T_eval 格），缺岗、提前大修次数、窗口决策数与 rost5 eval 相同。
结果写进 results/T20/<阶段>/jobs/<方法>__<实例>__s<种子>.json（原子写入，已存在就跳过；运行时加 fcntl 锁）。
用法：python3 code/t20_job.py <阶段> <方法> <实例> <种子>"""
import os, sys, json, time, subprocess, fcntl, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_inst, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import t13_methods as TM  # noqa: E402
import t13_job  # noqa: E402
from t11_job import write_geno  # noqa: E402
import rost5_ref as ref  # noqa: E402

R5 = os.environ.get("ROST5_BIN", os.path.join(ROOT, "code", "cpp", "rost5"))     # 环境变量只用于安装前的核对
TLIM = {2: 300, 3: 300, 4: 300, 6: 600, 8: 1800, 10: 1800, 15: 3600, 20: 3600}
SPLIT = "test"
OURS = ("ours", "ours_w0.25", "ours_w1", "ours_w2", "conn_only", "noconn", "shorteval")
SEARCH_W = {"ours": 0.5, "ours_w0.25": 0.25, "ours_w1": 1.0, "ours_w2": 2.0, "noconn": 0.5, "shorteval": 0.5}
T_SHORT = 360          # 消融"去掉长期可行评价"：720 h


def c_json(args):
    """运行 rost5，返回 JSON。"""
    return json.loads(subprocess.run([R5] + [str(a) for a in args], capture_output=True, text=True, check=True).stdout)


def gphh_tree(seed):
    """GPHH 的第 seed 条规则（T14 第 seed 次训练的最好树）。"""
    return os.path.join(ROOT, "results", "T14_gphh", f"run{seed}", "best.tree")


def py_check(I, name, rec, sp):
    """本方法一类：Python 参考实现重跑最终解（T_eval 格），与 rost5 eval 比较缺岗、提前大修次数、窗口决策数。"""
    with tempfile.TemporaryDirectory(prefix="t20py_") as tmpd:
        gp = os.path.join(tmpd, "best.geno"); write_geno(gp, rec["W"], rec["ndev"], rec["qon"], rec["genes"], 1)
        _, txt = inst_paths(name, SPLIT)
        e = c_json(["eval", txt, rec["Teval"], "0"] + rec["Nk"] + ["|"] + rec["Mj"] + ["|", "|"] + sp["eval_opts"] + [f"rdd={gp}"])
    D = ref.decode_geno(rec["genes"], rec["W"], rec["ndev"], rec["qon"], True)
    p = ref.run(I, rec["Nk"], rec["Mj"], rec["Teval"], False, rdd=D, **sp["py"])
    return dict(check_py=int(p["gaps"] == e["gaps"] and p["resets"] == e["resets"] and p["ndec_w"] == e["ndec_w"]),
                check_py_detail=dict(cpp=[e["gaps"], e["resets"], e["ndec_w"]], py=[p["gaps"], p["resets"], p["ndec_w"]]))


def short_checks(rec, txt, sp, dump, tmpd, Wd):
    """消融 shorteval：作业内核对（finite 判定重判、抽样缺岗），然后按长期口径判定最终解，不可行时补救。
    结果字段：R_finite（只看 360 格的资源数）、long_feasible、R（长期口径下的最终资源数：可行时等于 R_finite，否则为补救后的）、rescue_long。"""
    out = dict(R_finite=rec["R"], finite_status=rec["status"])
    if rec["R"] < 0: return out
    gp = os.path.join(tmpd, "best.geno"); write_geno(gp, rec["W"], rec["ndev"], rec["qon"], rec["genes"], 1)
    base = ["judge", txt]; fl = [",".join(map(str, rec["Nk"])), ",".join(map(str, rec["Mj"]))]
    tail = sp["opts"] + [f"W={Wd}", "resetopt", f"rdd={gp}"]
    j = c_json(base + [T_SHORT, T_SHORT] + fl + tail + ["finite"])
    out["check_judge"] = int(j["R"] == rec["R"] and j["feasible"] == 1 and j["status"] == rec["status"] and j["verify_gaps"] == rec["verify_gaps"])
    n_ok = n_all = 0
    for line in open(dump):
        Nk, Mj, g, genes = line.rstrip("\n").split("\t")
        gq = os.path.join(tmpd, "s.geno"); write_geno(gq, rec["W"], rec["ndev"], rec["qon"], [float(x) for x in genes.split(",")], 1)
        e = c_json(["eval", txt, rec["Teval"], "0"] + json.loads(Nk) + ["|"] + json.loads(Mj) + ["|", "|"] + sp["eval_opts"] + [f"rdd={gq}"])
        n_all += 1; n_ok += int(e["gaps"] == int(g))
    out.update(check_dump_n=n_all, check_dump_same=n_ok)
    L = c_json(base + [60000, 600000] + fl + tail)
    out["long"] = L; out["long_feasible"] = L["feasible"]
    if L["feasible"] == 1:       # 长期口径可行：判定字段换成长期口径的
        rec.update(status=L["status"], verify_gaps=L["verify_gaps"], verify_cert=L["verify_cert"], verify_period=L["verify_period"],
                   exact=L["exact"], resets=L["resets"], last_reset=L["last_reset"])
        return out
    Rr = c_json(base + [60000, 600000] + fl + tail + ["rescue"])
    out["rescue_long"] = Rr
    # 最终记录改成长期口径下补救后的（补救失败时 R = −1，记"没有解"）
    rec.update(Nk=Rr["Nk"], Mj=Rr["Mj"], status=Rr["status"], verify_gaps=Rr["verify_gaps"], verify_cert=Rr["verify_cert"], verify_period=Rr["verify_period"],
               exact=Rr["exact"], resets=Rr["resets"], last_reset=Rr["last_reset"])
    rec["R"] = Rr["R"] if Rr["feasible"] == 1 else -1
    return out


def main():
    """跑一个作业，作业内核对，原子写入。"""
    stage, method, name, seed = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    out = os.path.join(ROOT, "results", "T20", stage, "jobs"); os.makedirs(out, exist_ok=True)
    jn = f"{method}__{name}__s{seed}"; path = os.path.join(out, jn + ".json")
    if os.path.exists(path): print("已完成，跳过", jn); return
    lk = open(path + ".lock", "w")
    try: fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError: print("另一个进程正在做", jn); sys.exit(4)
    if os.path.exists(path): print("已完成，跳过", jn); return
    I = load_inst(name, SPLIT); S = S_max(I); W = 2 * S; _, txt = inst_paths(name, SPLIT); m = I.meta["m"]; TL = TLIM[m]
    t0 = time.time(); pre_sec = 0.0; extra = {}
    common = ["365", "60000", "600000"]
    if method in ("pb", "gphh"):
        tree = gphh_tree(seed) if method == "gphh" else None
        opts = TM.opts(method, name, S, tree) + ["rescue", f"hardtl={TL}"]   # 11-06 审核：定台数受作业总预算约束
        rec = c_json(["size", txt] + common + opts); rec["cmd"] = " ".join([R5, "size", txt] + common + opts)
        if tree: extra["tree"] = os.path.relpath(tree, ROOT)
    elif method == "rhmip":
        d = os.path.join(ROOT, "results", "T20", stage, "rhmip"); os.makedirs(d, exist_ok=True)
        args = [sys.executable, os.path.join(ROOT, "code", "t13_rhmip.py"), name, str(TL / 2 / 7), str(TL), d, SPLIT, str(TL / 2), f"rescue,hardtl={TL}"]   # hardtl：PB 热启动与最后定台数各自最多一个预算（RH-MIP 不截断，11-06 审核）
        subprocess.run(args, capture_output=True, text=True, check=True, env=os.environ)
        rec = json.load(open(os.path.join(d, f"{name}.json"))); rec["cmd"] = " ".join(args)
    elif method in ("esdw", "rk", "d5", "rddhand"):
        opts = TM.opts(method, name, S) + ["rescue"] + ([f"sizetl={TL / 12:g}"] if method == "d5" else [])
        with tempfile.TemporaryDirectory(prefix=f"t20_{jn}_") as tmpd:
            dump = os.path.join(tmpd, "dump.tsv")
            fdo = [f"W={W}", "stall=30", "swap=0.5", "fdsel=best", f"fddump={dump}"]
            if method == "rddhand":
                args = [R5, "rddfd", txt] + common + [str(6 * S), str(10 ** 9), str(TL), "4", "12", str(seed), "0.1", "0.4"] + opts + fdo + ["resetopt"]
            else:
                args = [R5, "fd", txt] + common + [str(6 * S), str(10 ** 9), str(TL), "4", "12", str(seed), TM.SEARCH[method], "0.1", "0.4"] + opts + fdo
            r = subprocess.run(args, capture_output=True, text=True)
            try: rec = json.loads(r.stdout)
            except Exception: print("rost5 输出无法解析", r.stdout[-300:], r.stderr[-300:]); sys.exit(2)
            rec["Teval"] = rec.get("Teval", 6 * S)
            rec.update(t13_job.checks(rec, txt, opts, dump, tmpd, method, W)); rec["cmd"] = " ".join(args[:-1])
    elif method in OURS:
        import t20_ours
        pstage = stage[:-2] if stage.endswith("_e") else stage   # (e) 上线后的队尾重跑（阶段 ours_e）共用 ours 的预处理缓存
        pre = t20_ours.prepare(pstage, name)                    # 按实例缓存；pre["wall"] 是预处理用时（扣出预算）
        pre_sec = pre["wall"]; extra["pre"] = {k: pre[k] for k in pre if k not in ("thetas",)}
        t0 = time.time()                                        # 第一个作业在这里算了预处理，不重复计入作业用时
        if pre.get("sigma_R", -1) < 0: rec = dict(R=-1, status="fail", note="σ 定不出台数")
        else:
            sp = t20_ours.spec(pstage, name, I, S, method); extra["spec_fallback"] = sp["fallback"]
            if method == "conn_only":
                args = [R5, "size", txt] + common + sp["opts"] + ["rescue", f"hardtl={max(1.0, TL - pre_sec):.1f}"]
                rec = c_json(args[1:]); rec["cmd"] = " ".join(args)
            else:
                Wd = int(SEARCH_W[method] * S); tl = max(1.0, TL - pre_sec)
                short = method == "shorteval"
                hor = [str(T_SHORT)] * 4 if short else common + [str(6 * S)]
                with tempfile.TemporaryDirectory(prefix=f"t20_{jn}_") as tmpd:
                    dump = os.path.join(tmpd, "dump.tsv")
                    args = [R5, "rddfd", txt] + hor + [str(10 ** 9), f"{tl:.1f}", "4", "12", str(seed), "0.1", "0.4"] + sp["opts"] + \
                           ["rescue", f"W={Wd}", "resetopt", "stall=30", "swap=0.5", "fdsel=best", f"fddump={dump}"] + (["finite"] if short else [])
                    r = subprocess.run(args, capture_output=True, text=True)
                    try: rec = json.loads(r.stdout)
                    except Exception: print("rost5 输出无法解析", r.stdout[-300:], r.stderr[-300:]); sys.exit(2)
                    rec["Teval"] = rec.get("Teval", 6 * S)
                    if short: rec.update(short_checks(rec, txt, sp, dump, tmpd, Wd))
                    else: rec.update(t13_job.checks(rec, txt, sp["opts"], dump, tmpd, "ours", Wd))
                    rec["cmd"] = " ".join(args[:-1])
                    if rec["R"] >= 0: rec.update(py_check(I, name, rec, sp))
                extra["search_tl"] = round(tl, 1)
    else:
        raise ValueError(method)
    rec.update(extra)
    wall_job = time.time() - t0
    # 方法用时（对照预算）：rost5 自己的计时（size / fd / rddfd 的 sec，含起点定台数）+ 预处理；RH-MIP 取整个子进程的用时。不含作业内核对。
    msec = (wall_job if method == "rhmip" else (rec.get("sec") or 0.0)) + pre_sec
    rec.update(stage=stage, method=method, inst=name, m=m, family=I.meta["family"], seed=seed, S=S, W_reset=W, TL=TL,
               pre_sec=round(pre_sec, 1), wall_job=round(wall_job, 1), wall=round(wall_job + pre_sec, 1), method_sec=round(msec, 1))
    rec["proven"] = int(rec.get("verify_cert") == 1 and rec.get("exact") == 1)
    rec["feasible"] = int(rec.get("R", -1) >= 0 and (rec["proven"] == 1 or rec.get("verify_gaps") == 0))
    rec["window_ok"] = int(rec.get("last_reset", -1) < W)
    tmp = path + ".tmp"
    with open(tmp, "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    print(json.dumps({k: rec.get(k) for k in ("method", "inst", "seed", "R", "status", "feasible", "window_ok", "wall", "check_judge", "check_dump_same", "check_py",
                                              "rescue_add", "size_timeouts")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
