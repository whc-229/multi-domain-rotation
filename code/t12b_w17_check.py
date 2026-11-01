"""T12b 补充核对（rost5 改动：es 输出多一个 w17 字段 = 全精度权重；原来的 w 只保留 4 位小数）。
起因：T12b 的 216 个严格窗口 ES 作业里，有 9 个用输出的 w（4 位小数）重新判定，结果与 es 自己的判定不同
（缺岗数不同，其中 3 个 es 判可行、重新判定有缺岗）。es 内部用的是全精度权重，舍入后决策的并列与次序会变。
这里核对：
(1) 新旧程序的 es 输出除 w17 与用时外逐字段相同（仿真逻辑没动）：实例 m02-K1-03、m03-K12-00、m04-K1-01，种子 0–2，
    2 代，参数与 T12b 相同（含 planW = 2·S_max）；
(2) 用 w17 重新判定（rost5 eval，同一组选项 + planW，60 万格判环）：缺岗、判环结果、周期与 es 的输出相同；最后一次提前大修 < 2·S_max；
(3) 另记下用 4 位小数的 w 重新判定是否相同（只记录，不算不一致）。
输出 results/T12b/check_w17.jsonl，打印不一致数。用法：python3 code/t12b_w17_check.py <旧 rost5> <新 rost5>"""
import os, sys, json, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402

T3R = os.path.join(ROOT, "results", "T3", "rlp")


def c_json(r5, args):
    """运行 rost5，解析一行 JSON。"""
    return json.loads(subprocess.run([r5] + [str(a) for a in args], capture_output=True, text=True, check=True).stdout)


def es_opts(n):
    """T6 / T12b es_lpgplan 的规则选项。"""
    return [f"lpg={T3R}/{n}.lpg.txt", f"planf={T3R}/{n}__fresh_rlp.planf"]


def rejudge(r5, txt, e, w, n, W):
    """用给定权重重新判定 es 的最终解：返回 (60 万格判环结果, 3000 格的最后一次提前大修)。"""
    ev = e["Nk"] + ["|"] + e["Mj"] + ["|"] + w + ["|"] + es_opts(n) + ["curdef", f"planW={W}"]
    return c_json(r5, ["eval", txt, 600000, 1] + ev), c_json(r5, ["eval", txt, 3000, 0] + ev)["last_reset"]


def main():
    """逐个 (实例, 种子) 跑新旧 es 并重新判定，写核对结果。"""
    old, new = sys.argv[1], sys.argv[2]
    out = []; bad = 0
    for n in ("m02-K1-03", "m03-K12-00", "m04-K1-01"):
        I = load_train(n); W = 2 * S_max(I); _, txt = inst_paths(n, "train")
        for seed in (0, 1, 2):
            args = ["es", txt, "certs", 365, 2, 4, 12, seed, 60000, 600000] + es_opts(n) + [f"planW={W}"]
            a = c_json(old, args); b = c_json(new, args)
            w17 = b.pop("w17"); a.pop("sec", None); b.pop("sec", None)
            same = a == b
            row = dict(inst=n, seed=seed, R=b["R"], same_old_new=int(same))
            ok = same
            if b["R"] >= 0:
                v, lr = rejudge(new, txt, b, w17, n, W)
                rep = v["gaps"] == b["verify_gaps"] and v["cert"] == b["verify_cert"] and (b["verify_cert"] != 1 or v["period"] == b["verify_period"])
                v4, _ = rejudge(new, txt, b, b["w"], n, W)
                rep4 = v4["gaps"] == b["verify_gaps"] and v4["cert"] == b["verify_cert"]
                row.update(rejudge_w17=int(rep), last_reset=lr, W=W, rejudge_w4=int(rep4), gaps_es=b["verify_gaps"], gaps_w17=v["gaps"], gaps_w4=v4["gaps"])
                ok = ok and rep and lr < W
            bad += (not ok); row["ok"] = int(ok); out.append(row)
            print(json.dumps(row, ensure_ascii=False), flush=True)
    with open(os.path.join(ROOT, "results", "T12b", "check_w17.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(json.dumps(dict(cases=len(out), mismatch=bad)))


if __name__ == "__main__":
    main()
