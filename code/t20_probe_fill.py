"""T20-0 冒烟发现大实例定不出台数后的量测（只用训练实例）：PB 与 RK 默认排序（全 0 基因）加 rescue 定台数，
记能否定出、用了哪一档起点、补救加了几个、用时。结果用来判断问题有多普遍，不用来调任何参数。
输出每行一个 JSON 到标准输出（报告里的表由 results/T20/probe/fill_m2-6.jsonl、fill_m8-20.jsonl 汇总）。
用法：python3 code/t20_probe_fill.py <m 列表，逗号分隔，例如 8,10,15,20>"""
import os, sys, json, subprocess
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names, load_train, inst_paths  # noqa: E402
from freshlb import S_max  # noqa: E402
import t13_methods as TM  # noqa: E402

R5 = os.environ.get("ROST5_BIN", os.path.join(ROOT, "code", "cpp", "rost5"))


def one(a):
    """一个 (训练实例, 方法)：rost5 size 加 rescue（6 档起点 + 最多 10 次补救），返回汇总字段。"""
    n, meth = a
    I = load_train(n); S = S_max(I); _, txt = inst_paths(n, "train")
    o = TM.opts(meth, n, S) + (["pd=rk", f"W={2 * S}"] if meth == "rk" else []) + ["rescue"]
    r = json.loads(subprocess.run([R5, "size", txt, "365", "60000", "600000"] + o, capture_output=True, text=True).stdout)
    return dict(inst=n, method=meth, R=r["R"], status=r["status"], factor=r["start_factor"], vgaps=r["verify_gaps"],
                rescue_add=r.get("rescue_add"), sec=r["sec"])


def main():
    """按给定 m 列出训练实例 × (pb, rk)，3 个并行。"""
    ms = [int(x) for x in sys.argv[1].split(",")]
    jobs = [(n, m) for n in train_names(ms) for m in ("pb", "rk")]
    with ThreadPoolExecutor(3) as ex:
        for r in ex.map(one, jobs): print(json.dumps(r, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
