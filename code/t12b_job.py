"""T12b：严格窗口下重算 fresh_rlp_w 在 m ≤ 6 上的结果（m ≥ 8 已在 T12 的 base 作业里算过）。
设置与 T12 的 base 作业相同：rost5 size（365 格补足时域、60000 格判环、60 万格复核）+ fresh_rlp_w 规则文件与 planf（results/T6/rules/）
+ planW = 2·S_max（计划提前大修只在窗口内执行）+ hsfb = 8·S_max（T10 的后备补足时域）。
逐实例写 results/T12b/jobs/<实例>__base.json（原子写入，已存在就跳过）；一个进程顺序跑完 m ∈ {2, 3, 4, 6} 的 48 个训练实例。
用法：python3 code/t12b_job.py"""
import os, sys, json, time, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, load_train, inst_paths, train_names  # noqa: E402
from freshlb import S_max  # noqa: E402
from t12_job import rule_opts, R5  # noqa: E402

OUT = os.path.join(ROOT, "results", "T12b", "jobs")


def main():
    """逐实例定台数并复核。"""
    os.makedirs(OUT, exist_ok=True)
    for n in sorted(train_names([2, 3, 4, 6]), key=lambda n: (int(n[1:3]), n)):
        path = os.path.join(OUT, f"{n}__base.json")
        if os.path.exists(path): continue
        I = load_train(n); _, txt = inst_paths(n, "train"); S = S_max(I)
        args = [R5, "size", txt, "365", "60000", "600000"] + rule_opts(n, S)
        t0 = time.time()
        rec = json.loads(subprocess.run(args, capture_output=True, text=True, check=True).stdout)
        rec.update(inst=n, m=I.meta["m"], family=I.meta["family"], kind="base", S=S, W=2 * S, wall=round(time.time() - t0, 1), cmd=" ".join(args))
        rec["proven"] = int(rec.get("verify_cert") == 1 and rec.get("exact") == 1)
        rec["window_ok"] = int(rec.get("last_reset", -1) < 2 * S)
        with open(path + ".tmp", "w") as f: f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        os.replace(path + ".tmp", path)
        print(n, rec["R"], rec["status"], "last_reset", rec.get("last_reset"), "W", 2 * S, f"{rec['wall']} s", flush=True)


if __name__ == "__main__":
    main()
