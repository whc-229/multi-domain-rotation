"""rost4 统计驱动：实例 × 配置 × 种子。配置写成 "名字=选项1+选项2"，例如 curdef=curdef、reset=reset、resetTw=reset+Tw=3000。
用法：python3 stats_run4.py <实例spec,...> <种子数> <代数> <输出.jsonl> <配置,...> [mode=certs] [H=365]"""
import subprocess, sys, json, os
from export_inst import build, export
specs = sys.argv[1].split(","); nseed = int(sys.argv[2]); gens = sys.argv[3]; out = sys.argv[4]
cfgs = [c.split("=", 1) for c in sys.argv[5].split(",")]
mode = sys.argv[6] if len(sys.argv) > 6 else "certs"; H = sys.argv[7] if len(sys.argv) > 7 else "365"
os.makedirs("inst", exist_ok=True)
done = set()
if os.path.exists(out):
    for l in open(out):
        d = json.loads(l); done.add((d["spec"], d["cfg"], d["seed"]))
for spec in specs:
    path = f"inst/{spec}.txt"
    if not os.path.exists(path): export(build(spec), path)
    for seed in range(nseed):
        for name, opts in cfgs:
            if (spec, name, seed) in done: continue
            args = ["cpp/rost4", "es", path, mode, H, gens, "4", "12", str(seed), "60000", "600000"] + [o for o in opts.split("+") if o]
            r = subprocess.run(args, capture_output=True, text=True)
            try: d = json.loads(r.stdout)
            except Exception: d = dict(error=r.stderr[-300:], inst=spec)
            d["spec"] = spec; d["cfg"] = name; d["seed"] = seed
            with open(out, "a") as f: f.write(json.dumps(d, ensure_ascii=False) + "\n")
