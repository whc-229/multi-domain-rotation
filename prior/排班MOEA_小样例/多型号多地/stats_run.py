"""统计实验驱动：实例列表 × 适应度 × 进化种子，调用 C++ 版 rost，逐行输出 JSON。
适应度：finite365（原问题 720h）、finite3000（约 8 个月）、cert（365 格定台数 + 6 万格内判环证书）。
复核：60 万格（Brent 判环）。用法：python3 stats_run.py <实例spec,...> <种子数> <代数> <输出文件> [gp]"""
import subprocess, sys, json, os
from export_inst import build, export
specs = sys.argv[1].split(","); nseed = int(sys.argv[2]); gens = sys.argv[3]; out = sys.argv[4]
gp = len(sys.argv) > 5 and sys.argv[5] == "gp"
MODES = [("finite", 365), ("finite", 3000), ("cert", 365)]
if len(sys.argv) > 6: MODES = [(m, int(h)) for m, h in (x.split(":") for x in sys.argv[6].split(","))]
os.makedirs("inst", exist_ok=True)
done = set()
if os.path.exists(out):
    for l in open(out):
        d = json.loads(l); done.add((d["spec"], d["mode"], d["H"], d["seed"], d["gp"]))
for spec in specs:
    path = f"inst/{spec}.txt"
    if not os.path.exists(path): export(build(spec), path)
    for seed in range(nseed):
        for mode, H in MODES:
            if (spec, mode, H, seed, int(gp)) in done: continue
            args = [os.environ.get("ROST", "cpp/rost2"), "es", path, mode, str(H), gens, "4", "12", str(seed), "60000", "600000"] + (["gp"] if gp else [])
            r = subprocess.run(args, capture_output=True, text=True)
            try:
                d = json.loads(r.stdout)
            except Exception:
                d = dict(error=r.stderr[-300:], inst=spec)
            d["spec"] = spec; d["mode"] = mode; d["H"] = H; d["seed"] = seed; d["gp"] = int(gp)
            with open(out, "a") as f: f.write(json.dumps(d, ensure_ascii=False) + "\n")
