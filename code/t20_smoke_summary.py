"""T20-0 冒烟汇总：只确认能跑与四处改动生效，不看资源数（11-06 要求）。
读 results/T20/smoke/jobs/*.json、results/T20/smoke/pre/*.json、results/T20/smoke/rhmip/*.json、results/T20/lb/jobs/（冒烟里顺带跑的 3 个 LP）。
逐作业列：可行与否、状态、作业用时（含预处理）/ 预算、起点放大档、补救加的个数、定台数超时档数、作业内核对、提前大修是否都在 [0, 2·S_max) 内；
另列 RH-MIP 各步限时与用时、接入构造的窗口与额外机器人、LP 是否解出。用法：python3 code/t20_smoke_summary.py"""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402

D = os.path.join(ROOT, "results", "T20", "smoke")


def yn(x):
    """1/0/None → 是/否/—。"""
    return "—" if x is None else ("是" if x else "否")


def main():
    """打印 markdown 表。"""
    rows = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(D, "jobs", "*.json")))]
    L = ["| 实例 | 方法 | 可行 | 状态 | 用时 / 预算（秒） | 起点档 | 补救加 | 超时档 | judge 重判 | 抽样 | Python 重跑 | 提前大修在窗口内 |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in sorted(rows, key=lambda r: (r["inst"], r["method"])):
        fac = r.get("start_factor", r.get("base_factor"))
        L.append(f"| {r['inst']} | {r['method']} | {yn(r['feasible'])} | {r.get('status')} | {r['wall']:.0f} / {r['TL']} | {fac} | {r.get('rescue_add')} | {r.get('size_timeouts')} | "
                 f"{yn(r.get('check_judge'))} | {str(r.get('check_dump_same')) + '/' + str(r.get('check_dump_n')) if r.get('check_dump_n') is not None else '—'} | "
                 f"{yn(r.get('check_py'))} | {yn(r['window_ok'])} |")
    print("\n".join(L))
    print("\n### RH-MIP 各步")
    for p in sorted(glob.glob(os.path.join(D, "rhmip", "*.json"))):
        r = json.load(open(p))
        st = "；".join(f"第 {s['step']} 步 {s.get('status')} {s.get('sec')} 秒" for s in r["steps"])
        print(f"- {r['inst']}：第 0 步限时 {r.get('tl_step0')}，其余各步 {r['tl_step']:.1f}，总限时 {r['tl_total']}；MIP 用时 {r['sec_mip']} 秒，含定台数 {r['wall']} 秒。{st}")
    print("\n### 本方法的预处理（接入构造，窗口 [0, 2·S_max)）")
    for p in sorted(glob.glob(os.path.join(D, "pre", "*.json"))):
        r = json.load(open(p)); c = r.get("conn") or {}
        late = None
        if c.get("plan"):     # 计划里"本班后提前大修"的条目落在窗口外的个数（应为 0）
            rows = [ln.split() for ln in open(os.path.join(ROOT, c["plan"])).readlines()[1:]]
            late = sum(1 for a in rows if int(a[4]) == 1 and int(a[1]) >= r["W_reset"])
        print(f"- {r['inst']}：σ 定台数 {r['sigma_status']}（起点档 {r['sigma_factor']}，补救加 {r['sigma_rescue_add']}），σ 可行 {yn(r['sigma_feasible'])}；"
              f"构造 {r.get('n_variants')} 个、可行 {r.get('n_feasible')} 个；选中 θ{c.get('theta')} {c.get('mode')}，额外机器人 {c.get('n_ext')}、预工作 {c.get('n_pre')}、"
              f"计划提前大修 {c.get('n_reset')}（落在窗口外 {late} 个）、T_plan / S_max = {c.get('Tplan', 0) / r['S']:.2f}；预处理用时 {r['wall']} 秒。")
    print("\n### LP")
    for p in sorted(glob.glob(os.path.join(ROOT, "results", "T20", "lb", "jobs", "*.json"))):
        r = json.load(open(p))
        ok = r.get("R_rate_g") is not None if r["kind"] == "rate" else r.get("obj") is not None
        print(f"- {os.path.basename(p)[:-5]}：解出 {yn(ok)}，用时 {r['wall']} 秒")


if __name__ == "__main__":
    main()
