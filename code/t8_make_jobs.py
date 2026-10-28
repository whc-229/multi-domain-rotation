"""T8 作业文件。
- sweep：m = 2 的 12 个实例 × {rlp, es, hand} × ρ ∈ {0.02, 0.1, 0.4} × 种子 0–2（按种子、实例交错，先拿到每个 ρ 的完整一轮）；
- full <ρ>：m ∈ {2, 3, 4, 6} 的 48 个实例 × 3 种规则 × 种子 0–4，在选定的 ρ 下；m = 2 的种子 0–2 已在扫描里做过（同名文件存在即跳过）。
  按 m 从小到大、同 m 内按种子排，先拿到 m ≤ 3（回答"能否低于 ES 饱和值"）。
用法：python3 code/t8_make_jobs.py sweep > results/T8/jobs_sweep.tsv
      python3 code/t8_make_jobs.py full 0.1 > results/T8/jobs_full.tsv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, train_names  # noqa: E402

J = os.path.join(ROOT, "results", "T8", "jobs")


def line(n, rule, rho, sd):
    """一行作业：完成标志、限时、命令。"""
    return f"{J}/{n}__{rule}__r{rho:g}__s{sd}.json\t{48 * 3600}\tcd {ROOT} && python3 code/t8_job.py {n} {rule} {rho:g} {sd}"


if __name__ == "__main__":
    if sys.argv[1] == "sweep":
        ns = sorted(train_names(m_list=[2]))
        for sd in range(3):
            for n in ns:
                for rule in ("rlp", "es", "hand"):
                    for rho in (0.02, 0.1, 0.4): print(line(n, rule, rho, sd))
    else:
        rho = float(sys.argv[2])
        for m in (2, 3, 4, 6):
            for sd in range(5):
                for n in sorted(train_names(m_list=[m])):
                    for rule in ("rlp", "es", "hand"): print(line(n, rule, rho, sd))
