"""T20 无人值守总队列的作业表（11-06 晚，TASKS"无人值守安排"）。每段一个文件，由 code/t20_master.sh 按顺序跑：
  q1：T20-1 剩余的 LP（就是 results/T20/lb/jobs.tsv，不另写）；
  q2：T20-2 m ≤ 6 的全部对手（pb、gphh、esdw、rk、d5、rddhand、rhmip）；
  q3：T20-2 m ≥ 8 的确定性方法：pb、gphh（5 条规则）、rhmip（1 次），以及 T20-3 的只用接入构造（conn_only，阶段目录仍是 ours）；
      RH-MIP 不是搜索型方法，(e) 不涉及它，所以放在这一段（技术选择，见报告）；
  q4：T20-3（阶段 ours）：先 m ≤ 6 再 m ≥ 8；本方法（种子 0–4）、conn_only（m ≤ 6，m ≥ 8 已在 q3）、消融 noconn / shorteval（m = 4–10，种子 0–2）、
      窗口扫描 ours_w0.25 / ours_w1 / ours_w2（m = 4、6，种子 0–2）；
  q5：T20-2 m ≥ 8 的搜索型对手（esdw、rk、d5、rddhand，种子 0–4）。
每段内按预算从大到小（同 code/t20_make_jobs.py）。行格式同 code/t20_run.sh。
用法：python3 code/t20_make_queue.py → results/T20/queue/q2.tsv … q5.tsv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT, test_names  # noqa: E402
from t20_make_jobs import job_lines  # noqa: E402

SMALL, LARGE = [2, 3, 4, 6], [8, 10, 15, 20]
SEARCH = ("esdw", "rk", "d5", "rddhand")


def ours_lines(ms):
    """T20-3 在给定 m 上的作业行。"""
    allt = test_names(ms); mid = test_names([m for m in ms if m in (4, 6, 8, 10)]); win = test_names([m for m in ms if m in (4, 6)])
    plan = [("ours", allt, range(5)), ("noconn", mid, range(3)), ("shorteval", mid, range(3)),
            ("ours_w0.25", win, range(3)), ("ours_w1", win, range(3)), ("ours_w2", win, range(3))]
    if ms == SMALL: plan.append(("conn_only", allt, [0]))
    return job_lines("ours", plan)


def main():
    """写 q2–q5。"""
    d = os.path.join(ROOT, "results", "T20", "queue"); os.makedirs(d, exist_ok=True)
    S, L = test_names(SMALL), test_names(LARGE)
    q = {
        "q2": job_lines("opp", [("pb", S, [0]), ("gphh", S, range(5)), ("esdw", S, range(5)), ("rk", S, range(5)),
                                ("d5", S, range(5)), ("rddhand", S, range(5)), ("rhmip", S, [0])]),
        "q3": job_lines("opp", [("pb", L, [0]), ("gphh", L, range(5)), ("rhmip", L, [0])]) + job_lines("ours", [("conn_only", L, [0])]),
        "q4": ours_lines(SMALL) + ours_lines(LARGE),
        "q5": job_lines("opp", [(m, L, range(5)) for m in SEARCH]),
    }
    for k, lines in q.items():
        with open(os.path.join(d, k + ".tsv"), "w") as f: f.write("\n".join(lines) + "\n")
        print(k, len(lines))


if __name__ == "__main__":
    main()
