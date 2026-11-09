"""后备上限核对 G1 的重新判定（11-09）。
code/t20_check_f0cap.py 的 G1 复用 t20_check_e.e1，那时 e1 只从新程序的输出里去掉 e_* 字段；
现装程序（b92f…）已带 (e)，它的输出也有 e_fixed / e_gaps0 / e_found_gen，于是每组都在这三个字段上记成"不同"。
check_f0cap.jsonl 里每组都存了完整的逐字段差异（diff 列出所有不同的字段），所以不必重跑，按 diff 重新判定：
- 相同 = diff 为空，或 diff 只含这三个字段，且旧程序在这三个字段上都是 −1（新程序那边已去掉）；
- 另外要求新程序 e_fixed = −1、f0_capped = 0（不触发 (e)、后备没被截）。
G2、G3 不受影响（G2 自己去掉了旧程序的 e_* 字段），原样保留判定。
原始结果 check_f0cap.jsonl 不动，重新判定写 results/T20/check_f0cap_regrade.jsonl，最后一行是汇总。
用法：python3 code/t20_check_f0cap_regrade.py"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402

E3 = {"e_fixed", "e_gaps0", "e_found_gen"}


def regrade(r):
    """一行 G1 的重新判定：返回 (相同, 通过)。"""
    d = r.get("diff") or {}
    same = not d or (set(d) <= E3 and all(v[0] == -1 and v[1] is None for v in d.values()))
    ok = same and r["e"].get("e_fixed", -1) == -1 and r["e"].get("f0_capped", 0) == 0
    return int(same), int(ok)


def main():
    """读原始结果，重判 G1，写 regrade 文件并打印汇总。"""
    T = os.path.join(ROOT, "results", "T20")
    rows = [json.loads(l) for l in open(os.path.join(T, "check_f0cap.jsonl"))]
    recs = [r for r in rows if "part" in r]
    out = []
    for r in recs:
        r = dict(r)
        if r["part"] == "G1": r["same_regrade"], r["ok_regrade"] = regrade(r)
        else: r["ok_regrade"] = r["ok"]
        out.append(r)
    summ = dict(G1=sum(r["part"] == "G1" for r in out), G1_bad=sum(r["part"] == "G1" and not r["ok_regrade"] for r in out),
                G2=sum(r["part"] == "G2" for r in out), G2_bad=sum(r["part"] == "G2" and not r["ok_regrade"] for r in out),
                G2_capped=sum(r["part"] == "G2" and r.get("capped", 0) == 1 for r in out),
                G3_bad=sum(r["part"] == "G3" and not r["ok_regrade"] for r in out))
    with open(os.path.join(T, "check_f0cap_regrade.jsonl"), "w") as f:
        for r in out: f.write(json.dumps(r, ensure_ascii=False) + "\n")
        f.write(json.dumps(summ) + "\n")
    print(json.dumps(summ))


if __name__ == "__main__":
    main()
