"""T20 超预算口径 (b)（用户 11-14：主口径 (a) 用最终解；(b) 保留备查）。
(b) = 只认预算到点前得到的解。每个作业一行写 results/T20/strict_budget.jsonl：{"method", "inst", "seed", "R_b", "R_a", "src"}。
- 作业集合与 code/t20_paper_tables.py 相同：opp、ours 的作业；opp_e / ours_e 里 (e) 实际触发（e_fixed ≥ 0）的替换同名作业。
- R_a = (a) 口径（最终解；status 不是 cert / uncert 记 null），备查对照用。
- 搜索型（结果有 traj）：traj 每条是完整判定通过后记下的 [资源数, 代, 评价数, 秒, 状态]，秒从 rost5 启动算。
  R_b = 秒 ≤ rost5 限时的最后一条可行（资源数 ≥ 0）的资源数；没有则 null。
  rost5 限时：对手 = 预算 TL；本方法一类 = TL − 预处理用时（作业记录的 search_tl）。
  (a) 口径无解的作业，R_b 也是 null（只有 shorteval 有 8 个：有限时域有解、长期补救后不可行）。
- shorteval：traj 是只看 360 格的解，最终解还要经搜索后的长期判定 / 补救（两种口径都不计入预算）。
  若到点时的有限时域解 = 最终有限时域解（R_finite），R_b = R_a；否则无法得到"到点解补救后的资源数"，记 null 并在 src 里注明。
- 没有 traj 的方法（PB、GPHH、接入构造、RH-MIP）：方法用时 method_sec ≤ TL + 1 秒（1 秒算计时误差）时 R_b = R_a，否则 null。
  RH-MIP 11-06 定为不截断（(a) 按实际用时报），所以 (b) 下超时后才得到的解都是 null。
用法：python3 code/t20_strict_budget.py；打印各方法 R_b 与 R_a 不同的个数。"""
import os, sys, json, glob
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrr import ROOT  # noqa: E402

T20 = os.path.join(ROOT, "results", "T20")


def r_a(rec):
    """(a) 口径：最终解；不可行记 None（与 t20_paper_tables.put 相同）。"""
    return rec["R"] if rec.get("status") in ("cert", "uncert") else None


def r_b(rec):
    """(b) 口径的资源数与来源说明。"""
    a = r_a(rec)
    if "traj" not in rec:
        over = (rec.get("method_sec") or 0.0) - rec["TL"]
        return (a, "无traj，预算内") if over <= 1.0 else (None, f"无traj，超出 {over:.0f} 秒")
    tl = rec.get("search_tl") or rec["TL"]
    ok = [t for t in rec["traj"] if t[0] >= 0 and t[3] <= tl]
    at = ok[-1][0] if ok else None
    if a is None: return None, "最终无解"
    if rec["method"] == "shorteval":
        return (a, "shorteval，到点解=最终解") if at == rec["R_finite"] else (None, "shorteval，到点解≠最终解")
    return at, ("traj，到点解=最终解" if at == a else "traj，到点前无解" if at is None else "traj，到点解更差")


def main():
    """收集作业（与论文表格同一集合），算 R_b，写 jsonl 并打印汇总。"""
    jobs = {}
    for d in ("opp", "ours"):
        for f in glob.glob(os.path.join(T20, d, "jobs", "*.json")):
            meth, inst, key = os.path.basename(f)[:-5].split("__"); jobs[(meth, inst, key)] = json.load(open(f))
    for d in ("opp_e", "ours_e"):
        for f in glob.glob(os.path.join(T20, d, "jobs", "*.json")):
            meth, inst, key = os.path.basename(f)[:-5].split("__"); rec = json.load(open(f))
            if rec.get("e_fixed", -1) >= 0: jobs[(meth, inst, key)] = rec      # (e) 实际触发才用新结果（11-09 审核）
    cnt = Counter(); out = os.path.join(T20, "strict_budget.jsonl"); tmp = out + ".tmp"
    with open(tmp, "w") as fo:
        for (meth, inst, key), rec in sorted(jobs.items()):
            b, src = r_b(rec); a = r_a(rec)
            fo.write(json.dumps(dict(method=meth, inst=inst, seed=int(key[1:]), R_b=b, R_a=a, src=src), ensure_ascii=False) + "\n")
            cnt[(meth, "相同" if a == b else src)] += 1
    os.replace(tmp, out)
    print("写入", out, "共", len(jobs), "行")
    for k in sorted(cnt): print(k, cnt[k])


if __name__ == "__main__":
    main()
