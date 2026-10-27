# atlas：超长周期排班第三轮的脚本与日志（2026-10-23）

配套文档：上一级目录的 `第三轮_定理与证据.md`。依赖 Python 3.10+、numpy、scipy（≥ 1.9，用 `scipy.optimize.milp` / HiGHS）。所有脚本都在本目录下直接运行，例如 `python cert.py`。

## 稳态（允许等待）
| 脚本 | 作用 | 输出 |
|---|---|---|
| `net.py` | 相位提升的周期网络 ILP（核心模块）：`solve(g, K, P, reset=False)` | — |
| `run_net.py` | 重点序列 K = 1…12 的稳态最少台数（`wait` / `reset` 模式） | `net_wait_K12.json`、`net_reset.log` |
| `atlas_random.py` | 80 个随机序列的图谱（网络 ILP、固定模式、铺砌判定） | `atlas_random_80.json` |
| `atlas_stats.py` | 图谱统计（c 分布、按 W 分类、固定 vs 等待） | 终端输出 |
| `cm_check.py`、`multitile.py` | 零余量：分圆判据与 K 重铺砌穷举、水平半群形状 | `multitile.log` |
| `verify_L2.py` | 定理 5（L = 2 闭式）核对 | `verify_L2.log` |
| `balanced_stretch.py` | 引理 6′：多超周期剩余均衡拉伸的最小拉伸量（动态规划） | 终端输出 |
| `mix_stretch.py` | 引理 6″ 与定理 6.5：遍历全部余数类，给出 c 的构造上界 | 终端输出 |
| `modP_stretch.py` | 引理 6‴：模 P = 2L 的混合构造（L = 5 例外类） | 终端输出 |
| `l3_table.py` | 定理 6.4 的证书表（L = 3） | 终端输出 |
| `bound_vs_atlas.py` | 构造上界与图谱真实值核对（628 组，0 违反，611 取等） | 终端输出 |
| `exc_check.py`、`c2_l3.py`、`c1_search.py` | 例外类、C2、随机定向搜索 | 对应 `.log` |
| `statecycle.py`、`statecycle_batch.py` | 稳态可行性的状态图精确判定（与周期上限无关），批量严格化 c = 1 | `statecycle_batch.log` |
| `c2_counter.py` | C2 反例 [5,2,2,11,3] 的周期扫描 | `c2_counter.log` |
| `uniform_test.py`、`validate_k2.py`、`cyc_ilp.py` | 统一拉伸假设（被否定）、K ≥ 2 复核、固定模式循环覆盖 | 对应 `.log` / `.json` |

## 全新启动
| 脚本 | 作用 | 输出 |
|---|---|---|
| `fresh.py` | 有限时域时间展开 ILP（下界） | — |
| `fresh_tail.py` | 末段强制周期、可无限延拓的排班（严格上界） | `fresh_tail_wait.log`、`fresh_tail_reset.log` |
| `fresh_tail_fix.py` | 上下界不一致处加大搜索 | `fresh_tail_fix.log` |
| `cert.py` | 流体下界的精确有理对偶证书 | `cert.json` |
| `weight_proof.py` | 定理 8 的显式权重核对（k ≤ 40） | 终端输出 |
| `fresh_lp*.py`、`fresh_L2_fit.py` | 流体值、各种灵活性、L = 2 拟合 | `fresh_lp*.json/log`、`fresh_L2_fit.json` |
| `c4_check.py` | 猜想 C4（区间铺砌 ⇔ 零启动代价）核对，478/478 | 终端输出 |
| `irup_check.py` | N_fr 与 max(⌈αK⌉, N) 的比较（`wait` / `reset`） | `irup_wait.log`、`irup_reset.log` |
| `run_fresh.py` | 早期的有限时域扫描（已被上面两者取代，保留日志） | `fresh_wait.log` |

## 随机维修、分组
| 脚本 | 作用 | 输出 |
|---|---|---|
| `stoch.py` | FIFO 策略（上一轮） | `stoch.log` |
| `stoch2.py` | 相位平衡 PB 策略，K = 1…128 | `stoch2.log`、`stoch2.json` |
| `stoch3.py` | 长时域（T = 60000）、3 个种子复核 | `stoch3.log` |
| `grouping.py` | 分组损失与上界 | 终端输出 |

依赖的旧模块（从 `../../启发式验证/` 等处复制而来）：`slotmodel.py`、`kpost.py`、`roster_sim.py`、`lb.py`、`conj2.py`。
