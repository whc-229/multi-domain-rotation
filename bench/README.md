# MDRR 实例集（MDRR-v1，2026-10-27 生成）

生成器：`code/gen_mdrr.py`（以 `prior/.../多型号多地/gen.py` 的 `make()` 为基础）。重跑 `python3 code/gen_mdrr.py` 得到字节相同的文件。

## 设定
- 机型：SIX = L1-0、L1-1、L1-2、L3-0、L3-1、L3-2（来自 `roster_sim.py` 的真实参数）。实例只保留至少一个站点允许的机型。
- 每站：允许机型随机 2–3 种；额外单程路程 δ_s ∈ {0, 1, 2, 4} h；需求 K_s。
- 族：`K1`（所有 K_s = 1）；`K12`（K_s 在 {1, 2} 中均匀抽取）。
- 站数 m ∈ {2, 3, 4, 6, 8, 10, 15, 20}。
- 共享工具池：4 种工具（`tools_model.ToolInst`，参数取自 `roster_sim`）。机型单价 c_k = 1，资源数 = 机器人数 + 工具数。
- 网格 2 h。

## 划分（生成时固定，写在 `gen_mdrr.py` 里）
- 每个 (m, 族) 16 个实例，序号 0–15。序号 0–5 为训练，6–15 为测试。
- 共 256 个：训练 96、测试 160。
- 种子 = `int(md5("MDRR-v1|m|族|序号")[:8], 16)`。
- 第一批诊断只用训练实例；测试实例已冻结，不跑、不看、不统计。

## 文件
- `instances/train/`、`instances/test/`：每个实例两份。
  - `<名>.json`：完整定义（机型参数、站点、单价、工具参数、`split`、`seed` 等），键排序，可用 `code/mdrr.py` 的 `load_json` 读回 `Inst`。
  - `<名>.txt`：C++ 文本，格式同 `prior/.../export_inst.export`；每站一行 `δ lead K`，末两行是速率下界的整数台数（定台数起点）。
- `split.json`：全部实例的 (名, m, 族, 序号, 种子, 划分)。
- 实例名：`m{站数两位}-{族}-{序号两位}`，例如 `m06-K12-03`。

## 核对（T1）
- 重跑生成器，513 个文件 md5 全部相同。
- C++ 读入（`code/cpp/inst_dump`，读入代码 `code/cpp/mdrr_inst.h` 与 rost4 相同）与 Python 由 JSON 重算的值逐项比对，256 个文件 0 不一致。
- K1 训练实例上 Python 仿真 `certsim.run` 与 `rost4 eval` 比对 288 组，0 不一致。
