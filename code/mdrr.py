"""MDRR 公共模块：prior/ 依赖的导入、实例 JSON 的读写。

prior/排班MOEA_小样例/多型号多地/ 下的脚本互相用同目录导入（mtms → roster_sim、lb；tools_model → roster_sim、mtms；
freshlb → mtms、tools_model、export_inst → gen），而且不依赖 ../启发式验证/ 等其他目录，所以只需把这一个目录加进 sys.path。
prior/ 下的文件不改。"""
import os, sys, json

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
PRIOR = os.path.join(ROOT, "prior", "排班MOEA_小样例", "多型号多地")
if PRIOR not in sys.path:
    sys.path.insert(0, PRIOR)

from mtms import Inst, real_types, U                          # noqa: E402,F401
from tools_model import ToolInst, tool_params, robot_tool_params  # noqa: E402,F401

BENCH = os.environ.get("MDRR_BENCH", os.path.join(ROOT, "bench"))   # 环境变量只用于"重跑生成器比 md5"
INST_DIR = os.path.join(BENCH, "instances")


def inst_paths(name, split):
    """实例名 → (JSON 路径, C++ 文本路径)。"""
    d = os.path.join(INST_DIR, split)
    return os.path.join(d, name + ".json"), os.path.join(d, name + ".txt")


def load_json(path):
    """读实例 JSON，重建 Inst；并核对 JSON 里存的工具参数与 roster_sim 内置参数一致（ToolInst 从 roster_sim 取参数）。"""
    d = json.load(open(path))
    I = Inst(d["types"], d["sites"], cost=d["cost"], name=d["name"])
    assert d["tools"] == tool_params(), "工具参数与 roster_sim 不一致"
    RP = robot_tool_params()
    assert all(d["robot_tools"][n] == RP[n] for n in I.tn), "机型工具参数与 roster_sim 不一致"
    I.meta = {k: d[k] for k in ("m", "family", "idx", "seed", "split")}
    return I


def train_names(m_list=None, fam_list=None):
    """bench/split.json 里的训练实例名（按 m、族过滤）。只返回训练实例。"""
    sp = json.load(open(os.path.join(BENCH, "split.json")))
    out = []
    for e in sp["instances"]:
        if e["split"] != "train": continue
        if m_list is not None and e["m"] not in m_list: continue
        if fam_list is not None and e["family"] not in fam_list: continue
        out.append(e["name"])
    return out


def load_train(name):
    """只允许按名字加载训练实例；测试实例在第一批中不得读取。"""
    jp, _ = inst_paths(name, "train")
    I = load_json(jp)
    assert I.meta["split"] == "train"
    return I
