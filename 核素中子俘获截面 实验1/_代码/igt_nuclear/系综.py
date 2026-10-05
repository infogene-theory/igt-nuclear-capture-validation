# -*- coding: utf-8 -*-
"""
IGT 核 τ 谱生成器 — 实现系综模块
================================================
把完整管线串成单次实现，并支持蒙特卡洛多次以形成统计置信带。

单次管线：
  门态骨架 → 门态聚类背景共振（或均匀 GOE 零模型）
  → 宽度赋值（PT 中子、常数 γ）→ MLBW 截面 → 多普勒展宽

两种共振位置模式：
  "门态"：背景共振在门态周围聚类（本实验工作假设：隐藏序，非 IGT 既有术语）；
  "均匀"：纯 GOE 等密度（无门态统计零模型）。
"""

import math
import numpy as np

from 门态 import 生成门态骨架
from 统计谱 import 门态聚类能级, 生成GOE能级
from 宽度 import 赋值宽度
from 截面 import MLBW截面, 多普勒展宽快


class IGT参数:
    """3 个全局标量参数（默认值对应 Au-197 标定）。"""
    def __init__(self, alpha_tau=1.0, Cgamma=1.0, V2_eV=25.0,
               每门态背景数=4, Cgamma全局=1.0):
        self.alpha_tau = alpha_tau      # 能级密度整体标度
        self.Cgamma = Cgamma           # 辐射宽度归一
        self.V2 = V2_eV               # 门态-背景耦合（聚类宽度）
        self.每门态背景数 = 每门态背景数
        self.Cgamma全局 = Cgamma全局   # 强外推 γ 全局归一（Au标定，外推共用）


def 单次实现(拓扑, 能量网格_eV, 参数, rng, 模式="门态"):
    """
    生成一次完整的截面实现。

    返回字典：截面_b、共振位置、J、Gn0、Gg、门态位置（若有）
    """
    Emax = 能量网格_eV[-1]
    # 有效平均间距（alpha_tau 标度）
    D0 = 拓扑.D0实验 / 参数.alpha_tau
    S0 = 拓扑.S0实验
    Gg = 拓扑.Gg实验 * 参数.Cgamma
    I靶 = 拓扑.基态自旋

    # 近阈值锚共振（固定，非统计）；统计共振从锚之后由 WD 自然接续
    # （不强制锚+D0，允许抽中真实的大 gap，如 Au 4.89→46）
    锚起始 = max([a[0] for a in 拓扑.锚共振], default=0.0)
    统计起点 = 锚起始

    门态位置 = None
    if 模式 == "门态":
        k = 参数.每门态背景数
        门间距 = k * D0
        门态位置 = 生成门态骨架(Emax, 门间距, rng)
        共振, 归属 = 门态聚类能级(
            门态位置, k, D0, 聚类宽度_eV=参数.V2, rng=rng)
        共振 = 共振[(共振 > 锚起始 + 0.1 * D0) & (共振 <= Emax)]
    else:
        # GOE：从锚起始，第一个间距由 WD 自然抽样（可大可小）
        共振 = 生成GOE能级(统计起点, D0, int(Emax / D0) + 5, rng)
        共振 = 共振[(共振 > 锚起始 + 0.1 * D0) & (共振 <= Emax)]

    # 正能共振宽度赋值（统计）
    w = 赋值宽度(共振, S0, D0, Gg, I靶, rng)

    def 固定态(表):
        """把 (ER, J, 物理GN, GG) 表转成数组，等效约化宽度。"""
        if not 表:
            return None
        ER = np.array([t[0] for t in 表])
        J = np.array([t[1] for t in 表])
        Gn0 = np.array([t[2] / math.sqrt(abs(t[0])) * 1000.0 for t in 表])
        Gg0 = np.array([t[3] * 1000.0 * 参数.Cgamma for t in 表])
        return ER, J, Gn0, Gg0

    锚 = 固定态(拓扑.锚共振)
    束 = 固定态(拓扑.束缚态)

    # 合并：统计正能 + 锚 + 束缚态
    全ER, 全J, 全Gn0, 全Gg = 共振, w["J"], w["Gn0_meV"], w["Gg_meV"]
    for 块 in [锚, 束]:
        if 块 is not None:
            全ER = np.concatenate([全ER, 块[0]])
            全J = np.concatenate([全J, 块[1]])
            全Gn0 = np.concatenate([全Gn0, 块[2]])
            全Gg = np.concatenate([全Gg, 块[3]])

    # MLBW 截面
    σ = MLBW截面(能量网格_eV, 全ER, 全J, 全Gn0, 全Gg, I靶)
    # 多普勒（室温）
    σ = 多普勒展宽快(能量网格_eV, σ, 拓扑.A, 温度_K=300.0)

    return {
        "截面_b": σ,
        "共振位置": 共振,
        "J": w["J"],
        "Gn0_meV": w["Gn0_meV"],
        "Gg_meV": w["Gg_meV"],
        "门态位置": 门态位置,
    }


def 蒙特卡洛系综(拓扑, 能量网格_eV, 参数, n次=200, 模式="门态",
              种子起点=1000):
    """
    生成 n次 实现，返回截面数组 (n次, nE) 与中位/分位数。
    """
    nE = len(能量网格_eV)
    所有 = np.zeros((n次, nE))
    for i in range(n次):
        rng = np.random.default_rng(种子起点 + i)
        r = 单次实现(拓扑, 能量网格_eV, 参数, rng, 模式)
        所有[i] = r["截面_b"]
    return 所有


def 置信带(所有截面, 分位数=(5, 25, 50, 75, 95)):
    """返回各分位数曲线字典。"""
    return {q: np.percentile(所有截面, q, axis=0) for q in 分位数}
