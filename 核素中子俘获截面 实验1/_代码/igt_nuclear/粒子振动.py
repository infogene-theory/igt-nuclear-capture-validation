# -*- coding: utf-8 -*-
"""
IGT 核 τ 谱生成器 — 粒子-振动门态第一性生成器
================================================
强判据（隐藏拓扑序）核心：复合核门态 = 中子单粒子轨道 ⊗ 靶核声子。

能量零点：靶核基态 + 自由中子阈值（Ex=Bn）。
  门态中子能量当量：E_door = ε_j（弱束缚轨道，相对连续）+ Σ nλ Eλ
  门态 Jπ：j_轨道 ⊗ (⊗ 声子角动量)，宇称相乘。

- 0 声子：弱束缚轨道（束缚门态，E<0）；
- 1、2 声子：低能集体声子组合抵消弱束缚能，给出近阈值（keV 精度）门态。

门态位置精度 ~keV（单粒子 ε 的 WS 误差），故检验在 keV 尺度：
门态附近实验共振是否更密、fast 光滑曲线 keV 起伏（如 Au 9.3 keV 隆起）是否对应。
"""

import itertools
import numpy as np

from 伍德萨克森 import WoodsSaxon求解器


def 两体耦合(J1, J2):
    """两角动量耦合允许的总 J 集合。"""
    lo = abs(J1 - J2)
    hi = J1 + J2
    n = int(round(hi - lo)) + 1
    return [round(lo + k, 6) for k in range(n)]


class 粒子振动生成器:
    def __init__(self, 拓扑, V0校准=-45.2, 声子最大数=2,
                 轨道窗口=(-3.0, 0.0)):
        self.拓扑 = 拓扑
        ws = WoodsSaxon求解器(拓扑.A)
        ws.V0 = V0校准
        # 近阈值弱束缚中子轨道：(ε, l, j)
        轨道 = ws.费米能级附近轨道(l最大=6, 窗口=轨道窗口)
        self.轨道 = [(e, l, j) for (e, l, j, 符号) in 轨道]
        # 声子库：靶核离散集体态 (E, J, π)
        self.声子 = list(getattr(拓扑, "离散态", []))
        self.声子最大数 = 声子最大数

    def _声子配置(self, E目标_MEV, 容差_MEV=4.0):
        """
        枚举 0..N 声子组合（可重复，玻色），返回 (E声子, J集合, π)。
        仅保留总能量在目标窗口内（抵消弱束缚轨道能）。
        """
        结果 = [(0.0, [0.0], +1)]
        for n in range(1, self.声子最大数 + 1):
            # 可重复组合（玻色声子）
            for 组合 in itertools.combinations_with_replacement(
                    range(len(self.声子)), n):
                E = sum(self.声子[i][0] for i in 组合)
                π = int(np.prod([self.声子[i][2] for i in 组合]))
                # 角动量递归耦合
                Js = [0.0]
                for i in 组合:
                    Jp = self.声子[i][1]
                    Js = sorted(set(
                        j2 for j1 in Js for j2 in 两体耦合(j1, Jp)))
                if abs(E - E目标_MEV) <= 容差_MEV:
                    结果.append((E, Js, π))
        return 结果

    def 生成门态(self, E窗口_MEV=(-0.2, 0.05), 容差_MEV=4.0):
        """
        生成落在中子能量 E窗口（MeV）内的粒子-振动门态。
        返回 [(E_door MeV, J门, π门), ...]。
        """
        门态 = []
        for (ε, l, j轨) in self.轨道:
            π轨 = (-1) ** l
            # 需要声子能量 ≈ E窗口中心 − ε
            E目标 = 0.5 * (E窗口_MEV[0] + E窗口_MEV[1]) - ε
            for (E声, Js声, π声) in self._声子配置(
                    E目标, 容差_MEV):
                E门 = ε + E声
                if E窗口_MEV[0] <= E门 <= E窗口_MEV[1]:
                    π门 = π轨 * π声
                    for J声 in Js声:
                        for J门 in 两体耦合(j轨, J声):
                            门态.append((round(E门, 6), J门, π门))
        # 去重、排序
        门态 = sorted(set(门态))
        return 门态


if __name__ == "__main__":
    from 拓扑 import 核素拓扑
    拓扑 = 核素拓扑(79, 118)
    gen = 粒子振动生成器(拓扑)
    # 近阈值 0–50 keV 门态（对应 Au 9.3 keV 隆起）
    门 = gen.生成门态(E窗口_MEV=(-0.05, 0.05), 容差_MEV=3.5)
    print("近阈值（−50–50 keV）粒子-振动门态：")
    for (E, J, π) in 门:
        print("  E=%8.3f keV  J=%.1f%s" %
              (E * 1000, J, "+" if π > 0 else "-"))
