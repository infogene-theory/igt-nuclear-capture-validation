# -*- coding: utf-8 -*-
"""
IGT 核 τ 谱生成器 — 光滑平均截面模块（unresolved / fast 区）
================================================
在 resolved 共振区之上（Au: >5 keV），截面经能级密度 self-averaging 后
光滑化，用完整角动量耦合的 Hauser-Feshbach 形式以少量全局标量预测。

对复合核 (J,π)：
  σγ = (π/k²) g_J · T_in Tγ / (T_in + T_inel + Tγ)
其中：
- T_in：入射中子透射（Σ_l T_l，l 满足 Jπ 耦合）；
- T_inel：(n,n') 非弹总透射 = 离散靶态逐道 + 连续能级密度积分，
  出射 l' 严格按 J 三角与宇称筛选；
- Tγ：γ 强度函数（Kopecky-Uhl/GDR）积分，Cγ 在阈值自动锚定 GG。

中子透射（强度函数归一，a 为通道半径）：
  T_0 = 2π S0 √E
  T_1 = 2π S1 √E (ka)²/3
  T_2 = 2π S2 √E (ka)^4/5
"""

import math
import numpy as np

from 统计谱 import 能级密度
from γ强度 import γ强度函数


class 光滑截面预测:
    def __init__(self, 拓扑, 参数, 通道半径_fm=7.4, 校准α=None,
                 l最大=2, 全局强外推=False):
        self.拓扑 = 拓扑
        self.参数 = 参数
        self.a通道 = 通道半径_fm
        self.l最大 = l最大
        self.全局模式 = 全局强外推
        self.I靶 = 拓扑.基态自旋

        self.能级密度 = 能级密度(
            拓扑.Z, 拓扑.N,
            壳修正_MEV=拓扑.壳修正, Bn_MEV=拓扑.Bn,
            形变类=拓扑.形变类)
        if 校准α is not None:
            self.能级密度.alpha = 校准α
        self.γ = γ强度函数(拓扑, 参数)
        self.γ.能级密度.alpha = self.能级密度.alpha

        if 全局强外推:
            # —— 强外推：全局光学势给透射，TRK 绝对 γ，零逐核参数 ——
            from 光学势 import 全局光学势
            self.光学 = 全局光学势(拓扑.Z, 拓扑.N)
            # 预计算透射表（对数能量 0.3 keV – 7 MeV），插值加速
            self._E表 = np.logspace(np.log10(300.0), np.log10(7.0e6), 70)
            self._T表 = {ll: np.array(
                [self.光学.透射(ll, Ee) for Ee in self._E表])
                for ll in range(l最大 + 1)}
            # γ 全局归一：在基准核 Au 标定一次的全局常数，外推核共用
            self.Cγ = getattr(参数, "Cgamma全局", 1.0)
        else:
            self.S0 = 拓扑.S0实验
            self.S1 = getattr(拓扑, "S1实验", 0.12e-4)
            # d 波强度函数系统学（RIPL 典型，比 S0 略大）
            self.S2 = 2.5e-4
            self.GG = 拓扑.Gg实验
            self._自动锚定Cγ()

    # —— 中子透射 ——
    def T_l(self, l, E_eV):
        """单 l 中子透射系数（0≤T≤1）。E_eV 为中子动能。"""
        if E_eV <= 0:
            return 0.0
        if self.全局模式:
            x = np.interp(math.log10(min(max(E_eV, 300.0), 7.0e6)),
                          np.log10(self._E表), self._T表[l])
            return float(x)
        ka = 2.19696e-4 * self.a通道 * math.sqrt(E_eV)
        if l == 0:
            T = 2 * math.pi * self.S0 * math.sqrt(E_eV)
        elif l == 1:
            T = 2 * math.pi * self.S1 * math.sqrt(E_eV) * ka ** 2 / 3.0
        elif l == 2:
            T = 2 * math.pi * self.S2 * math.sqrt(E_eV) * ka ** 4 / 5.0
        else:
            T = 0.0
        return min(T, 1.0)

    @staticmethod
    def _J可达(J, J态, l):
        """
        三角条件（中子为自旋 1/2 核子）：J 在 |J态−l−1/2|..J态+l+1/2，
        且半整数一致性。
        """
        l道 = l + 0.5
        Jmin = abs(J态 - l道)
        if Jmin - 1e-9 <= J <= J态 + l道 + 1e-9:
            return abs((J - J态 - l道) - round(J - J态 - l道)) < 1e-6
        return False

    def _入射透射(self, J, 宇称, E_eV):
        """入射中子到复合核 (J,π) 的总透射（靶基态 I,π=+1）。"""
        T = 0.0
        for l in range(self.l最大 + 1):
            if 宇称 == (+1) * (-1) ** l and self._J可达(J, self.I靶, l):
                T += self.T_l(l, E_eV)
        return T

    def _非弹透射(self, J, 宇称, E_eV):
        """(n,n') 非弹总透射：离散态 + 连续态，出射 l' 严格 Jπ 筛选。"""
        T总 = 0.0

        # 离散靶态
        for (Eμ, Jμ, πμ) in getattr(self.拓扑, "离散态", []):
            Er = E_eV - Eμ * 1e6
            if Er <= 0:
                continue
            for lp in range(self.l最大 + 1):
                if 宇称 == πμ * (-1) ** lp and self._J可达(J, Jμ, lp):
                    T总 += self.T_l(lp, Er)

        # 连续靶态（能级密度积分，U 从最高离散态到 E）
        起点 = (max(t[0] for t in self.拓扑.离散态)
                if self.拓扑.离散态 else 0.1)
        if E_eV > 起点 * 1e6:
            Us = np.linspace(起点, E_eV * 1e-6, 25)
            for i in range(len(Us) - 1):
                U = 0.5 * (Us[i] + Us[i + 1])
                Er = E_eV - U * 1e6
                if Er <= 0:
                    continue
                dU = Us[i + 1] - Us[i]
                for Jμ in np.arange(0, 12.5, 0.5):
                    for πμ in (+1, -1):
                        rho = self.能级密度.rho(U, Jμ, πμ)
                        if rho <= 0:
                            continue
                        for lp in range(self.l最大 + 1):
                            if 宇称 == πμ * (-1) ** lp and self._J可达(
                                    J, Jμ, lp):
                                T总 += self.T_l(lp, Er) * rho * dU
        return min(T总, 100.0)

    def _Tγ原始(self, J, 宇称, Ex):
        return self.γ.Tγ原始(J, 宇称, Ex)

    def _自动锚定Cγ(self):
        Bn = self.拓扑.Bn
        J锚 = self.I靶 + 0.5
        D_J = 1.0 / max(self.能级密度.rho(Bn, J锚, +1) * 1e-6, 1e-30)
        目标 = 2 * math.pi * (self.GG / 1000.0) / D_J
        self.Cγ = 目标 / self._Tγ原始(J锚, +1, Bn)

    def Tγ(self, J, 宇称, Ex):
        return min(self._Tγ原始(J, 宇称, Ex) * self.Cγ, 1.0)

    def 截面(self, E网格_eV):
        """返回光滑平均俘获截面（barn）。"""
        out = np.zeros_like(E网格_eV, dtype=float)
        for iE, E in enumerate(E网格_eV):
            Ex = self.拓扑.Bn + E * 1e-6
            k = 2.19696e-4 * math.sqrt(E)
            pref = math.pi / k ** 2 / 100.0
            总 = 0.0
            for J in np.arange(0, 12.5, 0.5):
                for 宇称 in (+1, -1):
                    Tin = self._入射透射(J, 宇称, E)
                    if Tin <= 0:
                        continue
                    Tinel = self._非弹透射(J, 宇称, E)
                    Tg = self.Tγ(J, 宇称, Ex)
                    g = (2 * J + 1) / (2.0 * (2 * self.I靶 + 1))
                    总 += g * Tin * Tg / (Tin + Tinel + Tg)
            out[iE] = pref * 总
        return out
