# -*- coding: utf-8 -*-
"""
IGT 核 τ 谱生成器 — γ 强度函数与辐射透射模块
================================================
【术语标注·审计整改】
  "普适耗散函数 F_IGT"为本模块工作用词（非 IGT 既有术语，待入术语表）；
  其形式——Kopecky-Uhl E1 强度函数 + GDR 洛伦兹 + 末态能级密度卷积——
  系核物理标准（RIPL-3）；本模块在 IGT τ 谱视角下重构并命名。

实现辐射透射系数 Tγ：
  用 Kopecky-Uhl 广义洛伦兹 E1 γ 强度函数（基于 GDR），
  对末态能级密度积分，给出辐射透射系数 Tγ。

关键：Tγ 由 γ 强度函数积分得到，独立于复合核间距 D，
      物理上 0 ≤ Tγ ≤ 1（不会因 D 减小而发散）。

  Tγ(E*) = 2π Cγ ∫_0^{E*} Eγ³ f_E1(Eγ,T) ρ总(E*−Eγ) dEγ

γ 强度函数与 GDR 光核截面关系：
  f_E1 [MeV^-3] = σ_GDR [fm²] · μc/(ℏc)² / (3π² Eγ)
"""

import math
import numpy as np

from 集体模 import GDR中心能量, GDR宽度, GDR峰截面_mb_v2
from 统计谱 import 能级密度


class γ强度函数:
    def __init__(self, 拓扑, 参数):
        self.拓扑 = 拓扑
        self.参数 = 参数
        A = 拓扑.A
        self.GDR_E = GDR中心能量(A, 拓扑.beta2)
        self.GDR_Γ = GDR宽度(A, 拓扑.beta2)
        # GDR 峰截面（mb，TRK），fm² = mb/10
        self.GDR_峰fm2 = GDR峰截面_mb_v2(
            拓扑.Z, 拓扑.N, A, 拓扑.beta2) / 10.0
        self.能级密度 = 能级密度(
            拓扑.Z, 拓扑.N,
            壳修正_MEV=拓扑.壳修正, Bn_MEV=拓扑.Bn,
            形变类=getattr(拓扑, "形变类", None))
        # f = σ_fm² * k / (3π² Eγ)，k=μc/(ℏc)²
        self.换算 = 939.6 / 197.3 ** 2

    def _GDR截面_fm2(self, Eγ, T):
        """广义洛伦兹 GDR 光核截面（fm²），多支求和 + KU 能量/温度展宽。"""
        σ = 0.0
        n支 = len(self.GDR_E)
        for Ed in self.GDR_E:
            # Kopecky-Uhl 能量依赖宽度 + 温度项
            能量项 = (Eγ / Ed) ** 1.9
            温度项 = 0.7 * 2 * math.pi * T / self.GDR_Γ
            Γ = self.GDR_Γ * (能量项 + 温度项)
            峰 = self.GDR_峰fm2 / n支
            分母 = (Eγ ** 2 - Ed ** 2) ** 2 + Eγ ** 2 * Γ ** 2
            # 标准洛伦兹：在 E=Ed 处 = 峰
            σ += 峰 * Eγ ** 2 * Γ ** 2 / 分母
        return σ

    def f_E1(self, Eγ, T):
        """E1 γ 强度函数，MeV^-3。"""
        if Eγ <= 1e-6:
            return 0.0
        σ = self._GDR截面_fm2(Eγ, T)
        return σ * self.换算 / (3.0 * math.pi ** 2 * Eγ)

    def Tγ原始(self, J, 宇称, E激发_MeV, n点=60):
        """
        F_IGT 积分（不含 Cγ、不截断），对给定复合核 (J,π)。
        E1 跃迁：末态 J'=J-1,J,J+1（去重、≥0），宇称相反。
        """
        T核 = self.能级密度.温度(E激发_MeV)
        Js = sorted(set(j for j in (J - 1.0, J, J + 1.0) if j >= 0))
        Eg = np.linspace(0.05, E激发_MeV - 0.02, n点)
        Eg = Eg[Eg > 0]
        被积 = np.zeros_like(Eg)
        for k, e in enumerate(Eg):
            U = max(E激发_MeV - e, 1e-3)
            末态 = sum(self.能级密度.rho(U, jp, -宇称) for jp in Js)
            被积[k] = e ** 3 * self.f_E1(e, T核) * 末态
        return 2.0 * math.pi * np.trapezoid(被积, Eg)

    def Tγ(self, J, 宇称, E激发_MeV, n点=60):
        """辐射透射系数（含 Cγ，物理截断 ≤1）。"""
        T = self.Tγ原始(J, 宇称, E激发_MeV, n点) * self.参数.Cgamma
        return min(T, 1.0)
