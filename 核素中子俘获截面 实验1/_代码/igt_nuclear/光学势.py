# -*- coding: utf-8 -*-
"""
IGT 核 τ 谱生成器 — 全局中子光学势散射模块
================================================
强外推（跨核普适性）专用：不逐核输入实验强度函数 S0/S1，而是用一套
【全局固定】的复 Woods–Saxon 光学势（Becchetti–Greenlees 1969 型，
A>40、低能中子），数值求解各 l 波散射相移，直接给透射系数 T_l。

全局势参数不随核调整（仅随 Z,A 有显式系统学公式），故外推核零自由参数。

复径向方程（l 波，含离心位垒；低能忽略自旋轨道对强度函数的小修正）：
  u'' = [ l(l+1)/r² + (2μ/ℏ²)(V(r) + iW(r) − E) ] u
用复系数 Numerov 积分，远场匹配球 Hankel 得 S 矩阵，T_l = 1−|S|²。
"""

import math
import cmath
import numpy as np
from scipy.special import spherical_jn, spherical_yn

# 约化质量系数（中子-核，单位 MeV⁻¹·fm⁻²）：2μ/ℏ² ≈ 0.04776·A/(A+1)
def 质量系数(A):
    return 0.0477606 * A / (A + 1.0)


def ws(r, R, a, 深):
    """Woods–Saxon 型（深度可为复数）。"""
    return 深 / (1.0 + math.exp((r - R) / a))


class 全局光学势:
    """Becchetti–Greenlees 全局中子光学势（低能 E≈0–5 MeV）。"""

    def __init__(self, Z, N, E_MEV=0.0):
        self.Z, self.N, self.A = Z, N, Z + N
        self.E = E_MEV
        不对称 = (N - Z) / self.A
        # BG 全局参数（MeV / fm）
        self.V深 = 56.3 - 0.32 * E_MEV - 24.0 * 不对称
        self.W面 = max(13.0 - 0.25 * E_MEV - 12.0 * 不对称, 0.0)
        self.W体 = max(0.22 * E_MEV - 1.56, 0.0)
        self.r0, self.a0 = 1.17, 0.75
        self.rd, self.ad = 1.26, 0.58
        self.rv = 1.17

    def 势(self, r):
        """复光学势 V + i W（约定吸收为 −iΓ，此处取 −i|W|）。"""
        Rv = self.r0 * self.A ** (1.0 / 3.0)
        Rd = self.rd * self.A ** (1.0 / 3.0)
        Rb = self.rv * self.A ** (1.0 / 3.0)
        v = -ws(r, Rv, self.a0, self.V深)          # 吸引实部（负）
        w面 = -ws(r, Rd, self.ad, self.W面)        # 面吸收
        w体 = -ws(r, Rb, self.a0, self.W体)        # 体吸收
        return complex(v, w面 + w体)

    def 透射(self, l, E_eV, R外_fm=None, h=0.5):
        """解 l 波散射，返回 T_l=1−|S|² 与（s波）等效强度函数。"""
        E = E_eV * 1e-6
        k = math.sqrt(质量系数(self.A) * E)        # fm⁻¹（自由波数）
        if R外_fm is None:
            # 匹配半径：核半径 + 约 1.2 个约化波长
            R外_fm = self.r0 * self.A ** (1.0 / 3.0) + min(1.2 / k, 60.0)
        n = int(R外_fm / h) + 1
        R外 = n * h
        B = 质量系数(self.A)

        def f(r):
            离 = 0.0 if (r <= 1e-12 and l == 0) else (
                1e10 if r <= 1e-12 else l * (l + 1.0) / r ** 2)
            return 离 + B * (self.势(r) - E)

        # Numerov 积分（复），初值 u0=0, u1∝h^(l+1)
        u = [0j, complex(h ** (l + 1), 0.0)]
        for i in range(1, n):
            r = i * h
            fn, f0, fp = f(r - h), f(r), f(r + h)
            右 = (2 * u[-1] * (1 + 5 * h ** 2 * f0 / 12)
                  - u[-2] * (1 - h ** 2 * fn / 12))
            u.append(右 / (1 - h ** 2 * fp / 12))
        uR = u[-1]
        duR = (u[-1] - u[-3]) / (2 * h)
        beta = R外 * duR / uR                       # 无量纲对数导数

        rho = k * R外
        j = spherical_jn(l, rho)
        nyn = spherical_yn(l, rho)
        jp = spherical_jn(l, rho, derivative=True)
        np_ = spherical_yn(l, rho, derivative=True)
        hm = j - 1j * nyn
        hp = j + 1j * nyn
        hmp = (jp - 1j * np_)
        hpp = (jp + 1j * np_)
        S = -(rho * hmp - beta * hm) / (rho * hpp - beta * hp)
        T = 1.0 - abs(S) ** 2
        return min(max(T, 0.0), 1.0)


if __name__ == "__main__":
    # 在 Au-197 验证：s波 T/(2π√E) 应≈S0=1.9e-4；p波检查
    op = 全局光学势(79, 118)
    for E in (100.0, 1000.0, 10000.0):
        T0 = op.透射(0, E)
        S0等效 = T0 / (2 * math.pi * math.sqrt(E))
        T1 = op.透射(1, E)
        print("E=%6.0f eV  T0=%.4f  等效S0=%.2e  T1=%.4f"
              % (E, T0, S0等效, T1))
