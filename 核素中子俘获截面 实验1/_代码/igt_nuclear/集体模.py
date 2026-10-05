# -*- coding: utf-8 -*-
"""
IGT 核 τ 谱生成器 — 集体模标度模块
================================================
由 A 与形变 β2 计算巨偶极共振（GDR）、巨四极共振（GQR）等集体模参数。
全部使用全局常量标度，无逐核拟合参数。

物理依据：
- GDR 中心能量：流体力学标度 E = c1 A^-1/3 + c2 A^-1/6；
- 变形核 GDR 几何劈裂为两支；
- GDR 峰截面由 TRK 求和规则约束，∝ A^(4/3)/E_GDR。
"""

import math
from 常量 import 全局核常量, 核半径_fm


def GDR中心能量(A, beta2=0.0):
    """
    计算 GDR 中心能量（MeV）。
    球形核返回单峰能量；变形核返回两支能量 (E1, E2)。
    """
    c1 = 全局核常量.GDR系数1
    c2 = 全局核常量.GDR系数2
    E0 = c1 * A ** (-1.0 / 3.0) + c2 * A ** (-1.0 / 6.0)

    if beta2 <= 1e-6:
        return (E0,)
    # 变形核：长轴、短轴两支，几何劈裂
    E1 = E0 * (1.0 + beta2 / 3.0)
    E2 = E0 * (1.0 - beta2 / 3.0)
    return (E1, E2)


def GDR宽度(A, beta2=0.0):
    """
    GDR 宽度（MeV）。
    球形核取全局典型宽度；变形核由几何劈裂展宽。
    """
    g0 = 全局核常量.GDR典型宽度_MEV
    if beta2 <= 1e-6:
        return g0
    # 劈裂两支的间距 + 本征宽度，近似合成
    劈裂间距 = (2.0 * beta2 / 3.0) * GDR中心能量(A, 0.0)[0]
    return math.sqrt(g0 ** 2 + 劈裂间距 ** 2)


def GDR峰截面_mb(A, beta2=0.0):
    """
    GDR 峰截面（mb），由 TRK 求和规则约束。
    积分 TRK：∫σ dE ≈ 60 NZ/A  (MeV·mb)
    对洛伦兹：积分 = (π/2) σ_peak Γ
    => σ_peak = (2/π) * 60 NZ/A / Γ
    """
    Z = int(round(A * 0.4))   # 近似，实际应由外部传入
    return None  # 由带 Z,N 的函数计算


def GDR峰截面_mb_v2(Z, N, A, beta2=0.0):
    """由 TRK 求和规则计算 GDR 峰截面（mb）。"""
    g = GDR宽度(A, beta2)
    TRK积分 = 60.0 * Z * N / A          # MeV·mb
    sigma_peak = (2.0 / math.pi) * TRK积分 / g
    return sigma_peak


def GQR能量(A):
    """巨四极共振（同位旋标量）中心能量，MeV。"""
    return 全局核常量.GQR系数 * A ** (-1.0 / 3.0)


def 集体模字典(Z, N, beta2=0.0):
    """
    汇总靶核的集体模参数，供门态生成与 γ 强度函数使用。
    """
    A = Z + N
    能量 = GDR中心能量(A, beta2)
    return {
        "GDR能量": 能量,
        "GDR宽度": GDR宽度(A, beta2),
        "GDR峰截面_mb": GDR峰截面_mb_v2(Z, N, A, beta2),
        "GQR能量": GQR能量(A),
    }


if __name__ == "__main__":
    for (Z, N, b2, 名) in [(79, 118, 0.10, "Au-197"),
                           (82, 126, 0.0, "Pb-208"),
                           (63, 88, 0.28, "Eu-151")]:
        d = 集体模字典(Z, N, b2)
        Es = ["%.2f" % e for e in d["GDR能量"]]
        print("%s: GDR E=[%s] MeV, Γ=%.2f MeV, σpeak=%.0f mb, GQR=%.2f MeV"
              % (名, ", ".join(Es), d["GDR宽度"], d["GDR峰截面_mb"], d["GQR能量"]))
