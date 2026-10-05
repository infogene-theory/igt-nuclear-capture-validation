# -*- coding: utf-8 -*-
"""
IGT 核 τ 谱生成器 — 截面计算模块
================================================
多能级 Breit-Wigner（MLBW）俘获截面 + 多普勒展宽。

σγ(E) = (π/k²) Σ_J g_J Σ_i Γn_i(E) Γγ_i / [(E-Ei)² + (Γi/2)²]

  g_J = (2J+1) / [2(2I+1)]
  Γi  = Γn_i + Γγ_i
  π/k² 单位 fm²，×100 转 barn（1 b = 100 fm²）
"""

import math
import numpy as np
from 常量 import 中子波数


def 统计因子(J, 靶核自旋):
    return (2 * J + 1) / (2.0 * (2 * 靶核自旋 + 1))


def MLBW截面(能量网格_eV, 共振位置_eV, J数组, Gn0_meV, Gg_meV,
            靶核自旋, 宽度截断_eV=None):
    """
    向量化计算 MLBW 俘获截面（barn）。

    参数：
      Gn0_meV：约化中子宽度（meV），物理 Γn(E)=Gn0·√E/1000 (eV)
      Gg_meV ：辐射宽度（meV）
    """
    E = np.asarray(能量网格_eV, dtype=float)
    Er = np.asarray(共振位置_eV, dtype=float)
    nE = len(E)
    σ = np.zeros(nE)

    for i in range(len(Er)):
        J = J数组[i]
        gJ = 统计因子(J, 靶核自旋)
        # 物理中子宽度随 E（s波 √E）
        Γn = Gn0_meV[i] * np.sqrt(np.maximum(E, 0)) / 1000.0   # eV
        Γγ = Gg_meV[i] / 1000.0                                # eV
        Γi = Γn + Γγ
        # π/k²：k=2.197e-4√E fm^-1；fm² -> barn（1 b = 100 fm²）需 /100
        k = 2.19696e-4 * np.sqrt(np.maximum(E, 1e-12))
        pref = math.pi / k ** 2 / 100.0                        # fm² -> barn
        # Breit-Wigner
        峰 = pref * gJ * Γn * Γγ / ((E - Er[i]) ** 2 + (Γi / 2) ** 2)
        σ += 峰
    return σ


def 多普勒展宽(能量网格_eV, 截面_b, 靶质量A, 温度_K=300.0):
    """
    自由核高斯多普勒展宽（逐共振近似：在网格上做变宽度高斯卷积）。
    Δ_D(E) = sqrt(4 E k_B T / m)
    为稳健，采用能量相关宽度的逐段卷积。
    """
    E = np.asarray(能量网格_eV)
    out = np.zeros_like(截面_b)
    kB_per_m = 8.617333e-5 / (A * 931.494e6)   # (kT/m)，单位 eV/eV
    # 每个网格点的多普勒宽度（eV）
    累积 = np.zeros_like(E)
    for j, E0 in enumerate(E):
        if E0 <= 0:
            out[j] = 截面_b[j]
            continue
        Δ = math.sqrt(4 * E0 * 8.617333e-5 * 温度_K / (A * 931.494))
        # 局部高斯窗
        窗 = np.exp(-((E - E0) ** 2) / max(Δ ** 2, 1e-12))
        窗 /= 窗.sum()
        out[j] = np.sum(截面_b * 窗)
    return out


def 多普勒展宽快(能量网格_eV, 截面_b, 靶质量A, 温度_K=300.0):
    """
    快速多普勒：在对数网格上宽度变化慢，用代表性宽度做单次 FFT 高斯卷积。
    """
    from scipy.ndimage import gaussian_filter1d
    E = np.asarray(能量网格_eV)
    # 用能量区间中部的典型多普勒宽度换算成网格 sigma
    典型E = np.exp(np.mean(np.log(np.maximum(E, 1e-6))))
    Δ = math.sqrt(4 * 典型E * 8.617333e-5 * 温度_K / (靶质量A * 931.494))
    if E[-1] > E[0]:
        dE = np.gradient(E)
        sigma网格 = np.median(Δ / dE)
    else:
        sigma网格 = 0
    return gaussian_filter1d(截面_b, sigma网格)


if __name__ == "__main__":
    # 单共振解析验证：4.9 eV 共振
    E = np.linspace(2, 8, 6001)
    Er = np.array([4.906])
    J = np.array([2.0])
    Gn0 = np.array([0.5])     # meV
    Gg = np.array([128.0])    # meV
    σ = MLBW截面(E, Er, J, Gn0, Gg, 1.5)
    ip = np.argmax(σ)
    print("共振峰位 E=%.3f eV，峰截面=%.0f b" % (E[ip], σ[ip]))
