# -*- coding: utf-8 -*-
"""
IGT 核 τ 谱生成器 — 宽度赋值模块
================================================
为每个共振赋值中子部分宽度 Γn 与辐射宽度 Γγ。

中子宽度：
- 约化宽度 Γn° 服从 Porter-Thomas（χ²,ν=1）；
- 均值 <Γn°> = S0 × D0（S0 强度函数，D0 平均间距）；
- 物理宽度（s波穿透因子 ∝√E）：Γn(E) = Γn° × √(E[eV])。

辐射宽度：
- 退激通道极多（ν≫1），中心极限 → Γγ 近似常数 Gg，仅小涨落。
"""

import numpy as np


def 分配自旋(共振位置, 靶核自旋, rng):
    """
    s波复合核自旋 J = I±1/2，按统计权重 g_J 随机分配。
    返回与共振等长的 J 数组。
    """
    J候选 = sorted(set([靶核自旋 + 0.5, abs(靶核自旋 - 0.5)]))
    # 统计权重 ∝ (2J+1)
    权重 = np.array([2 * J + 1 for J in J候选], dtype=float)
    权重 /= 权重.sum()
    return rng.choice(J候选, size=len(共振位置), p=权重)


def 赋值宽度(共振位置_eV, S0, D0_eV, Gg_meV, 靶核自旋,
           rng, gamma涨落=0.1):
    """
    为共振赋值宽度。

    返回字典：
      J：自旋数组
      Gn0：约化中子宽度（meV，能量无关）
      Gn：物理中子宽度函数，调用  Gn(E, i)
      Gg：辐射宽度（meV）
    """
    n = len(共振位置_eV)
    J = 分配自旋(共振位置_eV, 靶核自旋, rng)

    # 平均约化中子宽度（meV）= S0 × D0
    均值Gn0_meV = S0 * D0_eV * 1000.0
    # Porter-Thomas 抽样
    x = rng.standard_normal(n)
    Gn0_meV = 均值Gn0_meV * x ** 2

    # 辐射宽度：常数 + 小涨落（多通道平均，ν大）
    Gg_数组 = Gg_meV * np.clip(rng.normal(1.0, gamma涨落, n), 0.5, 1.5)

    def 物理Gn(E_eV, i):
        return Gn0_meV[i] * np.sqrt(max(E_eV, 0.0)) / 1000.0  # 返回 eV

    return {
        "J": J,
        "Gn0_meV": Gn0_meV,
        "Gg_meV": Gg_数组,
        "物理Gn": 物理Gn,
        "均值Gn0_meV": 均值Gn0_meV,
    }


if __name__ == "__main__":
    rng = np.random.default_rng(1)
    # Au 参数
    S0, D0, Gg, I = 1.90e-4, 15.5, 128, 1.5
    均值 = S0 * D0 * 1000
    print("Au 平均约化中子宽度 <Γn°> = %.3f meV" % 均值)
    # 生成几个共振位置演示
    pos = np.array([4.9, 45, 100, 200, 500.0])
    w = 赋值宽度(pos, S0, D0, Gg, I, rng)
    print("共振 J：", w["J"])
    print("Γn°(meV)：", np.round(w["Gn0_meV"], 3))
    print("Γγ(meV)：", np.round(w["Gg_meV"], 1))
    for i, E in enumerate(pos):
        print("  E=%.1f eV 处 Γn=%.4f eV" % (E, w["物理Gn"](E, i)))
