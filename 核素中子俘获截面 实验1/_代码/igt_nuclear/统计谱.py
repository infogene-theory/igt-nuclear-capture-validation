# -*- coding: utf-8 -*-
"""
IGT 核 τ 谱生成器 — 统计谱模块
================================================
1. Back-Shifted Fermi Gas（BSFG）能级密度，含 Ignatyuk 壳修正；
2. Wigner-Dyson（GOE）能级序列生成（能级排斥）；
3. Porter-Thomas（χ², ν=1）中子约化宽度抽样；
4. 门态位置调制下的背景共振聚类生成。

统计规律无自由参数（GOE/PT 形式由对称性决定）。
"""

import math
import numpy as np
from 常量 import 全局核常量, 核半径_fm


# ===============================================================
# 一、能级密度（BSFG + Ignatyuk）
# ===============================================================
class 能级密度:
    def __init__(self, Z, N, 壳修正_MEV=0.0, Bn_MEV=None, 形变类=None):
        self.Z = Z
        self.N = N
        self.A = Z + N
        self.壳修正 = 壳修正_MEV
        self.形变类 = 形变类
        c = 全局核常量
        self.alpha = c.能级密度系数
        self.gamma壳 = c.壳阻尼
        # 配对能移
        self.Delta = 0.5 * (self._配对隙(N) + self._配对隙(Z))
        # 结合能（用于自旋截断的能量归一）
        self.Bn = Bn_MEV if Bn_MEV else 6.5
        # 结合能处的经验自旋截断（配对超流使有效转动惯量小于刚体）
        self.sigma结合 = 0.98 * self.A ** 0.29

    def 集体增强(self, E):
        """
        集体运动对能级密度的增强因子（全局规则，由形变类自动判定）：
        - 变形核：转动增强 R_rot = σ⊥² ≈ 0.0138 A^(5/3) T（K 自由度积分）；
        - 过渡核：小振动增强，保守取 ~1（Au 已在 α 标定中隐含）；
        - 球形核：1。
        """
        if self.形变类 == "变形核":
            T = self.温度(E)
            sig_perp2 = 0.0138 * self.A ** (5.0 / 3.0) * max(T, 0.1)
            return min(max(sig_perp2, 1.0), 200.0)
        return 1.0

    @staticmethod
    def _配对隙(数):
        if 数 % 2 == 0:
            return 12.0 / math.sqrt(max(数, 1))
        return 0.0

    def _渐近a(self):
        return self.alpha * self.A

    def a有效(self, U):
        """Ignatyuk 能量依赖能级密度参数 a(U)。"""
        if U <= 0:
            return self._渐近a()
        a_tilde = self._渐近a()
        修正 = 1.0 + (self.壳修正 / U) * (1.0 - math.exp(-self.gamma壳 * U))
        return a_tilde * 修正

    def 温度(self, E):
        U = max(E - self.Delta, 1e-3)
        return math.sqrt(U / self.a有效(U))

    def 自旋截断(self, E):
        """
        经验自旋截断：结合能处 σ=0.98 A^0.29（已含配对超流对转动惯量的压低），
        随能量按 sqrt(E/Bn) 平滑缩放。
        """
        比例 = math.sqrt(max(E, 0.3) / self.Bn)
        return max(self.sigma结合 * 比例, 0.5)

    def rho总(self, E):
        """总能级密度（所有 J,π），单位 MeV^-1。"""
        U = E - self.Delta
        if U <= 0:
            return 0.0
        a = self.a有效(U)
        sigma = self.自旋截断(E)
        分子 = math.exp(2.0 * math.sqrt(a * U))
        分母 = 12.0 * math.sqrt(2.0) * sigma * a ** 0.25 * U ** 1.25
        return (分子 / 分母) * self.集体增强(E)

    def rho(self, E, J, 宇称):
        """指定 (J,π) 的能级密度，MeV^-1。"""
        sigma = self.自旋截断(E)
        pJ = (2 * J + 1) / (2 * sigma ** 2) * \
             math.exp(-(J + 0.5) ** 2 / (2 * sigma ** 2))
        return self.rho总(E) * pJ * 0.5


# ===============================================================
# 二、GOE 能级序列（Wigner-Dyson 间距）
# ===============================================================
def wigner_dyson间距(n, rng=None):
    """
    生成 n 个服从 Wigner-Dyson（GOE）猜测分布的无量纲间距。
    P(s) = (π/2) s exp(-π s²/4)，均值归一为 1。
    用逆变换抽样。
    """
    if rng is None:
        rng = np.random.default_rng()
    u = rng.random(n)
    # CDF: 1 - exp(-π s²/4) = u  => s = sqrt(-4 ln(1-u)/π)
    s = np.sqrt(-4.0 * np.log(1.0 - u) / math.pi)
    return s


def 生成GOE能级(起始能量, 平均间距_eV, n能级, rng=None):
    """
    生成一条具有 Wigner-Dyson 间距的能级序列（单位 eV）。
    返回能级能量数组（从起始能量开始）。
    """
    s = wigner_dyson间距(n能级, rng)
    间距 = s * 平均间距_eV
    能级 = 起始能量 + np.cumsum(间距)
    return 能级


# ===============================================================
# 三、Porter-Thomas 中子约化宽度
# ===============================================================
def porter_thomas宽度(n, 均值, rng=None):
    """
    生成 n 个 Porter-Thomas（χ², ν=1）分布的约化宽度。
    若 X~N(0,1)，则 宽度 = 均值 · X²（因 <X²>=1）。
    """
    if rng is None:
        rng = np.random.default_rng()
    x = rng.standard_normal(n)
    return 均值 * x ** 2


# ===============================================================
# 四、门态聚类的背景共振生成
# ===============================================================
def 门态聚类能级(门态位置_eV, 每门态背景数, 背景平均间距_eV,
              聚类宽度_eV, rng=None):
    """
    在每个门态位置周围生成聚类的背景共振。

    物理模型：门态通过 V2 耦合展宽/碎裂，在门态能量 ±聚类宽度 内
    产生 GOE 统计的背景共振。聚类宽度即门态耗散宽度 Γ_door。

    返回：(所有背景共振能量数组, 每个共振所属门态索引数组)
    """
    if rng is None:
        rng = np.random.default_rng()
    能级列表 = []
    归属列表 = []
    for id门, 门能 in enumerate(门态位置_eV):
        # 该门态碎裂成 k 个背景共振
        k = max(每门态背景数, 1)
        # 背景共振在门态宽度内，先用高斯定位（展宽），再叠加 WD 排斥
        中心 = rng.normal(门能, 聚类宽度_eV * 0.5)
        s = wigner_dyson间距(k, rng)
        局域间距 = s * 聚类宽度_eV / max(k, 1) * 0.8
        局域能级 = 中心 + np.cumsum(局域间距) - np.mean(np.cumsum(局域间距))
        能级列表.extend(局域能级)
        归属列表.extend([id门] * k)
    能级 = np.array(能级列表)
    归属 = np.array(归属列表)
    # 全局排序
    order = np.argsort(能级)
    return 能级[order], 归属[order]


if __name__ == "__main__":
    # 验证 Au-197 在结合能处的能级密度
    ld = 能级密度(79, 118, 壳修正_MEV=1.87)
    # 复合核 ^198Au，结合能 Bn=6.512 MeV（相对复合核基态）
    E = 6.512
    # s波复合核 J=1,2，宇称+（l=0 保持靶核宇称）
    rho1 = ld.rho(E, 1, +1)
    rho2 = ld.rho(E, 2, +1)
    D1 = 1.0 / rho1 * 1e6      # MeV -> eV
    D2 = 1.0 / rho2 * 1e6
    print("复合核^198Au 在 Bn=%.3f MeV：" % E)
    print("  J=1 序列平均间距 D ≈ %.1f eV" % D1)
    print("  J=2 序列平均间距 D ≈ %.1f eV" % D2)
    混合间距 = 1.0 / (rho1 + rho2) * 1e6
    print("  J=1,2 混合观测间距 ≈ %.1f eV（实验 D0=15.5 eV）" % 混合间距)

    rng = np.random.default_rng(42)
    间距 = wigner_dyson间距(5000, rng)
    print("WD 间距均值 = %.3f（应≈1），最小=%.3f（应>0，能级排斥）"
          % (间距.mean(), 间距.min()))
    宽度 = porter_thomas宽度(5000, 1.0, rng)
    print("PT 宽度均值 = %.3f（应≈1），中位数=%.3f（应远<1）"
          % (宽度.mean(), np.median(宽度)))
