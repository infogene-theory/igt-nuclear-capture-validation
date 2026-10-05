# -*- coding: utf-8 -*-
"""
IGT 核 τ 谱生成器 — Woods-Saxon 单粒子谱求解器
================================================
数值求解中子在 Woods-Saxon 势（含自旋轨道耦合）中的束缚态能级。

方法：在径向坐标网格上以有限差分（三点）构造径向哈密顿矩阵，
      对给定 (l, j) 对角化，取负能量（束缚态）本征值。

径向方程（约化径向波函数 u）：
  [ -ℏ²/(2m) d²/dr² + ℏ² l(l+1)/(2m r²) + V_c(r)
      + V_so(r) · <l·s> ] u = E u

  <l·s> = [j(j+1)-l(l+1)-s(s+1)]/2，s=1/2
"""

import math
import numpy as np
from 常量 import 全局核常量, 中子静能_MEV


class WoodsSaxon求解器:
    def __init__(self, A, n网格=400, r最大_fm=20.0):
        self.A = A
        c = 全局核常量
        self.R = c.半径参数_fm * A ** (1.0 / 3.0)
        self.R_so = c.自旋轨道半径_fm * A ** (1.0 / 3.0)
        self.a = c.弥散_fm
        self.a_so = c.自旋轨道弥散_fm
        self.V0 = c.势深度_MEV
        # 自旋轨道强度（等效能量·fm² 量纲的组合常量）
        # 标准 WS：V_so≈6 MeV，(ℏ/mπc)²≈2.0 fm²
        self.Vso强度 = 6.0 * 2.0

        # 径向网格（避开 r=0 奇点）
        self.r = np.linspace(0.05, r最大_fm, n网格)
        self.dr = self.r[1] - self.r[0]
        self.n = n网格

    def _中心势(self, r):
        return self.V0 / (1.0 + np.exp((r - self.R) / self.a))

    def _自旋轨道势(self, r):
        """
        V_so(r) = Vso强度 · (1/r) d/dr [ 1/(1+exp((r-Rso)/a_so)) ]
        导数为负（表面内侧），整体给出正确的能级劈裂。
        """
        f = 1.0 / (1.0 + np.exp((r - self.R_so) / self.a_so))
        # 数值导数
        dfdr = np.gradient(f, r)
        return self.Vso强度 * dfdr / r

    def 求能级(self, l, j):
        """
        对给定轨道角动量 l、总角动量 j，求束缚态能级（MeV，负值）。
        返回负能量本征值升序列表。
        """
        r = self.r
        n = self.n
        dr = self.dr

        # 动能系数 ℏ²/(2m)，单位 MeV·fm²
        动能系数 = 197.327 ** 2 / (2.0 * 中子静能_MEV)   # ≈ 20.72 MeV·fm²

        # 离心项
        离心 = 动能系数 * l * (l + 1) / r ** 2
        # 中心势
        Vc = self._中心势(r)
        # 自旋轨道 <l·s>
        ls = 0.5 * (j * (j + 1) - l * (l + 1) - 0.5 * 1.5)
        Vso = self._自旋轨道势(r) * ls

        # 对角势
        对角 = 2.0 * 动能系数 / dr ** 2 + 离心 + Vc + Vso
        # 三对角：动能 -ℏ²/2m d²/dr²
        非对角 = -动能系数 / dr ** 2

        H = np.diag(对角) + np.diag(np.full(n - 1, 非对角), 1) \
                        + np.diag(np.full(n - 1, 非对角), -1)

        # 边界：u(0)=0, u(rmax)=0（有限网格自然近似）
        本征值, _ = np.linalg.eigh(H)
        束缚 = [e for e in 本征值 if e < 0]
        return sorted(束缚)

    def 费米能级附近轨道(self, l最大=6, 窗口=(-12.0, -3.0)):
        """
        枚举所有 l<=l最大 的束缚态，返回落在能量窗口内的轨道列表。
        每项：(能量MeV, l, j, 谱学符号)
        """
        谱学字母 = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g", 5: "h", 6: "i", 7: "j"}
        结果 = []
        for l in range(l最大 + 1):
            # 去重的总角动量 j 列表
            js = sorted(set([l + 0.5, max(l - 0.5, 0.5)]))
            for jj in js:
                levels = self.求能级(l, jj)
                for 节点数, e in enumerate(levels):
                    if 窗口[0] <= e <= 窗口[1]:
                        符号 = "%d%s%d/2" % (节点数 + 1, 谱学字母[l], int(2 * jj))
                        结果.append((e, l, jj, 符号))
        return sorted(结果, key=lambda x: x[0])


if __name__ == "__main__":
    ws = WoodsSaxon求解器(197)
    print("R=%.2f fm, R_so=%.2f fm" % (ws.R, ws.R_so))
    轨道 = ws.费米能级附近轨道(窗口=(-12.0, -4.0))
    print("费米能（中子结合~6.5 MeV）附近的中子单粒子轨道：")
    for e, l, j, 符号 in 轨道:
        print("  %8.3f MeV   %s" % (e, 符号))
