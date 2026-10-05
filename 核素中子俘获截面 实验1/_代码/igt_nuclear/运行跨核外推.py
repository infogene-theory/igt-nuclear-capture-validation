# -*- coding: utf-8 -*-
"""
IGT 核 τ 谱生成器 — Phase 2 跨核普适性强外推
================================================
对 6 个核（Au-197 基准 + Ag-109 / Pb-208 / Mn-55 / Eu-151 / U-238）：
- 仅给 (Z,N) 身份边界（基态 Jπ、Bn，核素基本属性）；
- 中子透射：全局光学势（不喂实验 S0/S1）；
- γ 透射：TRK 求和规则绝对定标（Cγ=1，不锚实验 GG）；
- 能级密度：全局 α=0.086（不逐核校准）；
- 对比 ENDF 实验 fast 截面，给几何平均偏差因子。
外推核零自由参数。
"""

import numpy as np

from 拓扑 import 核素拓扑
from 系综 import IGT参数
from 光滑截面 import 光滑截面预测
from 本地配置 import 数据目录, 图目录
import os

核表 = [
    (79, 118, os.path.join(数据目录, "au_full.endf"), "Au-197"),
    (47, 62, os.path.join(数据目录, "endf_047109"), "Ag-109"),
    (82, 126, os.path.join(数据目录, "endf_082208"), "Pb-208"),
    (25, 30, os.path.join(数据目录, "endf_025055"), "Mn-55"),
    (63, 88, os.path.join(数据目录, "endf_063151"), "Eu-151"),
    (92, 146, os.path.join(数据目录, "endf_092238"), "U-238"),
]


def 解析MF3俘获(路径):
    """提取 ENDF MF3/MT102 的 (能量eV, 截面b)。"""
    行 = open(路径, errors="ignore").read().split("\n")
    段 = False
    i = 0
    xs, ys = [], []
    while i < len(行):
        l = 行[i]
        try:
            mf = int(l[70:72])
            mt = int(l[72:75])
        except Exception:
            i += 1
            continue
        if mf == 3 and mt == 102:
            # TAB1 行（本 HEAD 之后第一行含 NR,NP）
            t = 行[i + 1]
            NR = int(t[44:55])
            NP = int(t[55:66])
            j = i + 2
            # 跳过插值表（NR 对，每行 6 个整数）
            需整数 = 2 * NR
            行占 = int(np.ceil(需整数 / 6.0))
            j += 行占
            # 读 NP 对（每行 3 对 = 6 个浮点数）
            需 = NP
            while 需 > 0 and j < len(行):
                c = 行[j]
                for k in range(3):
                    try:
                        x = float(c[0 + k * 22:11 + k * 22])
                        y = float(c[11 + k * 22:22 + k * 22])
                        xs.append(x)
                        ys.append(y)
                        需 -= 1
                    except Exception:
                        break
                    if 需 <= 0:
                        break
                j += 1
            break
        i += 1
    return np.array(xs), np.array(ys)


def 单核(Z, N, 路径, 名, E网格):
    拓扑 = 核素拓扑(Z, N)
    参数 = IGT参数()
    模 = 光滑截面预测(拓扑, 参数, 全局强外推=True)
    sig = 模.截面(E网格)
    xe, ye = 解析MF3俘获(路径)
    return sig, xe, ye


if __name__ == "__main__":
    import time
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    from 本地配置 import 应用字体
    应用字体(plt, font_manager)

    E网格 = np.geomspace(5000.0, 5.0e6, 22)

    # 第一步：在基准核 Au 标定全局 Cγ（锚 Au GG），外推核共用
    tAu = 核素拓扑(79, 118)
    模Au = 光滑截面预测(tAu, IGT参数(Cgamma全局=1.0), 全局强外推=True)
    J锚 = tAu.基态自旋 + 0.5
    D锚 = 1 / max(模Au.能级密度.rho(tAu.Bn, J锚, +1) * 1e-6, 1e-30)
    目标 = 2 * np.pi * (tAu.Gg实验 / 1000.0) / D锚
    Cγ全局 = 目标 / 模Au.γ.Tγ原始(J锚, +1, tAu.Bn)
    print("Au 标定全局 Cγ =%.5f（D锚=%.1f eV）" % (Cγ全局, D锚))

    # 第二步：同一全局 Cγ 跑全部 6 核
    结果 = {}
    t0 = time.time()
    for (Z, N, 路径, 名) in 核表:
        t = 核素拓扑(Z, N)
        模 = 光滑截面预测(
            t, IGT参数(Cgamma全局=Cγ全局), 全局强外推=True)
        sig = 模.截面(E网格)
        xe, ye = 解析MF3俘获(路径)
        m = (xe >= E网格[0] * 0.8) & (xe <= E网格[-1] * 1.2) & (ye > 0)
        yexp = np.interp(E网格, xe[m], ye[m])
        比 = sig / yexp
        偏差因子 = np.exp(np.mean(np.abs(np.log(比))))
        中位比 = np.exp(np.median(np.log(比)))
        结果[名] = (xe, ye, sig, 偏差因子, 中位比)
        print("%-7s 中位比=%.2f 几何偏差因子=%.2f"
              % (名, 中位比, 偏差因子))
    print("总耗时 %.1f s" % (time.time() - t0))

    # 出图：2×3
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    for ax, (Z, N, 路径, 名) in zip(axes.ravel(), 核表):
        xe, ye, sig, 偏差因子, 中位比 = 结果[名]
        m = (xe > 3000) & (xe < 8e6) & (ye > 0)
        ax.plot(xe[m] / 1e6, ye[m], "k-", lw=1.3, label="ENDF 实验")
        ax.plot(E网格 / 1e6, sig, "r.-", ms=4, label="IGT 全局外推")
        ax.set_xscale("log")
        ax.set_yscale("log")
        标 = "基准" if 名 == "Au-197" else "外推"
        ax.set_title("%s（%s）中位比=%.2f，偏差因子=%.2f"
                     % (名, 标, 中位比, 偏差因子))
        ax.set_xlabel("能量 (MeV)")
        ax.set_ylabel("俘获截面 (b)")
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=8, loc="upper right")
    fig.suptitle("IGT Phase 2 跨核普适性强外推（全局光学势 + 全局Cγ，外推核零调参）",
                 fontsize=13)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(图目录, "Phase2_跨核普适外推.png")
    plt.savefig(out, dpi=130)
    print("已保存:", out)

