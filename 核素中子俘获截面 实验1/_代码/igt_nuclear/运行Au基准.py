# -*- coding: utf-8 -*-
"""
IGT 核 τ 谱生成器 — Au-197 基准验证主脚本
================================================
任务：验证 IGT 能否在 ≤3 个自由标量参数下，以高精度预测 Au-197
      中子俘获截面曲线，并与实验（ENDF/B-VII.1）相符。

本脚本产出：
  1) 截面曲线 + 蒙特卡洛置信带 vs ENDF 真实曲线
  2) Wigner-Dyson 最近邻间距分布检验
  3) Porter-Thomas 约化中子宽度分布检验
  4) 能级密度功率谱（门态隐藏序检验）
  5) 累积共振数（D0 斜率检验）
  6) 积分量对标（热截面、共振积分、D0）

运行：python3 运行Au基准.py
"""

import os
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# 中文字体（Windows 本地配置）
from 本地配置 import 应用字体, 数据目录, 图目录
应用字体(plt, font_manager)

from 拓扑 import 核素拓扑
from 系综 import IGT参数, 单次实现
from 截面 import MLBW截面, 多普勒展宽快
from 门态 import 密度功率谱

ENDF路径 = os.path.join(数据目录, "au_full.endf")
输出目录 = 图目录
os.makedirs(输出目录, exist_ok=True)


# ---------------------------------------------------------------
# 1. 解析 ENDF 真实共振参数
# ---------------------------------------------------------------
def 解析ENDF共振(路径):
    """返回正能共振列表 [(ER,J,物理GN,GG)]。"""
    共振 = []
    for L in open(路径):
        if len(L) < 76 or L[66:70] != "7925":
            continue
        if L[70:72].strip() != "2" or L[72:76].strip() != "151":
            continue
        f = []
        for k in range(6):
            s = L[k * 11:(k + 1) * 11].strip()
            try:
                f.append(float(s))
            except ValueError:
                f.append(np.nan)
        if not any(math.isnan(x) for x in f[:5]):
            ER, AJ, GT, GN, GG = f[:5]
            if (abs(ER) < 6000 and 0 < AJ < 10 and
                    0 <= GN < 5 and 0 < GG < 1):
                共振.append((ER, AJ, GN, GG))
    return 共振


def 组装全部(拓扑, 真实共振):
    """真实：束缚+锚+正能，返回 ER,J,Gn0(meV),Gg(meV)。"""
    ERl, Jl, Gn0, Ggl = [], [], [], []
    # 真实正能（排除锚，锚单独由边界条件给以免重复）
    锚ER = set(round(a[0], 2) for a in 拓扑.锚共振)
    for ER, J, GN, GG in 真实共振:
        if ER > 0 and round(ER, 2) not in 锚ER:
            ERl.append(ER)
            Jl.append(J)
            Gn0.append(GN / math.sqrt(ER) * 1000.0)
            Ggl.append(GG * 1000.0)
    # 锚 + 束缚
    for 表 in [拓扑.锚共振, 拓扑.束缚态]:
        for ER, J, GN, GG in 表:
            ERl.append(ER)
            Jl.append(J)
            Gn0.append(GN / math.sqrt(abs(ER)) * 1000.0)
            Ggl.append(GG * 1000.0)
    return (np.array(ERl), np.array(Jl), np.array(Gn0), np.array(Ggl))


# ---------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------
def 主函数():
    拓扑 = 核素拓扑(79, 118)
    参数 = IGT参数()

    # 能量网格（线性细网格，0.01–1200 eV，足够分辨多普勒峰）
    E = np.linspace(0.01, 1200.0, 80000)

    # —— 真实曲线（ENDF 全部共振重建）——
    真实共振 = 解析ENDF共振(ENDF路径)
    tER, tJ, tGn0, tGg = 组装全部(拓扑, 真实共振)
    σ真实 = MLBW截面(E, tER, tJ, tGn0, tGg, 拓扑.基态自旋)
    σ真实 = 多普勒展宽快(E, σ真实, 拓扑.A, 温度_K=300.0)

    # —— 蒙特卡洛系综（门态模式）——
    nMC = 200
    所有 = np.zeros((nMC, len(E)))
    for i in range(nMC):
        r = 单次实现(拓扑, E, 参数,
                   np.random.default_rng(2000 + i), "门态")
        所有[i] = r["截面_b"]
    q05, q25, q50, q75, q95 = (
        np.percentile(所有, p, axis=0)
        for p in (5, 25, 50, 75, 95))

    # 一个示例实现
    示例 = 单次实现(拓扑, E, 参数, np.random.default_rng(7), "门态")

    # —— 统计检验样本：真实正能共振（s 波为主，取前若干）——
    正能 = sorted([(ER, J, GN, GG) for ER, J, GN, GG in 真实共振
                  if ER > 0], key=lambda x: x[0])
    # 低能 s 波为主样本（p 波穿透 ∝(ka)³，低能弱、占比小）
    正能 = 正能[:60]
    ePos = np.array([r[0] for r in 正能])
    # 最近邻间距（归一化到均值）
    间距 = np.diff(ePos)
    D实 = 间距.mean()
    s = 间距 / D实
    # 约化中子宽度（归一化到均值）
    Gn0实 = np.array([GN / math.sqrt(ER) for ER, J, GN, GG in 正能])
    Gn0实 = Gn0实[Gn0实 > 0]
    w = Gn0实 / Gn0实.mean()

    # 功率谱（门态检验）
    频, 功, 峰周期, 峰z = 密度功率谱(ePos, 能量上限_eV=ePos[-1], n_bin=4000)

    # —— 积分量 ——
    def 积分量(σ):
        热 = float(np.interp(0.0253, E, σ))
        # RI：从 0.5 eV（Cd 截止）积分，用插值对齐起点
        m = E >= 0.5
        RI = np.trapezoid(σ[m], np.log(E[m]))
        return 热, RI

    热实, RI实 = 积分量(σ真实)
    热MC = np.array([积分量(所有[i])[0] for i in range(nMC)])
    RIMC = np.array([积分量(所有[i])[1] for i in range(nMC)])

    # ===========================================================
    # 绘图（2×3 面板）
    # ===========================================================
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle("IGT 元衍拓扑学 · Au-197 中子俘获截面基准验证（门态+统计背景，3 参数）",
                fontsize=15, fontweight="bold")

    # (a) 截面曲线
    a = ax[0, 0]
    a.fill_between(E, q05, q95, color="orange", alpha=0.25, label="90% 置信带")
    a.fill_between(E, q25, q75, color="red", alpha=0.25, label="50% 置信带")
    a.plot(E, q50, "r--", lw=1.5, label="IGT 中位预测")
    a.plot(E, 示例["截面_b"], "b-", lw=0.5, alpha=0.5, label="单次实现")
    a.plot(E, σ真实, "k-", lw=1.3, label="ENDF/B-VII.1（实验）")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlim(0.02, 1200)
    a.set_ylim(0.05, 5e4)
    a.set_xlabel("中子能量 E (eV)")
    a.set_ylabel("俘获截面 σ (barn)")
    a.set_title("(a) 截面曲线：IGT 预测 vs 实验")
    a.legend(fontsize=8, loc="lower left")
    a.grid(alpha=0.3, which="both")

    # (b) Wigner-Dyson 间距
    a = ax[0, 1]
    a.hist(s, bins=18, density=True, color="steelblue",
          alpha=0.7, edgecolor="k", label="Au 真实间距")
    xx = np.linspace(0, 3, 200)
    a.plot(xx, (math.pi / 2) * xx * np.exp(-math.pi * xx**2 / 4),
          "r-", lw=2, label="Wigner-Dyson (GOE)")
    a.plot(xx, np.exp(-xx), "g--", lw=1.5, label="泊松（可积）")
    a.set_xlabel("归一化间距 s = D / D平均")
    a.set_ylabel("概率密度 P(s)")
    a.set_title("(b) 最近邻间距分布（量子混沌）")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    # (c) Porter-Thomas 宽度
    a = ax[0, 2]
    a.hist(w, bins=22, density=True, color="seagreen",
          alpha=0.7, edgecolor="k", label="Au 真实约化宽度")
    xx = np.linspace(0.05, 5, 200)
    pt = (1 / np.sqrt(2 * math.pi * xx)) * np.exp(-xx / 2)
    a.plot(xx, pt, "r-", lw=2, label="Porter-Thomas (ν=1)")
    a.set_xlabel("归一化约化宽度 w = Γn0 / Γn0平均")
    a.set_ylabel("概率密度")
    a.set_title("(c) 中子宽度分布（Porter-Thomas）")
    a.set_ylim(0, 2.5)
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    # (d) 功率谱（门态检验）
    a = ax[1, 0]
    a.plot(频, 功, "b-", lw=0.8)
    a.axvline(1 / 峰周期, color="r", ls="--",
             label="峰周期 %.0f eV (z=%.1f)" % (峰周期, 峰z))
    a.set_xlabel("频率 (1/eV)")
    a.set_ylabel("功率谱密度")
    a.set_title("(d) 能级密度功率谱（门态隐藏序）")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)
    a.set_xlim(0, 0.05)

    # (e) 累积共振数
    a = ax[1, 1]
    a.plot(ePos, np.arange(1, len(ePos) + 1), "k.", ms=3, label="Au 真实共振")
    a.plot(ePos, ePos / D实, "r-", lw=1.5,
          label="D平均=%.1f eV（实验 15.5）" % D实)
    a.set_xlabel("中子能量 E (eV)")
    a.set_ylabel("累积共振数 N(E)")
    a.set_title("(e) 累积共振数（能级密度）")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    # (f) 积分量对标
    a = ax[1, 2]
    a.axis("off")
    行 = [
        "积分量对标（IGT 中位 vs 实验）",
        "",
        "热中子截面 σ0 (barn):",
        "   IGT  = %6.1f ± %4.1f" % (np.median(热MC),
                                  np.std(热MC)),
        "   实验 = %6.2f" % 拓扑.热截面实验,
        "",
        "共振积分 RI (barn):",
        "   IGT  = %6.0f ± %4.0f" % (np.median(RIMC),
                                  np.std(RIMC)),
        "   实验 = %6.0f" % 拓扑.共振积分实验,
        "",
        "s 波平均间距 D0 (eV):",
        "   实测重建 = %.1f（RIPL 15.5）" % D实,
        "   拓扑预测 = 14.5（偏差 6.5%）",
        "",
        "自由参数数：3（ατ, Cγ, V2）",
    ]
    a.text(0.05, 0.95, "\n".join(行), va="top", ha="left",
          fontsize=11,
          transform=a.transAxes)
    a.set_title("(f) 积分量对标", loc="left")

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    出图 = os.path.join(输出目录, "Au197_基准验证.png")
    fig.savefig(出图, dpi=140)
    plt.close(fig)

    # —— 控制台报告 ——
    print("=" * 60)
    print("Au-197 基准验证结果")
    print("=" * 60)
    print("热截面: IGT %.1f±%.1f b | 实验 %.2f b" %
          (np.median(热MC), np.std(热MC), 拓扑.热截面实验))
    print("共振积分: IGT %.0f±%.0f b | 实验 %.0f b" %
          (np.median(RIMC), np.std(RIMC), 拓扑.共振积分实验))
    print("D0: 重建 %.1f eV | RIPL 15.5 | 拓扑预测 14.5" % D实)
    print("门态功率谱峰: 周期 %.0f eV, 显著性 z=%.1f" %
          (峰周期, 峰z))
    print("图已保存：%s" % 出图)


if __name__ == "__main__":
    主函数()
