# -*- coding: utf-8 -*-
"""
IGT 核 τ 谱生成器 — 门态强判据检验图（Au-197, 3.5–84 keV）
================================================
面板 (a)：GELINA 实测 bin 平均截面（去趋势相对起伏）+
          粒子-振动第一性门态位置（竖线）。
面板 (b)：每 bin 门态数 vs 该 bin 相对截面（Spearman 相关）。
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.stats import spearmanr

import os
from 本地配置 import 应用字体, 中间目录, 图目录
应用字体(plt, font_manager)

Ebin = np.load(os.path.join(中间目录, "gel_E.npy"))
Sbin = np.load(os.path.join(中间目录, "gel_S.npy"))
N门 = np.load(os.path.join(中间目录, "gel_N.npy"))
rel = np.load(os.path.join(中间目录, "gel_rel.npy"))
门E = np.load(os.path.join(中间目录, "gel_doors.npy"))

rho, p = spearmanr(N门, rel)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8),
                               gridspec_kw={"height_ratios": [1.4, 1]})
# (a) 去趋势截面 + 门态
ax1.plot(Ebin / 1000, rel, "ko-", ms=5, label="GELINA 实测 bin（去趋势）")
for d in 门E:
    ax1.axvline(d / 1000, color="red", alpha=0.25, lw=1)
ax1.axhline(1, color="gray", ls=":")
ax1.text(门E[0] / 1000, 1.25, "红竖线：粒子-振动第一性门态",
         color="red", fontsize=9)
ax1.set_xscale("log")
ax1.set_xlabel("中子能量 (keV)"); ax1.set_ylabel("相对截面 σ/基线")
ax1.set_title("Au-197 门态强判据：第一性门态 vs 实测截面起伏（3.5–84 keV）")
ax1.legend(loc="upper right", fontsize=9); ax1.grid(alpha=0.25, which="both")

# (b) 散点
ax2.scatter(N门, rel, c="tab:blue", s=45, alpha=0.8)
ax2.axhline(1, color="gray", ls=":")
ax2.set_xlabel("bin 内第一性门态数"); ax2.set_ylabel("相对截面 σ/基线")
ax2.set_title("门态数 vs 截面：Spearman ρ=%.2f，p=%.2f（无显著正相关）"
              % (rho, p))
ax2.grid(alpha=0.25)

plt.tight_layout()
out = os.path.join(图目录, "Au197_门态强判据检验.png")
plt.savefig(out, dpi=130)
print("已保存:", out)
print("Spearman rho=%.3f p=%.3f" % (rho, p))
