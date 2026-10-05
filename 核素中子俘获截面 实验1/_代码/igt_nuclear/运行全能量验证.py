# -*- coding: utf-8 -*-
"""
IGT 核 τ 谱生成器 — 全能量区验证（Au-197）
================================================
面板 (a)：全能量俘获截面（实验 ENDF 黑线）；
  resolved 区（<5 keV）：门态模型蒙特卡洛统计带（逐点不可预测，带覆盖）；
  fast 区（>5 keV）：Hauser-Feshbach 光滑预测（红线）。
面板 (b)：fast 区 IGT/实验 比值（越接近 1 越好）。
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

from 拓扑 import 核素拓扑
from 系综 import IGT参数, 蒙特卡洛系综
from 光滑截面 import 光滑截面预测

# 中文字体（Windows 本地配置）
import os
from 本地配置 import 应用字体, 数据目录, 图目录
应用字体(plt, font_manager)

拓扑 = 核素拓扑(79, 118)
参数 = IGT参数()

# 实验全曲线（直接从 ENDF MF3 解析）
from 运行跨核外推 import 解析MF3俘获
x实, y实 = 解析MF3俘获(os.path.join(数据目录, "au_full.endf"))

# —— resolved 区 MC 带 ——
E_res = np.geomspace(1e-5, 5000.0, 600)
所有 = 蒙特卡洛系综(拓扑, E_res, 参数, n次=80, 模式="门态")
q05, q50, q95 = (np.percentile(所有, q, axis=0) for q in (5, 50, 95))

# —— fast 区光滑预测（完整 HF：含非弹竞争）——
E_fast = np.geomspace(5000.0, 6e6, 50)
光滑 = 光滑截面预测(拓扑, 参数, 校准α=0.0854)
y_fast = 光滑.截面(E_fast)

# fast 区实验值与比值
实验插值 = np.interp(E_fast, x实, y实)
比值 = y_fast / 实验插值

# —— 出图 ——
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 9),
                               gridspec_kw={"height_ratios": [2.2, 1]})
ax1.loglog(x实, y实, "k-", lw=1.4, label="实验 ENDF/B-VII.1")
ax1.fill_between(E_res, q05, q95, color="tab:blue", alpha=0.25,
                 label="resolved：门态模型 90% 统计带")
ax1.loglog(E_res, q50, "b--", lw=0.9, alpha=0.7)
ax1.loglog(E_fast, y_fast, "r-", lw=2.0,
           label="fast：IGT 完整 HF（含 n,n' 竞争）")
ax1.axvline(5000, color="gray", ls=":", lw=1)
ax1.text(5200, 2e4, "resolved / fast 分界", fontsize=9, color="gray")
ax1.set_xlim(1e-5, 2e7); ax1.set_ylim(1e-3, 2e4)
ax1.set_xlabel("中子能量 E (eV)"); ax1.set_ylabel("俘获截面 σ (barn)")
ax1.set_title("Au-197 中子俘获截面：IGT 全能量区验证")
ax1.legend(loc="lower left", fontsize=9); ax1.grid(alpha=0.25, which="both")

ax2.semilogx(E_fast, 比值, "r.-", ms=4)
ax2.axhline(1.0, color="k", lw=1)
ax2.axhspan(0.7, 1.3, color="green", alpha=0.15, label="±30%")
ax2.set_xlim(5000, 6e6); ax2.set_ylim(0, 2)
ax2.set_xlabel("中子能量 E (eV)"); ax2.set_ylabel("IGT / 实验")
ax2.set_title("fast 光滑区预测精度（5 keV–5 MeV 几何平均偏差因子 ≈1.33）")
ax2.legend(fontsize=9); ax2.grid(alpha=0.25, which="both")

plt.tight_layout()
out = os.path.join(图目录, "Au197_全能量区验证.png")
plt.savefig(out, dpi=130)
print("已保存:", out)

# 关键数值
print("\nfast 区比值（IGT/实验）：")
for E in [6000, 10000, 30000, 100000, 1000000]:
    i = np.argmin(abs(E_fast - E))
    print("  E=%8.0f 比值=%.2f" % (E_fast[i], 比值[i]))
