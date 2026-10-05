# -*- coding: utf-8 -*-
"""Pu-240 预测出图：fast 截面 IGT vs ENDF + 分段比值。"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

import os
from 本地配置 import 应用字体, 中间目录, 图目录
应用字体(plt, font_manager)

E = np.load(os.path.join(中间目录, "pu_E.npy"))
sig = np.load(os.path.join(中间目录, "pu_sig.npy"))
xe = np.load(os.path.join(中间目录, "pu_xe.npy"))
ye = np.load(os.path.join(中间目录, "pu_ye.npy"))

m = (xe > 700) & (xe < 6e6) & (ye > 0)
yexp = np.interp(E, xe[m], ye[m])
比 = sig / yexp

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8.5),
                               gridspec_kw={"height_ratios": [1.5, 1]})
ax1.plot(xe[m] / 1e6, ye[m], "k-", lw=1.4, label="ENDF/B 实验")
ax1.plot(E / 1e6, sig, "r.-", ms=5, label="IGT 全局强外推（零逐核参数）")
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel("中子能量 (MeV)")
ax1.set_ylabel("俘获截面 (b)")
ax1.set_title("Pu-240 中子俘获截面预测（锕系强形变核，Bn=5.242 MeV）")
ax1.legend(loc="upper right")
ax1.grid(alpha=0.25, which="both")

ax2.plot(E / 1e6, 比, "bo-", ms=4)
ax2.axhline(1, color="gray", ls=":")
ax2.set_xscale("log")
ax2.set_xlabel("中子能量 (MeV)")
ax2.set_ylabel("IGT / ENDF")
ax2.set_title("比值：10keV–1MeV 偏差因子 1.06–1.4；低能略高、高能(缺(n,2n))偏低")
ax2.set_ylim(0, 3.2)
ax2.grid(alpha=0.25, which="both")

plt.tight_layout()
out = os.path.join(图目录, "Pu240_预测.png")
plt.savefig(out, dpi=130)
print("已保存:", out)
