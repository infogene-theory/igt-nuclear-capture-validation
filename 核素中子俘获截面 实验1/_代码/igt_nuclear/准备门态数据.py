# -*- coding: utf-8 -*-
"""
IGT — 门态强判据数据准备（Au-197, GELINA 3.5–84 keV）
================================================
1) 解析 EXFOR 23253 的 bin（E中心、平均截面）；
2) 粒子-振动第一性门态（0–500 keV）；
3) 每 bin 门态数；去趋势相对截面（log-log 二次基线）；
4) 保存 gel_*.npy 到中间目录，供 运行门态检验.py 出图。
"""

import os
import re
import numpy as np

from 拓扑 import 核素拓扑
from 粒子振动 import 粒子振动生成器
from 本地配置 import 数据目录, 中间目录

# ---------- 1. 解析 EXFOR bin ----------
模式 = re.compile(r"^\s*(\d+)\.\s+(\d+)\.\s+([\d.]+)")
E低, E高, S = [], [], []
for 行 in open(os.path.join(数据目录, "exfor_23253.txt"),
               errors="ignore"):
    m = 模式.match(行)
    if m:
        a, b, s = float(m.group(1)), float(m.group(2)), float(m.group(3))
        if 1000 < a < 200000 and s > 0:
            E低.append(a)
            E高.append(b)
            S.append(s)
E低, E高, S = map(np.array, (E低, E高, S))
E中心 = 0.5 * (E低 + E高)
print("解析 bin 数 =%d（%.0f–%.0f eV）"
      % (len(E中心), E中心[0], E中心[-1]))

# ---------- 2. 粒子-振动门态（0–500 keV）----------
拓扑 = 核素拓扑(79, 118)
gen = 粒子振动生成器(拓扑)
门 = gen.生成门态(E窗口_MEV=(0.0, 0.5), 容差_MEV=4.0)
门E = np.array(sorted(set(E * 1e6 for (E, J, π) in 门)))
print("0–500 keV 第一性门态数 =%d" % len(门E))

# ---------- 3. 每 bin 门态数 ----------
N门 = np.array([
    int(np.sum((门E >= E低[i]) & (门E < E高[i])))
    for i in range(len(E中心))])

# ---------- 4. 去趋势（log-log 二次基线）----------
lx, ls = np.log(E中心), np.log(S)
系数 = np.polyfit(lx, ls, 2)
基线 = np.exp(np.polyval(系数, lx))
rel = S / 基线
print("门态数范围 %d–%d；相对截面范围 %.2f–%.2f"
      % (N门.min(), N门.max(), rel.min(), rel.max()))

# ---------- 保存 ----------
np.save(os.path.join(中间目录, "gel_E.npy"), E中心)
np.save(os.path.join(中间目录, "gel_S.npy"), S)
np.save(os.path.join(中间目录, "gel_N.npy"), N门)
np.save(os.path.join(中间目录, "gel_rel.npy"), rel)
np.save(os.path.join(中间目录, "gel_doors.npy"), 门E)
print("已保存 gel_*.npy 到:", 中间目录)
