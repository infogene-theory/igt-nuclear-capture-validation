# -*- coding: utf-8 -*-
"""
IGT — Pu-240 预测（在 IGT 适用范畴内）
================================================
强外推：全局光学势 + 形变核转动增强 + Au 标定的全局 Cγ，
仅给身份 (Z,N,I,Bn)，不喂实验 D0/S0/GG。
事后与 ENDF/RIPL 对比，诚实标注适用范畴与偏差。
"""

import os
import math
import numpy as np

from 拓扑 import 核素拓扑
from 系综 import IGT参数
from 光滑截面 import 光滑截面预测
from 运行跨核外推 import 解析MF3俘获
from 本地配置 import 数据目录, 中间目录

Z, N = 94, 146
t = 核素拓扑(Z, N)

# 全局强外推（Cγ 用 Au 标定值 0.00250）
Cγ全局 = 0.00250
模 = 光滑截面预测(t, IGT参数(Cgamma全局=Cγ全局), 全局强外推=True)

# 诊断阈值量
J锚 = 0.5
D模 = 1 / max(模.能级密度.rho(t.Bn, J锚, +1) * 1e-6, 1e-30)
S0光 = 模.光学.透射(0, 1000.0) / (2 * math.pi * math.sqrt(1000.0))
print("Pu-240 诊断（括号为 RIPL 实验）:")
print("  Bn=%.3f  形变类=%s β2=%.2f" % (t.Bn, t.形变类, t.beta2))
print("  模型 D0=%.1f eV（实验13.0）" % D模)
print("  光学 S0=%.2e（实验1.07e-4）" % S0光)

# fast 截面
E网格 = np.geomspace(1000.0, 5.0e6, 30)
sig = 模.截面(E网格)
xe, ye = 解析MF3俘获(os.path.join(数据目录, "endf_094240"))
m = (xe >= 800) & (xe <= 6e6) & (ye > 0)
yexp = np.interp(E网格, xe[m], ye[m])
比 = sig / yexp
print("\nfast 区预测 vs ENDF:")
for i in (0, 4, 8, 14, 20, 26, 29):
    print("  E=%8.0f IGT=%.4f ENDF=%.4f 比=%.2f"
          % (E网格[i], sig[i], yexp[i], 比[i]))
# 分段偏差因子
for 名, 掩 in [("1–10keV", (E网格 >= 1000) & (E网格 <= 1e4)),
               ("10–100keV", (E网格 > 1e4) & (E网格 <= 1e5)),
               ("0.1–1MeV", (E网格 > 1e5) & (E网格 <= 1e6)),
               ("1–5MeV", (E网格 > 1e6) & (E网格 <= 5e6))]:
    r = 比[掩]
    if len(r):
        print("  %-9s 中位比=%.2f 偏差因子=%.2f"
              % (名, np.exp(np.median(np.log(r))),
                 np.exp(np.mean(np.abs(np.log(r))))))
np.save(os.path.join(中间目录, "pu_E.npy"), E网格)
np.save(os.path.join(中间目录, "pu_sig.npy"), sig)
np.save(os.path.join(中间目录, "pu_xe.npy"), xe)
np.save(os.path.join(中间目录, "pu_ye.npy"), ye)
