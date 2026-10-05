# -*- coding: utf-8 -*-
"""
Pu-240 三方对比：V2纯框架轨表  vs  ENDF真实共振/热截面  vs  我们IGT
"""
import os
import math
import numpy as np

from 本地配置 import 数据目录

路径 = os.path.join(数据目录, "endf_094240")
行 = open(路径, errors="ignore").read().split("\n")

# ---- 1. 解析 MF2/MT151 共振（ER,AJ,GT,GN,GG）----
共振 = []
for l in 行:
    try:
        mf = int(l[70:72]); mt = int(l[72:75])
    except Exception:
        continue
    if mf == 2 and mt == 151:
        f = []
        for k in range(6):
            s = l[k*11:(k+1)*11].strip()
            try: f.append(float(s))
            except Exception: f.append(np.nan)
        if not any(math.isnan(x) for x in f[:3]):
            ER, AJ, GT = f[:3]
            GN = f[3] if not math.isnan(f[3]) else 0
            GG = f[4] if len(f) > 4 and not math.isnan(f[4]) else np.nan
            if -50 < ER < 5000 and 0 < AJ < 20:
                共振.append((ER, AJ, GN, GG, GT))

正能 = sorted([r for r in 共振 if r[0] > 0], key=lambda x: x[0])
束缚 = [r for r in 共振 if r[0] <= 0]
print("束缚态(近阈):")
for (ER, AJ, GN, GG, GT) in 束缚:
    print("  ER=%9.3f J=%.1f GN=%.4f GG=%.4f" % (ER, AJ, GN, GG))
print("\n真实低能正能共振（0–200 eV）:")
低 = [r for r in 正能 if r[0] < 200]
for (ER, AJ, GN, GG, GT) in 低:
    print("  ER=%8.3f J=%.1f GN=%.4f GG=%.4f GT=%.4f" % (ER, AJ, GN, GG, GT))

# 间距统计
if len(低) > 2:
    e = np.array([r[0] for r in 低])
    d = np.diff(e)
    print("\n0–200eV 真实共振数=%d, 间距均值=%.2f, 间距范围 %.2f–%.2f"
          % (len(e), d.mean(), d.min(), d.max()))

# ---- 2. MF3/MT102 热截面 ----
def 解析MF3(路径):
    段=False; i=0; xs=[]; ys=[]
    while i < len(行):
        l=行[i]
        try: mf=int(l[70:72]); mt=int(l[72:75])
        except Exception: i+=1; continue
        if mf==3 and mt==102:
            t=行[i+1]; NR=int(t[44:55]); NP=int(t[55:66]); j=i+2
            j += int(np.ceil(2*NR/6.0)); 需=NP
            while 需>0 and j<len(行):
                c=行[j]
                for k in range(3):
                    try:
                        x=float(c[0+k*22:11+k*22]); y=float(c[11+k*22:22+k*22])
                        xs.append(x); ys.append(y); 需-=1
                    except Exception: break
                    if 需<=0: break
                j+=1
            break
        i+=1
    return np.array(xs),np.array(ys)

xe, ye = 解析MF3(路径)
热 = float(np.interp(0.0253, xe, ye))
print("\nENDF 热中子(0.0253eV)俘获截面 = %.2f b" % 热)
# RI
m = xe >= 0.5
RI = np.trapezoid(ye[m], np.log(xe[m]))
print("ENDF 共振积分(0.5eV起) = %.0f b" % RI)
