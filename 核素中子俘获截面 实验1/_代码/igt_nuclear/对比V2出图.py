# -*- coding: utf-8 -*-
"""
Pu-240 三方对比出图：
 (a) 低能 0.1–200 eV：ENDF 真实共振重建(黑) + V2等间距s梯(红竖线)
 (b) 关键标量：热截面 / RI，V2 vs 实验
"""
import os
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

from 本地配置 import 数据目录, 图目录, 应用字体
应用字体(plt, font_manager)

# ---- 解析真实共振（严格过滤）----
行 = open(os.path.join(数据目录, "endf_094240"), errors="ignore").read().split("\n")
共 = []
for l in 行:
    try:
        mf=int(l[70:72]); mt=int(l[72:75])
        if not (mf==2 and mt==151): continue
        f=[float(l[k*11:(k+1)*11].strip()) for k in range(5)]
    except Exception:
        continue
    ER,AJ,GT,GN,GG=f
    if AJ<0.1 or AJ>10: continue
    if not (0.01 < GG < 0.1): continue
    if -30 < ER < 400 and GN>=0:
        共.append((ER,AJ,GN,GG))

ERl,Jl,Gn0,Ggl=[],[],[],[]
for (ER,J,GN,GG) in 共:
    ERl.append(ER); Jl.append(J)
    Gn0.append(GN/math.sqrt(abs(ER))*1000); Ggl.append(GG*1000)
ERl=np.array(ERl); Jl=np.array(Jl); Gn0=np.array(Gn0); Ggl=np.array(Ggl)

from 截面 import MLBW截面
E=np.geomspace(0.05,200,40000)
sig=MLBW截面(E,ERl,Jl,Gn0,Ggl,0.0)

# V2 等间距 s 梯位置
V2梯=np.array([1.055+9.86*k for k in range(21)])

fig,(a,b)=plt.subplots(2,1,figsize=(11,9),gridspec_kw={"height_ratios":[1.6,1]})
a.plot(E,sig,"k-",lw=1,label="ENDF 真实共振重建")
for x in V2梯:
    a.axvline(x,color="red",alpha=.3,lw=1)
a.axvline(V2梯[0],color="red",alpha=.8,lw=1.5,label="V2 等间距 s 梯（D≈9.86 eV）")
a.set_xscale("log"); a.set_yscale("log")
a.set_xlim(.08,200); a.set_ylim(.05,2e5)
a.set_xlabel("中子能量 (eV)"); a.set_ylabel("俘获截面 (b)")
a.set_title("Pu-240 低能区：V2 仅首峰 1.056 eV 命中，其后等间距梯与真实不规则共振不符")
a.legend(loc="upper right"); a.grid(alpha=.25,which="both")

# 标量对比
项=["热截面 σ0 (b)","共振积分 RI (b)"]
V2=[217.83, np.nan]   # V2 "232keV" 口径不明，不参与
实验=[289.5,8452]
x=0; w=.35
b.bar(x-w/2,[V2[0]],w,label="V2 纯框架轨",color="salmon")
b.bar(x+w/2,[实验[0]],w,label="ENDF/IAEA 实验",color="steelblue")
b.text(x-w/2,V2[0]+6,"%.1f（低估25%%）"%V2[0],ha="center",fontsize=10)
b.text(x+w/2,实验[0]+6,"%.1f"%实验[0],ha="center",fontsize=10)
b.set_xticks([x]); b.set_xticklabels(["热中子截面 σ0"])
b.set_ylabel("barn"); b.set_ylim(0,340)
b.set_title("热截面：V2=217.8 b vs 实验 289.5 b（RI 实验 8452 b；V2“232keV”口径不明未列）")
b.legend(loc="upper left"); b.grid(alpha=.25,axis="y")

plt.tight_layout()
out=os.path.join(图目录,"Pu240_V2三方对比.png")
plt.savefig(out,dpi=130)
print("已保存:",out)
# 控制台
真实低能=sorted(set(round(x,2) for x in ERl if x>0.5 and x<60))
print("真实 0–60eV 共振:",真实低能)
print("V2 0–60eV 梯:",[round(x,2) for x in V2梯 if x<60])
