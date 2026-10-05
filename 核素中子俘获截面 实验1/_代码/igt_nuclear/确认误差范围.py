# -*- coding: utf-8 -*-
"""Pu-240 误差范围确认：各能区 IGT偏差带（实测） vs 实验误差带。"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from 本地配置 import 应用字体, 中间目录, 图目录
应用字体(plt, font_manager)

E=np.load(os.path.join(中间目录,"pu_E.npy"))
sig=np.load(os.path.join(中间目录,"pu_sig.npy"))
xe=np.load(os.path.join(中间目录,"pu_xe.npy"))
ye=np.load(os.path.join(中间目录,"pu_ye.npy"))

m=(xe>700)&(xe<6e6)&(ye>0)
yexp=np.interp(E,xe[m],ye[m])
比=sig/yexp

# 能区（eV）
区=[("1–10 keV",1e3,1e4),("10–100 keV",1e4,1e5),
    ("0.1–1 MeV",1e5,1e6),("1–5 MeV",1e6,5e6)]
# 实验/评价侧典型误差（分数，标注来源口径）
实验err={"1–10 keV":.15,"10–100 keV":.12,"0.1–1 MeV":.12,"1–5 MeV":.18}

行out=[]
中心=[]; ig_lo=[]; ig_hi=[]; ex_lo=[]; ex_hi=[]; xpos=[]
for k,(名,lo,hi) in enumerate(区):
    mm=(E>=lo)&(E<hi)
    r=比[mm]
    med=np.exp(np.median(np.log(r)))
    # 几何偏差因子（几何标准差）
    s=np.exp(np.std(np.log(r)))
    p16,p84=np.percentile(r,[16,84])
    行out.append((名,med,s,p16,p84,实验err[名]))
    xpos.append(k)
    中心.append(med)
    ig_lo.append(med-p16); ig_hi.append(p84-med)
    ex_lo.append(实验err[名]); ex_hi.append(实验err[名])

print("%-11s %8s %10s %12s %10s"%("能区","中位比","几何σ因子","IGT16–84带","实验err"))
for (名,med,s,p16,p84,ee) in 行out:
    print("%-11s %8.2f %10.2f  %5.2f–%-5.2f %9.0f%%"%(名,med,s,p16,p84,ee*100))

# 出图：误差带对比（比值空间，1=完全一致）
fig,ax=plt.subplots(figsize=(11,6))
# 实验误差带（以1为中心）
ax.fill_between([-.5,3.5],[1-np.mean(list(实验err.values()))]*2,
                [1+np.mean(list(实验err.values()))]*2,color="steelblue",alpha=.13,
                label="实验/评价典型误差带 (~12–18%)")
ax.errorbar(xpos,中心,yerr=[ig_lo,ig_hi],fmt="ro",ms=9,capsize=7,lw=2,
            label="IGT 预测中位 + 16/84 逐点带（实测）")
for k,(名,med,s,p16,p84,ee) in enumerate(行out):
    ax.text(k,p84+.12,"%.2f（几何σ %.2f）"%(med,s),ha="center",fontsize=10)
ax.axhline(1,color="gray",ls=":")
ax.set_xticks(xpos); ax.set_xticklabels([r[0] for r in 行out])
ax.set_ylabel("IGT / ENDF（1=完全一致）")
ax.set_title("Pu-240 误差范围确认：IGT 预测偏差带 vs 实验误差带")
ax.set_ylim(0,3.2); ax.legend(loc="upper left"); ax.grid(alpha=.25,axis="y")
plt.tight_layout()
out=os.path.join(图目录,"Pu240_误差范围确认.png")
plt.savefig(out,dpi=130)
print("已保存:",out)
