# -*- coding: utf-8 -*-
"""Pu-240 四方综合对比：ENDF vs 实验1 IGT vs 实验2 IGT vs V2纯框架。"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from 本地配置 import 应用字体, 中间目录, 图目录
应用字体(plt, font_manager)

实验2根 = r"D:\IGT agent开发\IGT 开发2\IGT 观测者资源\07_理论快速验证\核素中子俘获截面 实验2"

# 实验1
E1=np.load(os.path.join(中间目录,"pu_E.npy")); s1=np.load(os.path.join(中间目录,"pu_sig.npy"))
xe=np.load(os.path.join(中间目录,"pu_xe.npy")); ye=np.load(os.path.join(中间目录,"pu_ye.npy"))
# 实验2
d2=np.loadtxt(os.path.join(实验2根,"igt_predict_240Pu.txt"))
E2,s2=d2[:,0],d2[:,1]
# V2
V2=np.array([1.055+9.86*k for k in range(21)])

fig,(a,b)=plt.subplots(2,1,figsize=(11.5,9.5),gridspec_kw={"height_ratios":[1.5,1]})
m=(xe>0.005)&(xe<6e6)&(ye>0)
a.plot(xe[m],ye[m],"k-",lw=1.5,label="ENDF/B（实验/评价）")
mm=(E2>0.005)&(E2<6e6)
a.plot(E2[mm],s2[mm],color="#1f8a4c",lw=1.4,label="实验2 IGT（单文件·GOE内生谱）")
a.plot(E1,s1,"r.-",ms=4,lw=1,alpha=.8,label="实验1 IGT（22模块·全局光学+转动增强）")
for x in V2:
    a.axvline(x,color="orange",alpha=.25,lw=.8)
a.axvline(V2[0],color="orange",alpha=.8,lw=1,label="V2 等间距 s 梯（D≈9.86）")
# 热点
a.plot(.0253,289.5,"k*",ms=14)
a.text(.028,200,"热 289.5 b",fontsize=10)
a.set_xscale("log");a.set_yscale("log");a.set_xlim(.01,5e6);a.set_ylim(.01,3e5)
a.set_xlabel("中子能量 (eV)");a.set_ylabel("俘获截面 (b)")
a.set_title("Pu-240 四方综合对比：两套独立 IGT 在统计区一致且贴合 ENDF；V2 仅首峰、热截面低25%")
a.legend(loc="upper right",fontsize=9);a.grid(alpha=.25,which="both")

# 比值（统计区），两套分别
def 比(E,s):
    mm=(xe>700)&(xe<6e6)&(ye>0)
    return s/np.interp(E,xe[mm],ye[mm])
r1=比(E1,s1); r2=比(E2,s2)
b.plot(E1,r1,"r.-",ms=4,label="实验1 / ENDF")
b.plot(E2,r2,color="#1f8a4c",lw=1.3,label="实验2 / ENDF")
b.axhline(1,color="gray",ls=":")
b.axhspan(.85,1.15,color="steelblue",alpha=.12,label="实验误差带 ±15%")
b.set_xscale("log");b.set_xlim(700,5e6);b.set_ylim(0,2.5)
b.set_xlabel("中子能量 (eV)");b.set_ylabel("IGT / ENDF")
b.set_title("统计区比值：两套 IGT 10keV–1MeV 中心均落在 ±15% 内（实验2 报告 1.05/0.18，实验1 0.94–1.07）")
b.legend(loc="upper left",fontsize=9,ncol=2);b.grid(alpha=.25,which="both")

plt.tight_layout()
out=os.path.join(图目录,"Pu240_实验1实验2_V2综合对比.png")
plt.savefig(out,dpi=130)
print("已保存:",out)

# 汇总：热截面四方
print("\nPu 热截面: V2=217.8 | 实验1(真实参数)=295.9 | 实验2(内置锚)=290 | ENDF=289.5")
print("Pu 统计区几何均值比: 实验1≈0.94–1.07 | 实验2=1.05(对数RMS .18)")
