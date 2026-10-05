# -*- coding: utf-8 -*-
"""
IGT — Windows 本地路径与中文字体统一配置
================================================
所有运行脚本通过本模块获取数据/图/中间目录与字体，
避免硬编码 Linux 路径。
"""

import os

根目录 = r"D:\IGT agent开发\IGT 开发2\IGT 观测者资源\07_理论快速验证\核素中子俘获截面 实验1"
数据目录 = os.path.join(根目录, "数据")
图目录 = os.path.join(根目录, "图")
中间目录 = os.path.join(根目录, "_中间")
报告目录 = os.path.join(根目录, "报告")

字体路径 = r"C:\Windows\Fonts\msyh.ttc"
字体名 = "Microsoft YaHei"

for _d in (数据目录, 图目录, 中间目录, 报告目录):
    os.makedirs(_d, exist_ok=True)


def 应用字体(plt, font_manager):
    """配置 matplotlib 中文字体（微软雅黑）。"""
    font_manager.fontManager.addfont(字体路径)
    plt.rcParams["font.family"] = [字体名, "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
