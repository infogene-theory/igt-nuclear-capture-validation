# 核素中子俘获截面 — IGT 快速验证实验（实验 1）

**日期**：2026-10-02　**理论**：元衍拓扑学（Intrinsic Generative Topology, IGT）　**理论侧**：2 号

---

## 一、实验问题

> 若一套理论能在**不引入大量参数（≤3 个自由标量）**前提下，以足够高精度预测某核素中子俘获截面曲线、并与实验相符——IGT 能否做到？

## 二、目录结构

| 路径 | 内容 |
|---|---|
| `报告\` | 6 份阶段报告（Phase 0 / 1 / 1.5 / 1.6 / 2 / Pu-240）+ 综合对比/误差确认 |
| `报告\P-IGT-NUC-002_..._V1.0.md` | **正式论文**：两套独立实现互证与适用边界（本方向归档成果，综合实验1、实验2） |
| `图\` | 5 张验证图（PNG） |
| `数据\` | 7 核 ENDF/B 评价 + EXFOR 23253（GELINA） |
| `IGT核谱代码包.zip` | 全部 22 个 Python 模块（Windows 适配版） |
| `_代码\igt_nuclear\` | 解压源码（同 zip，可直接运行） |
| `_中间\` | 复现用中间数据（npy） |

## 三、复现方法（Windows）

1. 依赖：`python -m pip install numpy scipy matplotlib`（本机 Python 3.14，已装 numpy 2.5 / scipy 1.18 / matplotlib 3.11）；
2. 路径与中文字体由 `本地配置.py` 统一管理（微软雅黑，无需改动）；
3. 依次运行：
   - `python 运行Au基准.py` → 图/Au197_基准验证.png
   - `python 运行全能量验证.py` → 图/Au197_全能量区验证.png
   - `python 准备门态数据.py` 然后 `python 运行门态检验.py` → 图/Au197_门态强判据检验.png
   - `python 运行跨核外推.py` → 图/Phase2_跨核普适外推.png
   - `python 运行Pu240预测.py` 然后 `python 运行Pu240出图.py` → 图/Pu240_预测.png

## 四、核心结论（一句话总览）

1. **框架正确**：中子俘获本质是非厄米共振散射，复本征值 E_n−iΓ_n/2 直接对应 Breit-Wigner 共振；
2. **积分量精确**：Au-197 热截面、共振积分偏差 1–2%；
3. **fast 光滑形状正确**：5 keV–5 MeV 偏差因子 ~1.33；
4. **门态逐点预测不成立**：少参数粒子-振动无法第一性定位共振（ρ≈0，强判据不通过）；
5. **跨质量区少参数全局普适不成立**：每补一个全局物理修好一类核、却揭示新的逐核依赖（β2、壳、配对、微观 γ）；
6. **Pu-240 成功**：作为强形变核（β2≈0.28），在适用范畴内零逐核参数预测成功（D0 −6%，10keV–1MeV 偏差 1.06–1.4）。

**IGT"能做好"的范畴**：强形变核的统计平均（fast 光滑）截面；而非全域少参数普适或逐点共振预测。

## 五、关键数据修正（开发记录）

- 中子波数 **k=2.197×10⁻⁴√E_eV fm⁻¹**（曾错 10 倍，致截面小 100 倍）；
- fm²→barn **除以 100**；
- p 波 S1=0.12×10⁻⁴（非 1.75·S0）；
- 光学势 Numerov 中间系数必须 `(1+5h²f/12)`；
- numpy 2.x 用 `np.trapezoid`（非 np.trapz）。

## 六、数据来源

- ENDF/B：IAEA point2012（`https://www-nds.iaea.org/point2012/100eV/ZAxxxxxx`）
- RIPL-3：`https://www-nds.iaea.org/RIPL-3/resonances/resonances0.dat`
- GELINA：DOI 10.1140/epja/i2014-14124-8；EXFOR 23253（GitHub iaea-nds/exfor-entry-file）
