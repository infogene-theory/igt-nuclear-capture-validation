# IGT 核素中子俘获截面少参数预测验证

**理论口径**：元衍拓扑学（Intrinsic Generative Topology, IGT，V10.1）
**理论侧**：2 号（全息系统工程）
**状态**：已完成验证，作为 IGT"理论快速验证"一个方向的成果归档

---

## 1. 研究问题

一套理论能否在**不引入大量参数（自由标量 ≤ 3 个）**的前提下，以足够精度预测核素中子俘获截面 σ(n,γ) 随中子能量的完整曲线，并与实验/评价数据相符？

中子俘获本质是非厄米共振散射，复本征值 $E_n-i\Gamma_n/2$ 直接对应 Breit-Wigner 共振——这是 IGT 框架（核心算子 $L=H-i\Gamma$，弛豫时间 τ 为第一表征标尺）的天然适用场景。

## 2. 核心结论

1. **框架正确**：非厄米共振散射框架与 Hauser-Feshbach 统计模型自洽；
2. **积分量精确**：基准核 ¹⁹⁷Au 热截面、共振积分偏差 1–2%；
3. **fast 光滑区形状正确**：5 keV–5 MeV 截面跨约 270 倍的下降被定量复现，中位比 1.12；
4. **²⁴⁰Pu 事前（非拟合）预测成功**：仅给身份边界，阈值 D0 偏差 −6%，10 keV–1 MeV 光滑截面偏差因子 1.06–1.40；
5. **两条边界被两套独立实现坐实**：少参数框架**无法逐点定位共振/门态**（Spearman ρ≈0），**跨质量区少参数全局普适不成立**。

**IGT 能做好的范畴**：强形变（及重核）的统计平均（fast 光滑）截面，可信误差带 ±12–15%（10–100 keV 最准，±6–12%）；而非全域少参数普适或逐点共振预测。

## 3. 两套独立实现（互证）

本仓库包含**两套彼此独立、代码形态与核集均不同**的实现，其结论收敛，是本成果最有价值之处。

| 维度 | 实验 1 | 实验 2 |
|---|---|---|
| 代码形态 | 26 模块（约 48 KB） | 单文件 `igt_nuclear.py`（14.3 KB） |
| 复合核谱 | 解析 Wigner-Dyson/Porter-Thomas + 系综 | **直接构造 GOE 随机矩阵内生谱** |
| 中子透射 | 全局复光学势（Numerov 解相移） | 约化宽度系统学 |
| 特色物理 | 形变核转动增强、门态粒子-振动、KU γ 强度、GELINA 实测 | Dref 锚点机制、独立 MACS 锚点、Wigner-Poisson 检验 |
| 覆盖核 | Au / Ag / Pb / Mn / Eu / U / **Pu** | Au / **Pu** / In / Fe |

"跨质量区少参数全局普适不成立"这一负面结论由两套独立代码、不同核集、不同机制分别得到——实验 1 在 Ag/Mn/Pb（需 β2、闭壳 QRPA），实验 2 在 In/Fe（需逐核 Dref、壳修正）。结论稳健。

## 4. 目录结构

```
igt-nuclear-capture-validation/
├── 核素中子俘获截面 实验1/
│   ├── 报告/                # 9 份阶段报告 + 正式论文 P-IGT-NUC-002
│   ├── 图/                  # 8 张验证图（PNG）
│   ├── 数据/                # 7 核 ENDF 评价 + EXFOR 23253（GELINA）
│   ├── _代码/igt_nuclear/   # 26 个 Python 模块源码
│   ├── _中间/               # 复现用中间数据（npy）
│   ├── IGT核谱代码包.zip    # 代码包（Windows 适配版）
│   └── README_导航与复现.md
└── 核素中子俘获截面 实验2/
    ├── P-IGT-NUC-001_..._V1.1.md   # 实验 2 报告
    ├── igt_nuclear.py                # 单文件核心实现
    ├── parse_endf_mt102.py / calibrate_dref.py / compare_igt_endf.py / wigpois_test.py / anchor_compare.py
    ├── *.endf                         # 4 核 ENDF 评价数据
    ├── igt_predict_*.txt             # IGT 预测数据
    ├── endf_*.txt                     # ENDF 平滑数据
    └── 图/                            # 验证图（PNG）
```

## 5. 复现方法（Windows）

### 实验 1

```powershell
pip install numpy scipy matplotlib
cd "核素中子俘获截面 实验1\_代码\igt_nuclear"
python 运行Au基准.py        # → 图/Au197_基准验证.png
python 运行全能量验证.py    # → 图/Au197_全能量区验证.png
python 准备门态数据.py; python 运行门态检验.py   # → 图/Au197_门态强判据检验.png
python 运行跨核外推.py      # → 图/Phase2_跨核普适外推.png
python 运行Pu240预测.py; python 运行Pu240出图.py  # → 图/Pu240_预测.png
```

路径与中文字体（微软雅黑）由 `本地配置.py` 统一管理。

### 实验 2

```powershell
cd "核素中子俘获截面 实验2"
python igt_nuclear.py --nuclide 197Au --out igt_predict_197Au.txt --plot sigma_197Au.png
python parse_endf_mt102.py in115.endf endf_In.txt
python calibrate_dref.py 115In endf_In.txt 5000 200000
python compare_igt_endf.py 197Au
python wigpois_test.py 40 400
```

## 6. 正式论文

本方向的完整成果论文见：

- `核素中子俘获截面 实验1/报告/P-IGT-NUC-002_少参数预测核素中子俘获截面_两套独立实现互证_V1.0.md`

论文包含：摘要、理论框架（L=H−iΓ、三级门态、两次 Feshbach 投影）、两套实现与参数清单、完整结果（Au 基准、fast 形状、门态负面、跨核负面、Pu-240 事前预测、误差带）、讨论（含对"纯框架/万物理论"主张的独立性判据）、适用边界、4 个可证伪预言、结论、参考文献、附录（RIPL 常数表 + 数据来源 URL）。

实验 2 的独立报告为 `P-IGT-NUC-001_中子俘获截面少参数预测验证_V1.1.md`。

## 7. 数据来源

- **ENDF/B 评价**：IAEA point2012（`https://www-nds.iaea.org/point2012/100eV/ZAxxxxxx`）
- **RIPL-3**：`https://www-nds.iaea.org/RIPL-3/resonances/resonances0.dat`
- **GELINA 高分辨率 TOF**：Massimi et al., *Eur. Phys. J. A* **50**, 124 (2014)，DOI 10.1140/epja/i2014-14124-8；EXFOR 23253
- **IAEA 热中子俘获截面**：sgnucdat a5（`https://www-nds.iaea.org/sgnucdat/a5.htm`）
- **KADoNiS 标准截面**：Dillmann et al. (2014)

## 8. 许可证

代码部分采用 MIT License；文档与报告采用 CC BY 4.0。引用时请注明理论口径为元衍拓扑学（IGT）。

---

*本仓库为 IGT"理论快速验证"方向的归档成果。核物理非 IGT 主线（主线在 agent 工程），本方向用于检验 IGT 少参数框架在复杂耗散系统上的真实能力与边界。*
