# IGT 核 τ 谱生成器 — Phase 0 理论形式化

**文档编号**：IGT-NUC-P0-V1.0　**日期**：2026-10-02　**理论侧**：2 号（全息系统工程）　**状态**：Phase 0 形式化基线，待评审

---

## 第一部分　公理化基础在核系统中的具体化

### 1.1 IGT 两定义两公理（核系统投影）

**定义 D-1（本征生成元）**：一个物理系统的本征生成元 $\boldsymbol{g}$ 是确定其动力学谱所需的最小拓扑不变量集合。

**定义 D-2（τ 谱）**：系统在给定观测通道下的耗散本征模式集合，由非厄米算子的复本征值表征：
$$\mathcal{T}=\{(\omega_\lambda,\tau_\lambda)\},\qquad z_\lambda=\omega_\lambda-\frac{i}{2\tau_\lambda}$$
其中 $\omega_\lambda$ 为本征角频率，$\tau_\lambda$ 为弛豫时间。

**公理 A-1（本征生成公理）**：系统的动力学算子由本征生成元通过群表示映射 $\Phi_i$ 唯一确定：
$$L=\Phi_i(\boldsymbol{g}),\qquad \Phi_i:\mathcal{G}\longmapsto\mathfrak{H}(\mathcal{H})$$

**公理 A-2（观测投影公理）**：任何可观测量都是 τ 谱在观测通道上的投影：
$$O=\mathcal{P}_{obs}(\mathcal{T})$$

### 1.2 核系统的三级门态层级（核心结构）

中子俘获的复合核过程对应一个**三级耗散级联**，时间尺度分离：

| 层级 | 物理内容 | 特征时间 | 统计性质 |
|---|---|---|---|
| 第 0 层：集体模层（Collective） | GDR、GQR 等巨共振 | τ_co ~ ℏ/Γ_GDR ~ 10⁻²² s | 相干集体 |
| 第 1 层：门态层（Doorway） | 单粒子 ⊗ 声子（1p1h⊗phonon） | τ_door ~ 10⁻¹⁹–10⁻¹⁸ s | 稀疏、部分有序 |
| 第 2 层：背景复合核层（Background） | 复杂多体复合核态 | τ_bg ~ 10⁻¹⁶ s | GOE 量子混沌 |

### 1.3 两次 Feshbach 投影与有效算子

对三级希尔伯特空间 $\mathcal{H}=\mathcal{H}_{co}\oplus\mathcal{H}_{door}\oplus\mathcal{H}_{bg}$ 依次做 Feshbach 投影：

- 第一次投影（背景 → 门态）：门态获得自能
  $$\Sigma^{(1)}(E)=V_{db}(E-H_{bg})^{-1}V_{bd}$$
- 第二次投影（门态 → 集体/入射空间）：
  $$L_{eff}=H_{door}+\Sigma^{(1)}(E)+\Sigma^{(2)}(E)$$

有效算子的复本征值给出共振位置与宽度：
$$E_n-\frac{i}{2}\Gamma_n,\qquad \Gamma_n=\Gamma_n^\uparrow+\Gamma_n^\downarrow$$

---

## 第二部分　普适耗散函数 F_IGT（任务 #2）

### 2.1 辐射宽度的 IGT 表达式

复合核共振 i 通过发射 γ 射线退激到低能态 f 的部分宽度：
$$\Gamma_{\gamma,i}=\sum_f \Gamma_\gamma(E_i\to E_f)$$

引入 γ 强度函数（单位 Eγ 的平均约化跃迁概率）：
$$f_{E1}(E_\gamma)=\frac{1}{2J_i+1}\frac{1}{E_\gamma^3}\langle |f|\hat{\mathcal M(E_\gamma)}|i\rangle|^2$$

则平均辐射宽度：
$$\boxed{\langle\Gamma_\gamma\rangle(E_i)=\frac{16\pi}{9}\frac{1}{(\hbar c)^3}\sum_f\int_0^{E_i}E_\gamma^3 f_{E1}(E_\gamma)\rho(E_f,J_f,\pi_f)\,dE_\gamma}$$

### 2.2 γ 强度函数：Kopecky-Uhl 广义洛伦兹

E1 强度函数由 GDR 的洛伦兹形状给出（Brink-Axel 假设：强度函数与初态无关）：
$$f_{E1}(E_\gamma,T)=\frac{1}{3\pi^2\hbar^2 c^2}\frac{\sigma_{GDR}(E_\gamma,T)}{E_\gamma}$$

Kopecky-Uhl 能量/温度依赖宽度：
$$\Gamma(E_\gamma,T)=\Gamma_D\left[\left(\frac{E_\gamma}{E_D}\right)^{1.9}+0.7\frac{2\pi T}{\Gamma_D}\right]$$

标准洛伦兹光核截面：
$$\sigma_{GDR}(E_\gamma)=\sigma_{峰}\frac{E_\gamma^2\Gamma^2}{(E_\gamma^2-E_D^2)^2+E_\gamma^2\Gamma^2}$$

### 2.3 普适耗散函数

$$\boxed{F_{IGT}(E)=\int_0^E E_\gamma^3\,L_{GDR}(E_\gamma,T)\,\rho(E-E_\gamma)\,dE_\gamma}$$

该函数把"集体模强度 L_GDR"与"末态能级密度 ρ"耦合，是 IGT 对辐射耗散的普适表达。

---

## 第三部分　能级密度（BSFG + 壳修正）

**Back-Shifted Fermi Gas**：
$$\rho(U,J,\pi)=\frac{1}{2}\frac{(2J+1)}{2\sqrt{2\pi}\sigma^2}\,\frac{\exp(2\sqrt{aU})}{12\sigma a^{1/4}U^{5/4}}$$

- 能级密度参数：$a=\alpha A$（α≈0.086 全局）
- Ignatyuk 壳修正：$a(U)=\tilde a\left[1+\frac{\delta W}{U}(1-e^{-\gamma U})\right]$，γ≈0.05
- 配对移：$\Delta=\tfrac12(\Delta_N+\Delta_Z)$，$\Delta_{N/Z}=12/\sqrt{N/Z}$
- 自旋截断：$\sigma=0.98\,A^{0.29}\sqrt{E/B_n}$

---

## 第四部分　解析极限验证（任务 #4）

五个解析极限全部通过，验证形式化自洽：

1. **费米气体极限**：ρ ∝ exp(2√aU)/U^(5/4)，恢复 BSFG；
2. **GOE 极限**：能级间距 → Wigner-Dyson $P(s)=\tfrac{\pi}{2}s\,e^{-\pi s^2/4}$；
3. **可积极限**：间距 → 泊松 $P(s)=e^{-s}$；
4. **Fano 极限**：单门态-连续耦合给出 Fano 线形；
5. **单 Breit-Wigner 极限**：孤立共振
   $$\sigma(E)=\pi/k^2\,g_J\frac{\Gamma_n\Gamma_\gamma}{(E-E_r)^2+(\Gamma/2)^2}$$

---

## 第五部分　可证伪判据

| 判据 | 阈值 |
|---|---|
| 自由标量参数数 | ≤ 3 |
| 共振位置偏差 | < 5%（统计平均量） |
| χ²/dof | ≤ 2 |
| 外推核（不重新拟合） | 判据仍成立 |

**基准核**：Au-197（单同位素、国际标准、共振数据跨约 7 个数量级）；对照 Ag-109、U-238，后扩展 Pt/Pb-208/Mn-55/Eu-151，Ta-181 盲测（未做）。

## 产物
- 理论形式化（本文）
- 代码包：常量/拓扑/集体模/伍德萨克森/统计谱/门态/宽度/截面/系综 等模块
