# -*- coding: utf-8 -*-
"""
compare_igt_endf.py — 统计区对比 IGT 预测 vs ENDF/B-VI 平滑截面。
用法: python3 compare_igt_endf.py <nuclide>
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NUCLIDES = {
    "197Au": dict(igt="igt_predict_197Au.txt", endf="endf_Au.txt", lo=5e3,  hi=2e5,
                  label="197Au(n,γ)", color="#1f77b4"),
    "240Pu": dict(igt="igt_predict_240Pu.txt", endf="endf_Pu.txt", lo=5e3, hi=2e5,
                  label="240Pu(n,γ)", color="#d62728"),
    "115In": dict(igt="igt_predict_115In.txt", endf="endf_In.txt", lo=5e3, hi=2e5,
                  label="115In(n,γ)", color="#9467bd"),
    "56Fe":  dict(igt="igt_predict_56Fe.txt",  endf="endf_Fe.txt", lo=4e5, hi=2e6,
                  label="56Fe(n,γ)", color="#8c564b"),
}

def load(f):
    d = np.loadtxt(f)
    return d[:, 0], d[:, 1]

def main(nuclide):
    cfg = NUCLIDES[nuclide]
    Ei, si = load(cfg["igt"])
    Ee, se = load(cfg["endf"])

    # 统计区窗口: 取 IGT 与 ENDF 共同覆盖的能量段
    lo = cfg["lo"]
    hi = min(Ei.max(), Ee.max())
    mi = (Ei >= lo) & (Ei <= hi)
    me = (Ee >= lo) & (Ee <= hi)
    Ei, si = Ei[mi], si[mi]
    Ee, se = Ee[me], se[me]

    # 在 ENDF 能量点上内插 IGT，算比值
    sig_i_on_e = np.interp(np.log10(Ee), np.log10(Ei), np.log10(si))
    sig_i_on_e = 10 ** sig_i_on_e
    ratio = se / sig_i_on_e
    med = np.median(ratio)
    geo = np.exp(np.mean(np.log(ratio)))      # 几何平均比
    l10 = np.log10(ratio)
    rms = np.sqrt(np.mean(l10**2))            # 对数域均方根偏差

    print(f"=== {nuclide} 统计区对比 (E∈[{lo:.0e},{hi:.1e}] eV, ENDF {len(Ee)} 点) ===")
    print(f"IGT中位σ: {np.median(si):.3g} b | ENDF中位σ: {np.median(se):.3g} b")
    print(f"比值 ENDF/IGT: 几何均值 {geo:.2f} | 中位数 {med:.2f} | 对数RMS {rms:.2f}")
    print(f"说明: 比值≈1 表示IGT统计区量级吻合; RMS<0.5 表示能量依赖趋势接近")
    print(f"  (>1 偏大, <1 偏小; |log10|>0.5 即差>3倍)")

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True,
                             gridspec_kw={'height_ratios': [3, 1]})
    ax = axes[0]
    ax.loglog(Ei, si, lw=1.8, color=cfg["color"], label=f"IGT 预测 {cfg['label']} (统计区)")
    ax.loglog(Ee, se, 'o', ms=3.5, color="#2ca02c", alpha=0.85,
              label=f"ENDF/B-VI 平滑 {cfg['label']}")
    ax.set_ylabel("俘获截面 σ$_γ$ [barn]")
    ax.set_title(f"IGT vs ENDF 统计区对比 — {cfg['label']}（跨核素迁移检验）")
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend()

    ax2 = axes[1]
    ax2.semilogx(Ee, ratio, 'o', ms=3, color="#7f7f7f", label="ENDF/IGT 比值")
    ax2.axhline(1.0, color="k", lw=1)
    ax2.axhline(10**rms, color="r", ls="--", lw=1, label=f"±10^{rms:.2f}(对数RMS)")
    ax2.axhline(10**(-rms), color="r", ls="--", lw=1)
    ax2.set_ylabel("比值 ENDF/IGT")
    ax2.set_xlabel("入射中子能量 E$_n$ [eV]")
    ax2.set_yscale("log")
    ax2.grid(True, which="both", ls=":", alpha=0.4)
    ax2.legend(loc="best")
    fig.tight_layout()
    out = f"compare_{nuclide}.png"
    fig.savefig(out, dpi=150)
    print(f"对比图: {out}")

if __name__ == "__main__":
    main(sys.argv[1])
