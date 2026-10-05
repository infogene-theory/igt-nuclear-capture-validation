# -*- coding: utf-8 -*-
"""
calibrate_dref.py — 校准统计区能级间距锚点 Dref。
原理：扫描 Dref，使统计区 log(ENDF/IGT) 中位数 ≈ 0。
Dref 是每核素一个物理量级锚点（S_n 处能级间距 / MACS 单点标定），
能量依赖趋势由模型内生，此处只定绝对量级。
"""
import sys
import numpy as np

def load(f):
    d = np.loadtxt(f); return d[:,0], d[:,1]

def igt_sigma_stat(E, model):
    """复刻 IGT 统计区公式，便于扫描。"""
    g = model.spin_stat_factor()
    A = model.topo.A
    out = np.zeros_like(E, dtype=float)
    from igt_nuclear import compute_4pi_k2
    for i, e in enumerate(E):
        Gn = model.g2bar_eV * np.sqrt(max(e, 1e-9))
        Gt = Gn + model.Ggamma_eV
        k2inv = compute_4pi_k2(max(e, 1e-9), A)
        out[i] = np.pi * k2inv * g * (2.0 * np.pi * Gn * model.Ggamma_eV / (model.Dref_eV * Gt))
    return out

def calibrate(nuclide, endf_file, Emin, Emax):
    from igt_nuclear import NuclearTopology, IGTNuclearModel, default_nuclide_params
    p = default_nuclide_params(nuclide)
    topo = NuclearTopology(A=p["A"], kappa=0.5, Dbar_eV=p["D0"], E0_eV=p.get("E0", 4.9))
    Ee, se = load(endf_file)
    m = (Ee >= Emin) & (Ee <= Emax)
    Ee, se = Ee[m], se[m]
    best = None
    for Dref in np.geomspace(p["Dref"]/20, p["Dref"]*20, 120):
        model = IGTNuclearModel(topo=topo, Ggamma_eV=p["Gg"], g2bar_eV=p["g2"], I=p["I"],
                                Sn_MeV=p.get("Sn", 6.5), Dref_eV=Dref)
        Ei = np.geomspace(Ee.min(), Ee.max(), 400)
        si = igt_sigma_stat(Ei, model)
        sig_on_e = 10**np.interp(np.log10(Ee), np.log10(Ei), np.log10(si))
        logr = np.log10(se / sig_on_e)
        rms = np.sqrt(np.mean(logr**2))
        if best is None or rms < best[0]:
            best = (rms, Dref, np.median(logr))
    rms, Dref, med = best
    print(f"{nuclide}: Dref={Dref:.0f} eV | 对数RMS={rms:.3f} | 中位log(ENDF/IGT)={med:.3f}")

if __name__ == "__main__":
    # 用法: python3 calibrate_dref.py 115In endf_In.txt 5000 200000
    calibrate(sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]))
