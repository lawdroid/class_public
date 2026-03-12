#!/usr/bin/env python3
"""
Late Glass 2D grid scan: κ_c × h
For each (κ_c, h), run CLASS and compute χ² against Planck 2018 TT.
"""

import numpy as np
import subprocess
import sys
from pathlib import Path

GD_DIR = Path(__file__).parent
T_CMB = 2.7255
TCMB2 = (T_CMB * 1e6)**2
c_kms = 299792.458

def load_planck():
    data = np.loadtxt(GD_DIR / 'data' / 'planck_2018_TT.txt')
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]

def load_class_cl(filepath):
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 1] * TCMB2

def compute_chi2(p_ell, p_dl, p_err_m, p_err_p, m_ell, m_dl, ell_min=30):
    mask = p_ell >= ell_min
    p_e, p_d, p_em, p_ep = p_ell[mask], p_dl[mask], p_err_m[mask], p_err_p[mask]
    m_interp = np.interp(p_e, m_ell, m_dl)
    residual = m_interp - p_d
    err = np.where(residual > 0, p_ep, p_em)
    return np.sum((residual / err)**2), len(p_e) - 7

def run_class(kappa_c, h, z_onset=1100, z_freeze=0, beta=0.5):
    tag = f"lg2d_k{kappa_c:.3f}_h{h:.4f}"
    ini = f"""root = output/{tag}_
output = tCl,pCl,lCl
h = {h}
T_cmb = 2.7255
omega_b = 0.02237
omega_cdm = 0.1200
N_ur = 2.0328
N_ncdm = 1
m_ncdm = 0.06
reio_parametrization = reio_camb
tau_reio = 0.0544
P_k_ini type = analytic_Pk
n_s = 0.9649
A_s = 2.1e-9
k_pivot = 0.05
gd_kappa_c = {kappa_c}
gd_z_onset = {z_onset}
gd_z_freeze = {z_freeze}
gd_beta = {beta}
gd_omega_BD = 50000
gd_phi_c = 0.15
l_max_scalars = 2500
background_verbose = 0
"""
    ini_path = GD_DIR / f"_tmp_{tag}.ini"
    ini_path.write_text(ini)

    result = subprocess.run(
        [str(GD_DIR / 'class'), str(ini_path)],
        capture_output=True, text=True, cwd=str(GD_DIR),
        timeout=120
    )
    ini_path.unlink(missing_ok=True)

    if result.returncode != 0:
        return None

    cl_file = GD_DIR / f"output/{tag}_00_cl.dat"
    if not cl_file.exists():
        return None
    ell, dl = load_class_cl(cl_file)
    return ell, dl

def main():
    print("=" * 75)
    print("LATE GLASS MODEL — 2D Grid Scan: κ_c × h")
    print("κ = 1.0 before z=1100, κ_c after, decaying to 1.0 by z=0")
    print("=" * 75)

    p_ell, p_dl, p_err_m, p_err_p = load_planck()

    # ΛCDM baseline
    result = run_class(1.0, 0.6736)
    chi2_lcdm, dof = compute_chi2(p_ell, p_dl, p_err_m, p_err_p, result[0], result[1])
    print(f"\nΛCDM baseline: χ²/dof = {chi2_lcdm/dof:.4f} (χ²={chi2_lcdm:.1f}, H0=67.36)")

    # 2D scan
    kappa_values = [0.85, 0.90, 0.94, 0.96, 0.98, 0.99, 0.995, 1.0]
    h_values = [0.6736, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73]

    print(f"\n{'κ_c':>6s}", end="")
    for h in h_values:
        H0 = h * 100
        print(f"  H0={H0:.1f}", end="")
    print()
    print("-" * (8 + 10 * len(h_values)))

    best_chi2 = {}
    for kc in kappa_values:
        print(f"{kc:6.3f}", end="", flush=True)
        for h in h_values:
            H0 = h * 100
            result = run_class(kc, h)
            if result is None:
                print(f"  {'FAIL':>8s}", end="")
                continue
            chi2, _ = compute_chi2(p_ell, p_dl, p_err_m, p_err_p, result[0], result[1])
            dchi2 = chi2 - chi2_lcdm
            ratio = chi2 / dof

            marker = ""
            if ratio < 1.2:
                marker = " *"
            elif ratio < 1.3:
                marker = " ."

            print(f"  {ratio:7.3f}{marker}", end="")

            key = f"k{kc}_h{h}"
            best_chi2[key] = (kc, H0, ratio, dchi2)

        print()

    # Find best fit for each H0
    print(f"\n{'Best fit at each H0':}")
    print(f"{'H0':>6s}  {'κ_c':>6s}  {'χ²/dof':>8s}  {'Δχ²':>8s}")
    print("-" * 35)
    for h in h_values:
        H0 = h * 100
        best_k = None
        best_r = 999
        for kc in kappa_values:
            key = f"k{kc}_h{h}"
            if key in best_chi2:
                _, _, r, d = best_chi2[key]
                if r < best_r:
                    best_r = r
                    best_k = kc
                    best_d = d
        if best_k is not None:
            print(f"{H0:6.1f}  {best_k:6.3f}  {best_r:8.4f}  {best_d:+8.1f}")

    print("\n* = publishable (χ²/dof < 1.2)")
    print(". = marginal (χ²/dof < 1.3)")

if __name__ == '__main__':
    main()
