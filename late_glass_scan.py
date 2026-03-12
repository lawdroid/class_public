#!/usr/bin/env python3
"""
Late Glass grid scan: κ_c values from 0.85 to 1.0.
Computes χ² against Planck 2018 TT for each.
"""

import numpy as np
import subprocess
import sys
from pathlib import Path

GD_DIR = Path(__file__).parent
T_CMB = 2.7255
TCMB2 = (T_CMB * 1e6)**2

def load_planck():
    data = np.loadtxt(GD_DIR / 'data' / 'planck_2018_TT.txt')
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]

def load_class_cl(filepath):
    data = np.loadtxt(filepath)
    ell = data[:, 0]
    dl_tt = data[:, 1] * TCMB2
    return ell, dl_tt

def load_background(filepath):
    data = np.loadtxt(filepath)
    return {'z': data[:, 0], 'H': data[:, 3], 'sound_horizon': data[:, 7]}

def get_val_at_z(arr_z, arr_val, z_target):
    idx = np.argmin(np.abs(arr_z - z_target))
    return arr_val[idx]

def compute_chi2(p_ell, p_dl, p_err_m, p_err_p, m_ell, m_dl, ell_min=30):
    mask = p_ell >= ell_min
    p_e, p_d, p_em, p_ep = p_ell[mask], p_dl[mask], p_err_m[mask], p_err_p[mask]
    m_interp = np.interp(p_e, m_ell, m_dl)
    residual = m_interp - p_d
    err = np.where(residual > 0, p_ep, p_em)
    chi2 = np.sum((residual / err)**2)
    dof = len(p_e) - 7
    return chi2, dof

def run_class(kappa_c, z_onset=1100, z_freeze=0, beta=0.5):
    tag = f"lg_k{kappa_c:.3f}"
    ini = f"""root = output/{tag}_
write_background = yes
output = tCl,pCl,lCl
h = 0.6736
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
background_verbose = 1
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
        print(f"  CLASS FAILED for kappa_c={kappa_c}: {result.stderr[-200:]}")
        return None

    cl_file = GD_DIR / f"output/{tag}_00_cl.dat"
    bg_file = GD_DIR / f"output/{tag}_00_background.dat"

    ell, dl = load_class_cl(cl_file)
    bg = load_background(bg_file)

    rs = get_val_at_z(bg['z'], bg['sound_horizon'], 1089)
    H0 = bg['H'][np.argmin(np.abs(bg['z']))]

    return {'ell': ell, 'dl': dl, 'rs': rs, 'H0': H0}


def main():
    print("=" * 70)
    print("LATE GLASS MODEL — Grid Scan vs Planck 2018 TT")
    print("κ = 1.0 before recombination, κ_c after, decaying back to 1.0")
    print("=" * 70)

    p_ell, p_dl, p_err_m, p_err_p = load_planck()
    print(f"Planck data: {len(p_ell)} points, ℓ = {p_ell[0]:.0f}–{p_ell[-1]:.0f}")

    # ΛCDM baseline
    print("\nRunning ΛCDM baseline...")
    lcdm = run_class(1.0)
    if lcdm is None:
        sys.exit(1)
    chi2_lcdm, dof = compute_chi2(p_ell, p_dl, p_err_m, p_err_p, lcdm['ell'], lcdm['dl'])
    print(f"  ΛCDM: χ²/dof = {chi2_lcdm/dof:.4f}  (χ²={chi2_lcdm:.1f}, dof={dof})")
    print(f"  r_s = {lcdm['rs']:.2f} Mpc, H0 = {lcdm['H0']:.2f} km/s/Mpc")

    # Late Glass scan
    kappa_values = [0.85, 0.90, 0.92, 0.94, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999]

    print(f"\n{'κ_c':>8s}  {'H0':>8s}  {'r_s':>8s}  {'χ²/dof':>10s}  {'Δχ²':>8s}  {'Status':>12s}")
    print("-" * 70)
    print(f"{'1.000':>8s}  {lcdm['H0']:8.2f}  {lcdm['rs']:8.2f}  {chi2_lcdm/dof:10.4f}  {'0.0':>8s}  {'ΛCDM':>12s}")

    for kc in kappa_values:
        result = run_class(kc)
        if result is None:
            continue
        chi2, dof = compute_chi2(p_ell, p_dl, p_err_m, p_err_p, result['ell'], result['dl'])
        dchi2 = chi2 - chi2_lcdm

        status = ""
        if chi2/dof < 1.3:
            status = "publishable"
        elif chi2/dof < 1.5:
            status = "marginal"
        else:
            status = "excluded"

        print(f"{kc:8.3f}  {result['H0']:8.2f}  {result['rs']:8.2f}  {chi2/dof:10.4f}  {dchi2:+8.1f}  {status:>12s}")

    print("\n" + "=" * 70)
    print("NOTE: Late Glass preserves r_s (κ=1.0 pre-recombination).")
    print("H0 shift comes from post-recombination angular diameter distance change.")
    print("=" * 70)


if __name__ == '__main__':
    main()
