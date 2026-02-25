#!/usr/bin/env python3
"""
Compare GD-CLASS Kohlrausch transition with LCDM and Planck 2018 data.

Produces:
  1. TT power spectrum: LCDM vs KR-GD vs Planck data
  2. Residuals (model - LCDM) / LCDM
  3. Low-ell ISW zoom
  4. Sound horizon & H0 diagnostics
  5. Beta parameter scan (if multiple beta outputs exist)

Usage:
  python3 compare_planck.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess
import sys

# ----- Configuration -----
GD_DIR = Path(__file__).parent
T_CMB = 2.7255  # K
TCMB2 = (T_CMB * 1e6)**2  # (muK)^2 conversion factor

CONFIGS = {
    'LCDM': {
        'ini': 'lcdm_test.ini',
        'cl': 'output/lcdm_test_00_cl.dat',
        'bg': 'output/lcdm_test_00_background.dat',
        'color': 'black',
        'ls': '-',
        'label': r'$\Lambda$CDM (Planck 2018 best-fit)',
    },
    'Compliant': {
        'ini': 'gd_test_compliant.ini',
        'cl': 'output/gd_compliant_00_cl.dat',
        'bg': 'output/gd_compliant_00_background.dat',
        'color': 'blue',
        'ls': '-',
        'label': r'GD Compliant $\kappa_c=0.85$, $\beta=0.5$',
    },
}


def load_class_cl(filepath):
    """Load CLASS Cl output. Returns ell, D_l (muK^2)."""
    data = np.loadtxt(filepath)
    ell = data[:, 0]
    cl_tt = data[:, 1]  # dimensionless l(l+1)/2pi C_l
    dl_tt = cl_tt * TCMB2
    return ell, dl_tt


def load_planck_data(filepath):
    """Load Planck 2018 TT power spectrum. Returns ell, D_l, err_minus, err_plus."""
    data = np.loadtxt(filepath)
    ell = data[:, 0]
    dl = data[:, 1]
    err_m = data[:, 2]
    err_p = data[:, 3]
    return ell, dl, err_m, err_p


def load_background(filepath):
    """Load background data. Returns dict with z, conf_time, H, sound_horizon."""
    data = np.loadtxt(filepath)
    return {
        'z': data[:, 0],
        'proper_time': data[:, 1],
        'conf_time': data[:, 2],
        'H': data[:, 3],
        'sound_horizon': data[:, 7],
    }


def get_sound_horizon_at_z(bg, z_target):
    """Interpolate sound horizon at given redshift."""
    idx = np.argmin(np.abs(bg['z'] - z_target))
    return bg['sound_horizon'][idx]


def run_beta_scan(betas=[0.3, 0.5, 0.7, 0.9]):
    """Run CLASS with different beta values, return results."""
    results = {}
    for beta in betas:
        ini_name = f'gd_test_KR_beta{beta:.1f}.ini'
        root = f'output/gd_test_KR_beta{beta:.1f}_'

        # Create temporary ini
        ini_path = GD_DIR / ini_name
        with open(GD_DIR / 'gd_test_KR.ini', 'r') as f:
            content = f.read()
        content = content.replace('gd_beta = 0.5', f'gd_beta = {beta}')
        content = content.replace('root = output/gd_test_KR_', f'root = {root}')
        content = content.replace('background_verbose = 2', 'background_verbose = 0')
        with open(ini_path, 'w') as f:
            f.write(content)

        # Run CLASS
        print(f"  Running beta={beta}...")
        result = subprocess.run(
            ['./class', ini_name],
            cwd=str(GD_DIR),
            capture_output=True, text=True, timeout=120
        )
        if result.returncode != 0:
            print(f"  WARNING: CLASS failed for beta={beta}: {result.stderr[:200]}")
            continue

        cl_file = GD_DIR / f'{root}00_cl.dat'
        bg_file = GD_DIR / f'{root}00_background.dat'
        if cl_file.exists():
            ell, dl = load_class_cl(cl_file)
            bg = load_background(bg_file)
            rs = get_sound_horizon_at_z(bg, 1089)
            results[beta] = {'ell': ell, 'dl': dl, 'rs': rs, 'bg': bg}

        # Clean up temp ini
        ini_path.unlink(missing_ok=True)

    return results


def main():
    print("=" * 60)
    print("GD-CLASS: Kohlrausch Transition vs Planck 2018")
    print("=" * 60)

    # ----- Load data -----
    print("\nLoading data...")

    # Planck
    planck_file = GD_DIR / 'data' / 'planck_2018_TT.txt'
    if not planck_file.exists():
        print(f"ERROR: Planck data not found at {planck_file}")
        sys.exit(1)
    p_ell, p_dl, p_err_m, p_err_p = load_planck_data(planck_file)
    print(f"  Planck: {len(p_ell)} multipoles, l={p_ell[0]:.0f}-{p_ell[-1]:.0f}")

    # CLASS models
    models = {}
    for name, cfg in CONFIGS.items():
        cl_file = GD_DIR / cfg['cl']
        bg_file = GD_DIR / cfg['bg']
        if not cl_file.exists():
            print(f"  WARNING: {cl_file} not found, skipping {name}")
            continue
        ell, dl = load_class_cl(cl_file)
        bg = load_background(bg_file)
        models[name] = {'ell': ell, 'dl': dl, 'bg': bg, 'cfg': cfg}
        print(f"  {name}: l={ell[0]:.0f}-{ell[-1]:.0f}")

    if 'LCDM' not in models:
        print("ERROR: LCDM baseline not found")
        sys.exit(1)

    # ----- Diagnostics -----
    print("\n" + "-" * 60)
    print("DIAGNOSTICS")
    print("-" * 60)

    for name, m in models.items():
        bg = m['bg']
        rs_1089 = get_sound_horizon_at_z(bg, 1089)
        rs_1100 = get_sound_horizon_at_z(bg, 1100)
        print(f"\n{name}:")
        print(f"  Sound horizon r_s(z=1089) = {rs_1089:.4f} Mpc")
        print(f"  Sound horizon r_s(z=1100) = {rs_1100:.4f} Mpc")

    # Hubble tension diagnostics
    # Hubble tension diagnostics for any GD model present
    gd_models = {k: v for k, v in models.items() if k != 'LCDM'}
    if 'LCDM' in models and gd_models:
        rs_lcdm = get_sound_horizon_at_z(models['LCDM']['bg'], 1089)
        for gd_name, gd_m in gd_models.items():
            rs_gd = get_sound_horizon_at_z(gd_m['bg'], 1089)
            delta_rs = (rs_gd - rs_lcdm) / rs_lcdm
            h0_input = 67.36  # from .ini
            # H0 ~ 1/r_s at fixed angular scale: smaller r_s => larger H0
            h0_implied = h0_input * rs_lcdm / rs_gd
            print(f"\n  {gd_name}:")
            print(f"  Delta r_s / r_s = {delta_rs*100:.2f}%")
            print(f"  Input H0 = {h0_input} km/s/Mpc")
            print(f"  Implied H0 (from r_s ratio) = {h0_implied:.2f} km/s/Mpc")

    # ----- Run beta scan -----
    print("\n" + "-" * 60)
    print("BETA PARAMETER SCAN")
    print("-" * 60)
    beta_results = run_beta_scan([0.3, 0.5, 0.7, 0.9])
    for beta, res in sorted(beta_results.items()):
        print(f"  beta={beta:.1f}: r_s(z=1089) = {res['rs']:.4f} Mpc")

    # ----- Compute chi2 -----
    print("\n" + "-" * 60)
    print("CHI-SQUARED vs PLANCK (l=30-2500)")
    print("-" * 60)

    # Use symmetric error for chi2
    p_err_sym = (p_err_m + p_err_p) / 2.0

    for name, m in models.items():
        # Interpolate model onto Planck ell values
        mask = (p_ell >= 30) & (p_ell <= m['ell'][-1])
        p_ell_sub = p_ell[mask]
        p_dl_sub = p_dl[mask]
        p_err_sub = p_err_sym[mask]
        model_dl_interp = np.interp(p_ell_sub, m['ell'], m['dl'])
        chi2 = np.sum(((model_dl_interp - p_dl_sub) / p_err_sub) ** 2)
        ndof = len(p_ell_sub)
        print(f"  {name}: chi2 = {chi2:.1f} / {ndof} dof = {chi2/ndof:.3f}")

    # Also chi2 for low-ell (ISW region l=2-30)
    print("\nCHI-SQUARED vs PLANCK (l=2-30, ISW region)")
    for name, m in models.items():
        mask = (p_ell >= 2) & (p_ell <= 30) & (p_ell <= m['ell'][-1])
        p_ell_sub = p_ell[mask]
        p_dl_sub = p_dl[mask]
        p_err_sub = p_err_sym[mask]
        model_dl_interp = np.interp(p_ell_sub, m['ell'], m['dl'])
        chi2 = np.sum(((model_dl_interp - p_dl_sub) / p_err_sub) ** 2)
        ndof = len(p_ell_sub)
        print(f"  {name}: chi2 = {chi2:.1f} / {ndof} dof = {chi2/ndof:.3f}")

    # ----- PLOTS -----
    print("\nGenerating plots...")

    fig, axes = plt.subplots(3, 1, figsize=(14, 16), gridspec_kw={'height_ratios': [3, 1, 2]})

    # --- Panel 1: Full TT spectrum ---
    ax = axes[0]
    # Planck data
    ax.errorbar(p_ell, p_dl, yerr=[p_err_m, p_err_p],
                fmt='.', color='gray', alpha=0.3, markersize=2, elinewidth=0.5,
                label='Planck 2018 TT', zorder=1)

    # Models
    for name, m in models.items():
        cfg = m['cfg']
        ax.plot(m['ell'], m['dl'], color=cfg['color'], ls=cfg['ls'],
                lw=1.5, label=cfg['label'], zorder=2)

    # Beta scan
    beta_colors = {0.3: 'blue', 0.7: 'orange', 0.9: 'green'}
    for beta, res in sorted(beta_results.items()):
        if beta == 0.5:
            continue  # Already plotted as KR_0.5
        color = beta_colors.get(beta, 'purple')
        ax.plot(res['ell'], res['dl'], color=color, ls='--', lw=1.0,
                label=rf'GD KR $\beta={beta}$', zorder=2)

    ax.set_xlim(2, 2500)
    ax.set_xlabel(r'Multipole $\ell$')
    ax.set_ylabel(r'$\mathcal{D}_\ell^{TT}$ [$\mu K^2$]')
    ax.set_title('CMB TT Power Spectrum: GD Kohlrausch vs Planck 2018')
    ax.legend(loc='upper right', fontsize=8)
    ax.set_xscale('log')

    # --- Panel 2: Residuals ---
    ax = axes[1]
    lcdm_ell = models['LCDM']['ell']
    lcdm_dl = models['LCDM']['dl']

    for name, m in models.items():
        if name == 'LCDM':
            continue
        cfg = m['cfg']
        # Compute residual (model - LCDM) / LCDM
        residual = (m['dl'] - lcdm_dl) / lcdm_dl
        ax.plot(m['ell'], residual * 100, color=cfg['color'], ls=cfg['ls'],
                lw=1.5, label=cfg['label'])

    for beta, res in sorted(beta_results.items()):
        if beta == 0.5:
            continue
        color = beta_colors.get(beta, 'purple')
        lcdm_dl_interp = np.interp(res['ell'], lcdm_ell, lcdm_dl)
        residual = (res['dl'] - lcdm_dl_interp) / lcdm_dl_interp
        ax.plot(res['ell'], residual * 100, color=color, ls='--', lw=1.0,
                label=rf'$\beta={beta}$')

    ax.axhline(0, color='black', ls=':', lw=0.5)
    ax.set_xlim(2, 2500)
    ax.set_xlabel(r'Multipole $\ell$')
    ax.set_ylabel(r'$\Delta \mathcal{D}_\ell / \mathcal{D}_\ell^{\Lambda CDM}$ [%]')
    ax.set_title('Residuals vs LCDM')
    ax.legend(loc='upper right', fontsize=8)
    ax.set_xscale('log')

    # --- Panel 3: Low-ell ISW zoom ---
    ax = axes[2]
    # Planck data (low-ell)
    mask_low = p_ell <= 50
    ax.errorbar(p_ell[mask_low], p_dl[mask_low],
                yerr=[p_err_m[mask_low], p_err_p[mask_low]],
                fmt='o', color='gray', alpha=0.5, markersize=4, elinewidth=1,
                label='Planck 2018', zorder=1)

    for name, m in models.items():
        cfg = m['cfg']
        mask = m['ell'] <= 50
        ax.plot(m['ell'][mask], m['dl'][mask], color=cfg['color'], ls=cfg['ls'],
                lw=2, label=cfg['label'], zorder=2)

    for beta, res in sorted(beta_results.items()):
        if beta == 0.5:
            continue
        color = beta_colors.get(beta, 'purple')
        mask = res['ell'] <= 50
        ax.plot(res['ell'][mask], res['dl'][mask], color=color, ls='--', lw=1.5,
                label=rf'$\beta={beta}$', zorder=2)

    ax.set_xlim(2, 50)
    ax.set_xlabel(r'Multipole $\ell$')
    ax.set_ylabel(r'$\mathcal{D}_\ell^{TT}$ [$\mu K^2$]')
    ax.set_title(r'Low-$\ell$ ISW Region (l=2-50)')
    ax.legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    outfile = GD_DIR / 'compare_planck.png'
    plt.savefig(outfile, dpi=150)
    print(f"\nSaved: {outfile}")

    # ----- Summary table -----
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"{'Model':<30} {'r_s [Mpc]':>10} {'D_l(l=2)':>12} {'D_l(l=220)':>12} {'Implied H0':>12}")
    print("-" * 76)
    rs_lcdm = get_sound_horizon_at_z(models['LCDM']['bg'], 1089)
    for name, m in models.items():
        rs = get_sound_horizon_at_z(m['bg'], 1089)
        dl_2 = m['dl'][0]  # l=2
        idx_220 = np.argmin(np.abs(m['ell'] - 220))
        dl_220 = m['dl'][idx_220]
        h0_implied = 67.36 * rs_lcdm / rs
        print(f"{name:<30} {rs:>10.4f} {dl_2:>12.1f} {dl_220:>12.1f} {h0_implied:>12.1f}")

    for beta, res in sorted(beta_results.items()):
        name = f"KR beta={beta}"
        dl_2 = res['dl'][0]
        idx_220 = np.argmin(np.abs(res['ell'] - 220))
        dl_220 = res['dl'][idx_220]
        h0_implied = 67.36 * rs_lcdm / res['rs']
        print(f"{name:<30} {res['rs']:>10.4f} {dl_2:>12.1f} {dl_220:>12.1f} {h0_implied:>12.1f}")

    print("\nDone.")


if __name__ == '__main__':
    main()
