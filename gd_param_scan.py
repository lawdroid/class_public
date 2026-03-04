#!/usr/bin/env python3
"""
GD-CLASS Parameter Scan: Fit GD model to Planck 2018 TT data.

Strategy: background-only GD (Phase C disabled, justified by omega_BD=50000).
Scan over (h, omega_cdm) with GD kappa=0.85 quench transition.
Compare C_l to Planck data via chi-squared.

Usage:
  python3 gd_param_scan.py          # Run coarse + fine scan
  python3 gd_param_scan.py --fine   # Skip coarse, run fine around known minimum
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess
import sys
import tempfile
import os

GD_DIR = Path(__file__).parent
CLASS_BIN = GD_DIR / 'class'
PLANCK_FILE = GD_DIR / 'data' / 'planck_2018_TT.txt'
OUTPUT_DIR = GD_DIR / 'output'
T_CMB = 2.7255
TCMB2 = (T_CMB * 1e6)**2  # (muK)^2


def load_planck():
    """Load Planck 2018 TT. Returns ell, D_l, sigma."""
    data = np.loadtxt(PLANCK_FILE)
    ell = data[:, 0]
    dl = data[:, 1]
    err_m = data[:, 2]
    err_p = data[:, 3]
    sigma = (err_m + err_p) / 2.0
    return ell, dl, sigma


def load_class_cl(filepath):
    """Load CLASS Cl. Returns ell, D_l (muK^2)."""
    data = np.loadtxt(filepath)
    ell = data[:, 0]
    cl_tt = data[:, 1]  # dimensionless l(l+1)/2pi C_l
    dl_tt = cl_tt * TCMB2
    return ell, dl_tt


def make_ini_content(h, omega_b, omega_cdm, A_s, n_s, tau_reio,
                     kappa_c, z_freeze, z_onset, beta, root):
    """Generate .ini content for a GD CLASS run."""
    return f"""# GD Parameter Scan — auto-generated
root = {root}
write_background = yes
write_thermodynamics = no
output = tCl,pCl,lCl

h = {h}
T_cmb = 2.7255
omega_b = {omega_b}
omega_cdm = {omega_cdm}
N_ur = 2.0328
N_ncdm = 1
m_ncdm = 0.06

reio_parametrization = reio_camb
tau_reio = {tau_reio}

P_k_ini type = analytic_Pk
n_s = {n_s}
A_s = {A_s}
k_pivot = 0.05

gd_kappa_c = {kappa_c}
gd_z_freeze = {z_freeze}
gd_z_onset = {z_onset}
gd_beta = {beta}

l_max_scalars = 2500
background_verbose = 0
thermodynamics_verbose = 0
perturbations_verbose = 0
transfer_verbose = 0
primordial_verbose = 0
harmonic_verbose = 0
lensing_verbose = 0
output_verbose = 0
"""


def run_class(params, tag):
    """Run CLASS with given params. Returns (ell, dl) or None on failure."""
    root = f'output/scan_{tag}_'
    ini_content = make_ini_content(root=root, **params)
    ini_path = GD_DIR / f'_scan_{tag}.ini'

    with open(ini_path, 'w') as f:
        f.write(ini_content)

    try:
        result = subprocess.run(
            [str(CLASS_BIN), str(ini_path)],
            cwd=str(GD_DIR),
            capture_output=True, text=True, timeout=60
        )
        if result.returncode != 0:
            print(f"  FAIL [{tag}]: {result.stderr[:100]}")
            return None

        # Find the output cl file
        cl_file = GD_DIR / f'{root}00_cl.dat'
        if not cl_file.exists():
            # Try incrementing
            for i in range(20):
                cf = GD_DIR / f'{root}{i:02d}_cl.dat'
                if cf.exists():
                    cl_file = cf
                    break

        if not cl_file.exists():
            print(f"  FAIL [{tag}]: no output file")
            return None

        ell, dl = load_class_cl(cl_file)

        # Also get r_s from background file
        bg_file = str(cl_file).replace('_cl.dat', '_background.dat')
        rs = None
        if os.path.exists(bg_file):
            bg = np.loadtxt(bg_file)
            z_col = bg[:, 0]
            rs_col = bg[:, 7]  # sound horizon column
            idx = np.argmin(np.abs(z_col - 1089))
            rs = rs_col[idx]

        return ell, dl, rs

    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT [{tag}]")
        return None
    finally:
        ini_path.unlink(missing_ok=True)


def compute_chi2(model_ell, model_dl, p_ell, p_dl, p_sigma, l_min=30, l_max=2500):
    """Compute chi-squared between model and Planck."""
    mask = (p_ell >= l_min) & (p_ell <= l_max)
    ell_sub = p_ell[mask]
    dl_sub = p_dl[mask]
    sig_sub = p_sigma[mask]
    model_interp = np.interp(ell_sub, model_ell, model_dl)
    chi2 = np.sum(((model_interp - dl_sub) / sig_sub) ** 2)
    return chi2, len(ell_sub)


def run_scan(param_grid, fixed_params, p_ell, p_dl, p_sigma, label="scan"):
    """Run a grid scan. Returns results dict."""
    results = []
    total = len(param_grid)
    for i, point in enumerate(param_grid):
        params = {**fixed_params, **point}
        tag = f"{label}_{i:03d}"
        print(f"  [{i+1}/{total}] h={params['h']:.4f} omega_cdm={params['omega_cdm']:.5f}", end='')

        out = run_class(params, tag)
        if out is None:
            print(" -> FAILED")
            results.append({**point, 'chi2': np.inf, 'ndof': 0, 'rs': None})
            continue

        ell, dl, rs = out
        chi2, ndof = compute_chi2(ell, dl, p_ell, p_dl, p_sigma)
        chi2_red = chi2 / ndof if ndof > 0 else np.inf
        print(f" -> chi2/dof = {chi2_red:.3f}  (rs={rs:.2f} Mpc)" if rs else f" -> chi2/dof = {chi2_red:.3f}")
        results.append({**point, 'chi2': chi2, 'ndof': ndof, 'rs': rs, 'ell': ell, 'dl': dl})

    return results


def cleanup_scan_files():
    """Remove temporary scan output files."""
    for f in OUTPUT_DIR.glob('scan_*'):
        f.unlink()


def main():
    print("=" * 70)
    print("GD-CLASS Parameter Scan: Fitting to Planck 2018 TT")
    print("=" * 70)

    # Load Planck
    p_ell, p_dl, p_sigma = load_planck()
    print(f"Planck data: {len(p_ell)} multipoles (l={p_ell[0]:.0f}-{p_ell[-1]:.0f})")

    # Fixed GD parameters
    gd_params = {
        'kappa_c': 0.85,
        'z_freeze': 1040,
        'z_onset': 1140,
        'beta': 0.5,
    }

    # Fixed cosmological parameters (will vary h and omega_cdm)
    fixed_cosmo = {
        'omega_b': 0.02237,
        'n_s': 0.9649,
        'A_s': 2.1e-9,
        'tau_reio': 0.0544,
    }

    fixed_params = {**fixed_cosmo, **gd_params}

    # ========== PHASE 1: COARSE SCAN ==========
    print("\n" + "=" * 70)
    print("PHASE 1: COARSE SCAN (h x omega_cdm)")
    print("=" * 70)

    h_values = [0.68, 0.70, 0.72, 0.73, 0.74, 0.76]
    ocdm_values = [0.10, 0.11, 0.12, 0.13, 0.14, 0.15]

    coarse_grid = [{'h': h, 'omega_cdm': oc}
                   for h in h_values for oc in ocdm_values]

    coarse_results = run_scan(coarse_grid, fixed_params, p_ell, p_dl, p_sigma, "coarse")

    # Find best
    valid = [r for r in coarse_results if r['chi2'] < np.inf]
    if not valid:
        print("ERROR: All runs failed!")
        sys.exit(1)

    best_coarse = min(valid, key=lambda r: r['chi2'])
    print(f"\nBest coarse: h={best_coarse['h']:.4f}, omega_cdm={best_coarse['omega_cdm']:.5f}")
    print(f"  chi2/dof = {best_coarse['chi2']/best_coarse['ndof']:.4f}")
    if best_coarse['rs']:
        print(f"  r_s(z=1089) = {best_coarse['rs']:.2f} Mpc")

    # Print coarse grid table
    print("\nCoarse grid (chi2/dof):")
    print(f"{'h \\ omega_cdm':>14}", end='')
    for oc in ocdm_values:
        print(f"  {oc:>8.4f}", end='')
    print()
    for h in h_values:
        print(f"  {h:>12.4f}", end='')
        for oc in ocdm_values:
            r = next((x for x in coarse_results if x['h'] == h and x['omega_cdm'] == oc), None)
            if r and r['chi2'] < np.inf:
                print(f"  {r['chi2']/r['ndof']:>8.3f}", end='')
            else:
                print(f"  {'FAIL':>8}", end='')
        print()

    # ========== PHASE 2: FINE SCAN ==========
    print("\n" + "=" * 70)
    print("PHASE 2: FINE SCAN (around coarse minimum)")
    print("=" * 70)

    h_best = best_coarse['h']
    oc_best = best_coarse['omega_cdm']

    h_fine = [h_best + dh for dh in [-0.02, -0.01, -0.005, 0, 0.005, 0.01, 0.02]]
    oc_fine = [oc_best + doc for doc in [-0.02, -0.01, -0.005, 0, 0.005, 0.01, 0.02]]

    fine_grid = [{'h': h, 'omega_cdm': oc}
                 for h in h_fine for oc in oc_fine]

    fine_results = run_scan(fine_grid, fixed_params, p_ell, p_dl, p_sigma, "fine")

    valid_fine = [r for r in fine_results if r['chi2'] < np.inf]
    best_fine = min(valid_fine, key=lambda r: r['chi2'])
    print(f"\nBest fine: h={best_fine['h']:.4f}, omega_cdm={best_fine['omega_cdm']:.5f}")
    print(f"  chi2/dof = {best_fine['chi2']/best_fine['ndof']:.4f}")
    if best_fine['rs']:
        print(f"  r_s(z=1089) = {best_fine['rs']:.2f} Mpc")

    # Print fine grid
    print("\nFine grid (chi2/dof):")
    print(f"{'h \\ omega_cdm':>14}", end='')
    for oc in oc_fine:
        print(f"  {oc:>8.5f}", end='')
    print()
    for h in h_fine:
        print(f"  {h:>12.4f}", end='')
        for oc in oc_fine:
            r = next((x for x in fine_results if abs(x['h'] - h) < 1e-6 and abs(x['omega_cdm'] - oc) < 1e-6), None)
            if r and r['chi2'] < np.inf:
                print(f"  {r['chi2']/r['ndof']:>8.3f}", end='')
            else:
                print(f"  {'FAIL':>8}", end='')
        print()

    # ========== PHASE 3: SCAN omega_b, A_s, n_s ==========
    print("\n" + "=" * 70)
    print("PHASE 3: SCAN SECONDARY PARAMETERS")
    print("=" * 70)

    best_h = best_fine['h']
    best_oc = best_fine['omega_cdm']

    # Scan omega_b
    print("\n--- omega_b scan ---")
    ob_values = [0.0200, 0.0210, 0.0220, 0.02237, 0.0230, 0.0240, 0.0250]
    ob_grid = [{'h': best_h, 'omega_cdm': best_oc, 'omega_b': ob} for ob in ob_values]
    ob_results = run_scan(ob_grid, {k: v for k, v in fixed_params.items() if k != 'omega_b'},
                          p_ell, p_dl, p_sigma, "ob")
    valid_ob = [r for r in ob_results if r['chi2'] < np.inf]
    best_ob = min(valid_ob, key=lambda r: r['chi2'])
    print(f"Best omega_b = {best_ob['omega_b']:.5f} (chi2/dof = {best_ob['chi2']/best_ob['ndof']:.4f})")

    # Scan A_s
    print("\n--- A_s scan ---")
    As_values = [1.6e-9, 1.8e-9, 2.0e-9, 2.1e-9, 2.2e-9, 2.4e-9, 2.6e-9]
    As_grid = [{'h': best_h, 'omega_cdm': best_oc, 'omega_b': best_ob['omega_b'], 'A_s': As}
               for As in As_values]
    As_results = run_scan(As_grid, {k: v for k, v in fixed_params.items()
                                    if k not in ('omega_b', 'A_s')},
                          p_ell, p_dl, p_sigma, "As")
    valid_As = [r for r in As_results if r['chi2'] < np.inf]
    best_As = min(valid_As, key=lambda r: r['chi2'])
    print(f"Best A_s = {best_As['A_s']:.2e} (chi2/dof = {best_As['chi2']/best_As['ndof']:.4f})")

    # Scan n_s
    print("\n--- n_s scan ---")
    ns_values = [0.94, 0.95, 0.96, 0.9649, 0.97, 0.98, 0.99]
    ns_grid = [{'h': best_h, 'omega_cdm': best_oc, 'omega_b': best_ob['omega_b'],
                'A_s': best_As['A_s'], 'n_s': ns} for ns in ns_values]
    ns_results = run_scan(ns_grid, {k: v for k, v in fixed_params.items()
                                    if k not in ('omega_b', 'A_s', 'n_s')},
                          p_ell, p_dl, p_sigma, "ns")
    valid_ns = [r for r in ns_results if r['chi2'] < np.inf]
    best_ns = min(valid_ns, key=lambda r: r['chi2'])
    print(f"Best n_s = {best_ns['n_s']:.4f} (chi2/dof = {best_ns['chi2']/best_ns['ndof']:.4f})")

    # ========== FINAL BEST-FIT ==========
    print("\n" + "=" * 70)
    print("FINAL BEST-FIT RUN")
    print("=" * 70)

    final_params = {
        'h': best_h,
        'omega_b': best_ob['omega_b'],
        'omega_cdm': best_oc,
        'A_s': best_As['A_s'],
        'n_s': best_ns['n_s'],
        'tau_reio': 0.0544,
        **gd_params,
    }

    print(f"Parameters:")
    for k, v in final_params.items():
        if k in gd_params:
            continue
        print(f"  {k} = {v}")

    out = run_class(final_params, "final")
    if out is None:
        print("ERROR: Final run failed!")
        sys.exit(1)

    final_ell, final_dl, final_rs = out
    chi2_final, ndof_final = compute_chi2(final_ell, final_dl, p_ell, p_dl, p_sigma)
    print(f"\nFinal chi2/dof = {chi2_final/ndof_final:.4f} ({chi2_final:.1f} / {ndof_final})")
    if final_rs:
        print(f"r_s(z=1089) = {final_rs:.2f} Mpc")
        h0_implied = final_params['h'] * 100
        print(f"H0 = {h0_implied:.1f} km/s/Mpc")

    # ========== LCDM BASELINE ==========
    print("\nRunning LCDM baseline for comparison...")
    lcdm_params = {
        'h': 0.6736,
        'omega_b': 0.02237,
        'omega_cdm': 0.1200,
        'A_s': 2.1e-9,
        'n_s': 0.9649,
        'tau_reio': 0.0544,
        'kappa_c': 1.0,
        'z_freeze': 1040,
        'z_onset': 1140,
        'beta': 0.5,
    }
    lcdm_out = run_class(lcdm_params, "lcdm")
    if lcdm_out:
        lcdm_ell, lcdm_dl, lcdm_rs = lcdm_out
        chi2_lcdm, ndof_lcdm = compute_chi2(lcdm_ell, lcdm_dl, p_ell, p_dl, p_sigma)
        print(f"LCDM chi2/dof = {chi2_lcdm/ndof_lcdm:.4f} ({chi2_lcdm:.1f} / {ndof_lcdm})")

    # ========== COMPARISON PLOT ==========
    print("\nGenerating comparison plot...")

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), gridspec_kw={'height_ratios': [3, 1]})

    # Panel 1: TT spectrum
    ax = axes[0]
    ax.errorbar(p_ell, p_dl, yerr=p_sigma, fmt='.', color='gray', alpha=0.3,
                markersize=2, elinewidth=0.5, label='Planck 2018 TT', zorder=1)

    if lcdm_out:
        ax.plot(lcdm_ell, lcdm_dl, 'k-', lw=1.5,
                label=rf'$\Lambda$CDM (h={lcdm_params["h"]}, $\chi^2$/dof={chi2_lcdm/ndof_lcdm:.2f})',
                zorder=2)

    ax.plot(final_ell, final_dl, 'r-', lw=1.5,
            label=rf'GD best-fit (h={best_h}, $\chi^2$/dof={chi2_final/ndof_final:.2f})',
            zorder=3)

    ax.set_xlim(2, 2500)
    ax.set_xlabel(r'Multipole $\ell$')
    ax.set_ylabel(r'$\mathcal{D}_\ell^{TT}$ [$\mu K^2$]')
    ax.set_title(r'CMB TT: GD ($\kappa_c=0.85$, bg-only) best-fit vs Planck 2018')
    ax.legend(loc='upper right', fontsize=10)
    ax.set_xscale('log')

    # Panel 2: Residuals
    ax = axes[1]
    if lcdm_out:
        # GD residual vs Planck
        gd_interp = np.interp(p_ell, final_ell, final_dl)
        lcdm_interp = np.interp(p_ell, lcdm_ell, lcdm_dl)
        mask = p_ell >= 2
        gd_resid = (gd_interp[mask] - p_dl[mask]) / p_sigma[mask]
        lcdm_resid = (lcdm_interp[mask] - p_dl[mask]) / p_sigma[mask]

        # Binned residuals for clarity
        bin_size = 20
        n_bins = len(p_ell[mask]) // bin_size
        for label_name, resid, color in [('LCDM', lcdm_resid, 'black'), ('GD best-fit', gd_resid, 'red')]:
            ell_binned = []
            resid_binned = []
            for b in range(n_bins):
                s = b * bin_size
                e = s + bin_size
                ell_binned.append(np.mean(p_ell[mask][s:e]))
                resid_binned.append(np.mean(resid[s:e]))
            ax.plot(ell_binned, resid_binned, 'o-', color=color, markersize=3, lw=1,
                    label=label_name, alpha=0.7)

    ax.axhline(0, color='gray', ls=':', lw=0.5)
    ax.set_xlim(2, 2500)
    ax.set_xlabel(r'Multipole $\ell$')
    ax.set_ylabel(r'$(D_\ell^{model} - D_\ell^{Planck}) / \sigma$')
    ax.set_title('Residuals vs Planck (binned)')
    ax.legend(loc='upper right', fontsize=10)
    ax.set_xscale('log')

    plt.tight_layout()
    outfile = GD_DIR / 'gd_bestfit_vs_planck.png'
    plt.savefig(outfile, dpi=150)
    print(f"Saved: {outfile}")

    # ========== SUMMARY TABLE ==========
    print("\n" + "=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    print(f"{'Parameter':<20} {'LCDM Planck':>15} {'GD best-fit':>15}")
    print("-" * 50)
    print(f"{'h':<20} {0.6736:>15.4f} {best_h:>15.4f}")
    print(f"{'omega_b':<20} {0.02237:>15.5f} {best_ob['omega_b']:>15.5f}")
    print(f"{'omega_cdm':<20} {0.1200:>15.5f} {best_oc:>15.5f}")
    print(f"{'A_s':<20} {2.1e-9:>15.3e} {best_As['A_s']:>15.3e}")
    print(f"{'n_s':<20} {0.9649:>15.4f} {best_ns['n_s']:>15.4f}")
    print(f"{'tau_reio':<20} {0.0544:>15.4f} {0.0544:>15.4f}")
    print(f"{'kappa_c':<20} {'1.0':>15} {'0.85':>15}")
    print(f"{'H0 [km/s/Mpc]':<20} {67.36:>15.1f} {best_h*100:>15.1f}")
    if final_rs:
        print(f"{'r_s [Mpc]':<20} {lcdm_rs:>15.2f} {final_rs:>15.2f}" if lcdm_out and lcdm_rs else "")
    print(f"{'chi2/dof':<20} {chi2_lcdm/ndof_lcdm:>15.4f} {chi2_final/ndof_final:>15.4f}"
          if lcdm_out else f"{'chi2/dof':<20} {'N/A':>15} {chi2_final/ndof_final:>15.4f}")

    # Clean up scan outputs
    cleanup_scan_files()
    print("\nDone. Scan output files cleaned up.")


if __name__ == '__main__':
    main()
