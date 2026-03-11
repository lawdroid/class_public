#!/usr/bin/env python3
"""
gd_paper.py — Publication analysis pipeline for
"Glassy Dynamics and the Hubble Tension"

Produces all figures, tables, and numerical results for the paper.
Uses cached CLASS outputs to enable incremental development.

Usage:
  python3 gd_paper.py              # Full pipeline
  python3 gd_paper.py --figures    # Regenerate figures from cached data
  python3 gd_paper.py --scan-only  # Run kappa scan, save cache, skip figures
  python3 gd_paper.py --clean      # Clear cache and re-run everything
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess
import json
import hashlib
import argparse
import sys

# ============================================================
# SECTION 0: Constants, paths, configuration
# ============================================================

GD_DIR = Path(__file__).parent
CLASS_BIN = GD_DIR / 'class'
PLANCK_FILE = GD_DIR / 'data' / 'planck_2018_TT.txt'
OUTPUT_DIR = GD_DIR / 'output'
CACHE_DIR = GD_DIR / 'paper_cache'
FIGURES_DIR = GD_DIR / 'paper_figures'

T_CMB = 2.7255
TCMB2 = (T_CMB * 1e6)**2

# GD fixed parameters (recombination quench)
GD_TRANSITION = {'z_freeze': 1040, 'z_onset': 1140, 'beta': 0.5}

# Planck 2018 LCDM fiducial
LCDM_FIDUCIAL = {
    'h': 0.6736, 'omega_b': 0.02237, 'omega_cdm': 0.1200,
    'A_s': 2.1e-9, 'n_s': 0.9649, 'tau_reio': 0.0544,
}

# SH0ES measurement
H0_SHOES = 73.04
H0_SHOES_ERR = 1.04

# Kappa scan values
KAPPA_VALUES = np.round(np.arange(0.930, 1.001, 0.005), 3)


def setup_pub_style():
    """Publication matplotlib style."""
    try:
        plt.rcParams.update({'text.usetex': True, 'font.family': 'serif'})
        # Test if LaTeX works
        fig, ax = plt.subplots()
        ax.set_title(r'$\kappa$')
        fig.savefig('/dev/null', format='png')
        plt.close(fig)
    except Exception:
        plt.rcParams.update({'text.usetex': False, 'mathtext.fontset': 'cm'})

    plt.rcParams.update({
        'font.size': 12, 'axes.labelsize': 14, 'axes.titlesize': 14,
        'xtick.labelsize': 11, 'ytick.labelsize': 11, 'legend.fontsize': 10,
        'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    })


# ============================================================
# SECTION 1: Infrastructure
# ============================================================

def load_planck():
    """Load Planck 2018 TT. Returns (ell, D_l, sigma)."""
    data = np.loadtxt(PLANCK_FILE)
    ell = data[:, 0]
    dl = data[:, 1]
    sigma = (data[:, 2] + data[:, 3]) / 2.0
    return ell, dl, sigma


def load_class_cl(filepath):
    """Load CLASS Cl output. Returns (ell, D_l_muK2)."""
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 1] * TCMB2


def load_background(filepath):
    """Load CLASS background output. Returns dict with arrays."""
    data = np.loadtxt(filepath)
    result = {
        'z': data[:, 0],
        'conf_time': data[:, 2],
        'H': data[:, 3],
        'comov_dist': data[:, 4],
        'ang_diam_dist': data[:, 5],
        'sound_horizon': data[:, 7],
    }
    # GD columns (if present, after standard 23 columns)
    if data.shape[1] > 23:
        result['gd_kappa'] = data[:, 23]
    return result


def get_rs_at_z(bg, z_target):
    """Interpolate sound horizon at z_target."""
    return np.interp(z_target, bg['z'][::-1], bg['sound_horizon'][::-1])


def get_DA_at_z(bg, z_target):
    """Interpolate angular diameter distance at z_target."""
    return np.interp(z_target, bg['z'][::-1], bg['ang_diam_dist'][::-1])


def compute_chi2(model_ell, model_dl, p_ell, p_dl, p_sigma, l_min=30, l_max=2500):
    """Chi-squared between model and Planck. Returns (chi2, ndof)."""
    mask = (p_ell >= l_min) & (p_ell <= l_max)
    ell_sub = p_ell[mask]
    sig_sub = p_sigma[mask]
    model_interp = np.interp(ell_sub, model_ell, model_dl)
    chi2 = np.sum(((model_interp - p_dl[mask]) / sig_sub) ** 2)
    return chi2, int(np.sum(mask))


def compute_chi2_by_range(model_ell, model_dl, p_ell, p_dl, p_sigma):
    """Chi2 breakdown by l-range."""
    ranges = [(2, 29), (30, 800), (801, 2500), (2, 2500)]
    result = {}
    for l_min, l_max in ranges:
        chi2, ndof = compute_chi2(model_ell, model_dl, p_ell, p_dl, p_sigma, l_min, l_max)
        result[(l_min, l_max)] = (chi2, ndof)
    return result


def kappa_of_z(z, kappa_c, z_freeze=1040, z_onset=1140, beta=0.5):
    """Replicate kappa(z) from background.c exactly."""
    z = np.atleast_1d(z).astype(float)
    kappa = np.ones_like(z)
    for i, zi in enumerate(z):
        if zi >= z_onset:
            kappa[i] = kappa_c
        elif zi <= z_freeze:
            kappa[i] = 1.0
        else:
            t = (zi - z_freeze) / (z_onset - z_freeze)
            decay = np.exp(-(t / 0.5)**beta)
            kappa[i] = kappa_c + (1.0 - kappa_c) * decay
    return kappa


# --- Cache system ---

def cache_key(params):
    """Deterministic hash of parameter dict."""
    s = json.dumps(params, sort_keys=True, default=str)
    return hashlib.md5(s.encode()).hexdigest()[:12]


def make_ini_content(h, omega_b, omega_cdm, A_s, n_s, tau_reio,
                     kappa_c, z_freeze, z_onset, beta, root):
    """Generate CLASS .ini content."""
    return f"""root = {root}
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


def run_class_cached(cosmo_params, kappa_c):
    """Run CLASS with caching. Returns dict with ell, dl, rs, DA, bg_file or None."""
    full_params = {**cosmo_params, 'kappa_c': kappa_c, **GD_TRANSITION}
    key = cache_key(full_params)

    cache_file = CACHE_DIR / f'{key}.json'
    if cache_file.exists():
        with open(cache_file) as f:
            cached = json.load(f)
        return {
            'ell': np.array(cached['ell']),
            'dl': np.array(cached['dl']),
            'rs': cached['rs'],
            'DA': cached.get('DA', None),
            'bg_file': cached.get('bg_file', ''),
        }

    # Run CLASS
    tag = f"p_{key}"
    root = f'output/paper_{tag}_'
    ini_content = make_ini_content(
        root=root, kappa_c=kappa_c,
        z_freeze=GD_TRANSITION['z_freeze'],
        z_onset=GD_TRANSITION['z_onset'],
        beta=GD_TRANSITION['beta'],
        **cosmo_params
    )
    ini_path = GD_DIR / f'_paper_{tag}.ini'
    with open(ini_path, 'w') as f:
        f.write(ini_content)

    try:
        result = subprocess.run(
            [str(CLASS_BIN), str(ini_path)],
            cwd=str(GD_DIR), capture_output=True, text=True, timeout=60
        )
        if result.returncode != 0:
            return None

        # Find output files
        cl_file = GD_DIR / f'{root}00_cl.dat'
        bg_file = GD_DIR / f'{root}00_background.dat'
        if not cl_file.exists():
            return None

        ell, dl = load_class_cl(cl_file)
        bg = load_background(bg_file)
        rs = get_rs_at_z(bg, 1089)
        DA = get_DA_at_z(bg, 1089)

        # Cache result
        CACHE_DIR.mkdir(exist_ok=True)
        with open(cache_file, 'w') as f:
            json.dump({
                'params': {k: str(v) if isinstance(v, float) and v < 1e-6 else v
                           for k, v in full_params.items()},
                'ell': ell.tolist(), 'dl': dl.tolist(),
                'rs': float(rs), 'DA': float(DA),
                'bg_file': str(bg_file),
            }, f)

        return {'ell': ell, 'dl': dl, 'rs': rs, 'DA': DA, 'bg_file': str(bg_file)}

    except subprocess.TimeoutExpired:
        return None
    finally:
        ini_path.unlink(missing_ok=True)


# ============================================================
# SECTION 2: Kappa scan with per-kappa optimization
# ============================================================

def _eval(h, oc, ob, As, ns, kappa_c, p_ell, p_dl, p_sigma):
    """Evaluate chi2 for a parameter set. Returns (chi2, ndof, result) or None."""
    cosmo = {'h': h, 'omega_b': ob, 'omega_cdm': oc,
             'A_s': As, 'n_s': ns, 'tau_reio': 0.0544}
    out = run_class_cached(cosmo, kappa_c)
    if out is None:
        return None
    chi2, ndof = compute_chi2(out['ell'], out['dl'], p_ell, p_dl, p_sigma)
    return chi2, ndof, out


def optimize_for_kappa(kappa_c, p_ell, p_dl, p_sigma, prev_best=None):
    """Find best-fit cosmological params for given kappa_c.
    Iterates grid + secondary scans for convergence.
    Returns dict with all results."""

    # Seed from previous best or LCDM fiducial
    if prev_best is not None:
        best_h = prev_best['h']
        best_oc = prev_best['omega_cdm']
        best_ob = prev_best.get('omega_b', LCDM_FIDUCIAL['omega_b'])
        best_As = prev_best.get('A_s', LCDM_FIDUCIAL['A_s'])
        best_ns = prev_best.get('n_s', LCDM_FIDUCIAL['n_s'])
    else:
        best_h = LCDM_FIDUCIAL['h']
        best_oc = LCDM_FIDUCIAL['omega_cdm']
        best_ob = LCDM_FIDUCIAL['omega_b']
        best_As = LCDM_FIDUCIAL['A_s']
        best_ns = LCDM_FIDUCIAL['n_s']

    best_chi2 = np.inf
    best_ndof = 0

    # Iterate: grid search (h, ocdm) + secondary scans, repeat twice
    for iteration in range(2):
        h0, oc0 = best_h, best_oc

        # Phase 1: 5x5 grid in (h, omega_cdm) using current best secondary params
        h_grid = [h0 + dh for dh in [-0.02, -0.01, 0, 0.01, 0.02]]
        oc_grid = [oc0 + doc for doc in [-0.01, -0.005, 0, 0.005, 0.01]]

        for h in h_grid:
            if h < 0.5 or h > 0.9:
                continue
            for oc in oc_grid:
                if oc < 0.05 or oc > 0.25:
                    continue
                r = _eval(h, oc, best_ob, best_As, best_ns, kappa_c,
                          p_ell, p_dl, p_sigma)
                if r and r[0] < best_chi2:
                    best_chi2, best_ndof = r[0], r[1]
                    best_h, best_oc = h, oc

        # Phase 2: 3x3 refinement around best (h, omega_cdm)
        for dh in [-0.005, 0, 0.005]:
            for doc in [-0.0025, 0, 0.0025]:
                h = best_h + dh
                oc = best_oc + doc
                if h < 0.5 or oc < 0.05:
                    continue
                r = _eval(h, oc, best_ob, best_As, best_ns, kappa_c,
                          p_ell, p_dl, p_sigma)
                if r and r[0] < best_chi2:
                    best_chi2, best_ndof = r[0], r[1]
                    best_h, best_oc = h, oc

        # Phase 3: 1D scans of secondary params using current best (h, ocdm)
        # omega_b scan
        for ob in [0.020, 0.0205, 0.021, 0.0215, 0.022, 0.02237, 0.023, 0.024]:
            r = _eval(best_h, best_oc, ob, best_As, best_ns, kappa_c,
                      p_ell, p_dl, p_sigma)
            if r and r[0] < best_chi2:
                best_chi2, best_ndof = r[0], r[1]
                best_ob = ob

        # A_s scan
        for As in [1.9e-9, 2.0e-9, 2.05e-9, 2.1e-9, 2.15e-9, 2.2e-9]:
            r = _eval(best_h, best_oc, best_ob, As, best_ns, kappa_c,
                      p_ell, p_dl, p_sigma)
            if r and r[0] < best_chi2:
                best_chi2, best_ndof = r[0], r[1]
                best_As = As

        # n_s scan (denser near Planck value)
        for ns in [0.95, 0.955, 0.96, 0.9649, 0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 1.00]:
            r = _eval(best_h, best_oc, best_ob, best_As, ns, kappa_c,
                      p_ell, p_dl, p_sigma)
            if r and r[0] < best_chi2:
                best_chi2, best_ndof = r[0], r[1]
                best_ns = ns

    # Final run with best params
    cosmo_best = {'h': best_h, 'omega_b': best_ob, 'omega_cdm': best_oc,
                  'A_s': best_As, 'n_s': best_ns, 'tau_reio': 0.0544}
    out = run_class_cached(cosmo_best, kappa_c)

    return {
        'kappa_c': kappa_c,
        'h': best_h, 'omega_b': best_ob, 'omega_cdm': best_oc,
        'A_s': best_As, 'n_s': best_ns, 'tau_reio': 0.0544,
        'chi2': best_chi2, 'ndof': best_ndof,
        'chi2_dof': best_chi2 / best_ndof if best_ndof > 0 else np.inf,
        'H0': best_h * 100,
        'rs': out['rs'] if out else None,
        'DA': out['DA'] if out else None,
        'ell': out['ell'].tolist() if out else None,
        'dl': out['dl'].tolist() if out else None,
        'bg_file': out['bg_file'] if out else None,
    }


def run_kappa_scan(kappa_values, p_ell, p_dl, p_sigma):
    """Run full kappa scan with sequential optimization."""
    print(f"\nScanning {len(kappa_values)} kappa values...")
    print(f"{'kappa':>8} {'h':>7} {'ocdm':>8} {'ob':>8} {'ns':>7} {'chi2/dof':>10} {'H0':>6} {'rs':>8}")
    print("-" * 70)

    results = []
    prev_best = None

    # Scan from kappa=1.00 downward
    for kappa_c in sorted(kappa_values, reverse=True):
        r = optimize_for_kappa(kappa_c, p_ell, p_dl, p_sigma, prev_best)
        results.append(r)
        prev_best = r
        rs_str = f"{r['rs']:.2f}" if r['rs'] else "N/A"
        print(f"{kappa_c:>8.3f} {r['h']:>7.4f} {r['omega_cdm']:>8.5f} "
              f"{r['omega_b']:>8.5f} {r['n_s']:>7.4f} {r['chi2_dof']:>10.4f} "
              f"{r['H0']:>6.1f} {rs_str:>8}")

    # Sort by kappa ascending for output
    results.sort(key=lambda r: r['kappa_c'])
    return results


def save_master_results(results):
    """Save scan results to JSON."""
    CACHE_DIR.mkdir(exist_ok=True)
    outfile = CACHE_DIR / 'master_results.json'
    with open(outfile, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"Saved master results: {outfile}")


def load_master_results():
    """Load cached scan results."""
    infile = CACHE_DIR / 'master_results.json'
    with open(infile) as f:
        results = json.load(f)
    # Convert ell/dl back to numpy
    for r in results:
        if r.get('ell'):
            r['ell'] = np.array(r['ell'])
            r['dl'] = np.array(r['dl'])
    return results


# ============================================================
# SECTION 3: Figures
# ============================================================

def figure_1_theory(scan_results):
    """Figure 1: kappa(z) profile + H(z)/H_LCDM(z) ratio."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    z_arr = np.logspace(1, 4, 1000)
    colors = {0.93: '#D32F2F', 0.96: '#FF9800', 0.99: '#1976D2', 1.00: 'black'}
    labels = {0.93: r'$\kappa_c = 0.93$', 0.96: r'$\kappa_c = 0.96$',
              0.99: r'$\kappa_c = 0.99$', 1.00: r'$\Lambda$CDM ($\kappa = 1$)'}

    # Panel (a): kappa(z)
    for kc in [1.00, 0.99, 0.96, 0.93]:
        kappa = kappa_of_z(z_arr, kc)
        ax1.plot(z_arr, kappa, color=colors[kc], lw=2, label=labels[kc])

    ax1.axvline(1089, color='gray', ls='--', lw=0.8, alpha=0.7)
    ax1.axvspan(GD_TRANSITION['z_freeze'], GD_TRANSITION['z_onset'],
                alpha=0.08, color='blue')
    ax1.text(1089, 1.008, r'$z_{\rm rec}$', ha='center', fontsize=10, color='gray')
    ax1.set_xscale('log')
    ax1.set_xlabel(r'Redshift $z$')
    ax1.set_ylabel(r'Spacetime stiffness $\kappa(z)$')
    ax1.set_title('(a) GD stiffness profile')
    ax1.set_xlim(10, 1e4)
    ax1.set_ylim(0.92, 1.015)
    ax1.legend(loc='lower left', fontsize=9)

    # Panel (b): H(z)/H_LCDM(z)
    # Load background files for key models
    lcdm_r = next((r for r in scan_results if abs(r['kappa_c'] - 1.0) < 0.001), None)
    if lcdm_r and lcdm_r.get('bg_file') and Path(lcdm_r['bg_file']).exists():
        bg_lcdm = load_background(lcdm_r['bg_file'])

        for kc in [0.99, 0.96, 0.93]:
            r = next((r for r in scan_results if abs(r['kappa_c'] - kc) < 0.001), None)
            if r and r.get('bg_file') and Path(r['bg_file']).exists():
                bg_gd = load_background(r['bg_file'])
                # Interpolate H_gd onto LCDM z grid
                z_common = bg_lcdm['z']
                mask = (z_common >= 10) & (z_common <= 1e4)
                H_lcdm = bg_lcdm['H'][mask]
                H_gd = np.interp(z_common[mask], bg_gd['z'][::-1], bg_gd['H'][::-1])
                ratio = H_gd / H_lcdm
                ax2.plot(z_common[mask], ratio, color=colors[kc], lw=2, label=labels[kc])

        ax2.axhline(1.0, color='black', ls='-', lw=1)
        ax2.axvline(1089, color='gray', ls='--', lw=0.8, alpha=0.7)
    else:
        ax2.text(0.5, 0.5, 'Background files\nnot available',
                 transform=ax2.transAxes, ha='center', va='center')

    ax2.set_xscale('log')
    ax2.set_xlabel(r'Redshift $z$')
    ax2.set_ylabel(r'$H_{\rm GD}(z) / H_{\Lambda{\rm CDM}}(z)$')
    ax2.set_title('(b) Expansion rate ratio')
    ax2.set_xlim(10, 1e4)
    ax2.legend(loc='upper left', fontsize=9)

    plt.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(FIGURES_DIR / f'fig1_theory.{ext}')
    plt.close(fig)
    print("  Saved fig1_theory")


def figure_2_cmb(scan_results, p_ell, p_dl, p_sigma):
    """Figure 2: TT power spectrum + residuals."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7),
                                    gridspec_kw={'height_ratios': [3, 1]})
    fig.subplots_adjust(hspace=0.05)

    # Planck data
    ax1.errorbar(p_ell, p_dl, yerr=p_sigma, fmt='.', color='lightgray',
                 alpha=0.4, markersize=1.5, elinewidth=0.3,
                 label='Planck 2018 TT', zorder=1)

    # Models
    models = [
        (1.00, 'black', '-', 2.0),
        (0.99, '#1976D2', '-', 2.0),
        (0.97, '#FF9800', '--', 1.5),
    ]

    for kc, color, ls, lw in models:
        r = next((r for r in scan_results if abs(r['kappa_c'] - kc) < 0.001), None)
        if r and r.get('ell') is not None:
            ell = np.array(r['ell']) if not isinstance(r['ell'], np.ndarray) else r['ell']
            dl = np.array(r['dl']) if not isinstance(r['dl'], np.ndarray) else r['dl']
            if kc == 1.0:
                label = rf'$\Lambda$CDM ($H_0={r["H0"]:.1f}$, $\chi^2$/dof={r["chi2_dof"]:.2f})'
            else:
                label = rf'GD $\kappa={kc}$ ($H_0={r["H0"]:.1f}$, $\chi^2$/dof={r["chi2_dof"]:.2f})'
            ax1.plot(ell, dl, color=color, ls=ls, lw=lw, label=label, zorder=3)

    ax1.set_xlim(2, 2500)
    ax1.set_ylim(0, 7000)
    ax1.set_xscale('log')
    ax1.set_ylabel(r'$\mathcal{D}_\ell^{TT}$ [$\mu$K$^2$]')
    ax1.set_title('CMB TT Power Spectrum: Glassy Dynamics vs Planck 2018')
    ax1.legend(loc='upper right', fontsize=9)
    ax1.set_xticklabels([])

    # Residuals panel
    bin_size = 30
    mask = p_ell >= 2

    for kc, color, ls, lw in [(1.00, 'black', '-', 1.5), (0.99, '#1976D2', '-', 1.5)]:
        r = next((r for r in scan_results if abs(r['kappa_c'] - kc) < 0.001), None)
        if r and r.get('ell') is not None:
            ell = np.array(r['ell']) if not isinstance(r['ell'], np.ndarray) else r['ell']
            dl = np.array(r['dl']) if not isinstance(r['dl'], np.ndarray) else r['dl']
            interp = np.interp(p_ell[mask], ell, dl)
            resid = (interp - p_dl[mask]) / p_sigma[mask]
            n_bins = len(resid) // bin_size
            ell_b, res_b = [], []
            for b in range(n_bins):
                s = b * bin_size
                ell_b.append(np.mean(p_ell[mask][s:s+bin_size]))
                res_b.append(np.mean(resid[s:s+bin_size]))
            lbl = r'$\Lambda$CDM' if kc == 1.0 else rf'GD $\kappa={kc}$'
            ax2.plot(ell_b, res_b, 'o-', color=color, markersize=2.5, lw=1,
                     label=lbl, alpha=0.8)

    ax2.axhline(0, color='gray', ls=':', lw=0.5)
    for v in [-1, 1, -2, 2]:
        ax2.axhline(v, color='gray', ls='--', lw=0.3, alpha=0.4)
    ax2.set_xlim(2, 2500)
    ax2.set_ylim(-3, 3)
    ax2.set_xscale('log')
    ax2.set_xlabel(r'Multipole $\ell$')
    ax2.set_ylabel(r'$(D_\ell^{\rm model} - D_\ell^{\rm Planck})/\sigma$')
    ax2.legend(loc='upper left', fontsize=9)

    for ext in ['png', 'pdf']:
        fig.savefig(FIGURES_DIR / f'fig2_cmb_tt.{ext}')
    plt.close(fig)
    print("  Saved fig2_cmb_tt")


def figure_3_constraint(scan_results):
    """Figure 3: H0 vs kappa + chi2/dof vs kappa."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    kappas = [r['kappa_c'] for r in scan_results]
    h0s = [r['H0'] for r in scan_results]
    chi2s = [r['chi2_dof'] for r in scan_results]

    # Panel (a): H0 vs kappa
    ax1.plot(kappas, h0s, 'o-', color='#1976D2', lw=2, markersize=6, zorder=3)
    ax1.axhspan(H0_SHOES - H0_SHOES_ERR, H0_SHOES + H0_SHOES_ERR,
                alpha=0.15, color='red', label=rf'SH0ES $H_0 = {H0_SHOES} \pm {H0_SHOES_ERR}$')
    ax1.axhspan(67.36 - 0.54, 67.36 + 0.54, alpha=0.15, color='gray',
                label=r'Planck $\Lambda$CDM $H_0 = 67.4 \pm 0.5$')
    ax1.set_xlabel(r'$\kappa_c$')
    ax1.set_ylabel(r'$H_0$ [km/s/Mpc]')
    ax1.set_title(r'(a) Hubble constant vs $\kappa_c$')
    ax1.set_xlim(0.925, 1.005)
    ax1.legend(loc='upper left', fontsize=8)
    ax1.grid(True, alpha=0.2)

    # Panel (b): chi2/dof vs kappa
    ax2.plot(kappas, chi2s, 'o-', color='#1976D2', lw=2, markersize=6, zorder=3)

    # 95% CL threshold (Delta chi2 = 3.84 for 1 extra parameter)
    lcdm_r = next((r for r in scan_results if abs(r['kappa_c'] - 1.0) < 0.001), None)
    if lcdm_r:
        chi2_threshold = lcdm_r['chi2'] + 3.84
        chi2dof_threshold = chi2_threshold / lcdm_r['ndof']
        ax2.axhline(chi2dof_threshold, color='red', ls='--', lw=1.5,
                     label=rf'95\% CL ($\Delta\chi^2 = 3.84$)')
        ax2.axhline(lcdm_r['chi2_dof'], color='black', ls=':', lw=1,
                     label=rf'$\Lambda$CDM $\chi^2$/dof = {lcdm_r["chi2_dof"]:.3f}')

        # Find kappa at 95% CL bound
        for r in sorted(scan_results, key=lambda x: x['kappa_c'], reverse=True):
            if r['chi2'] > chi2_threshold:
                # Interpolate
                kappa_bound = r['kappa_c']
                break

    ax2.set_xlabel(r'$\kappa_c$')
    ax2.set_ylabel(r'$\chi^2$/dof vs Planck TT')
    ax2.set_title(r'(b) Fit quality vs $\kappa_c$')
    ax2.set_xlim(0.925, 1.005)
    ax2.legend(loc='upper left', fontsize=8)
    ax2.grid(True, alpha=0.2)

    plt.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(FIGURES_DIR / f'fig3_constraint.{ext}')
    plt.close(fig)
    print("  Saved fig3_constraint")


def figure_4_background(scan_results):
    """Figure 4: Sound horizon accumulation + r_s/D_A vs kappa."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))

    colors = {1.00: 'black', 0.99: '#1976D2', 0.97: '#FF9800', 0.95: '#D32F2F'}
    labels = {1.00: r'$\Lambda$CDM', 0.99: r'$\kappa_c = 0.99$',
              0.97: r'$\kappa_c = 0.97$', 0.95: r'$\kappa_c = 0.95$'}

    # Panel (a): r_s(z) accumulation
    for kc in [1.00, 0.99, 0.97, 0.95]:
        r = next((r for r in scan_results if abs(r['kappa_c'] - kc) < 0.001), None)
        if r and r.get('bg_file') and Path(r['bg_file']).exists():
            bg = load_background(r['bg_file'])
            mask = (bg['z'] >= 100) & (bg['z'] <= 1e5)
            ax1.plot(bg['z'][mask], bg['sound_horizon'][mask],
                     color=colors[kc], lw=2, label=labels[kc])
            # Annotate final r_s at z=1089
            rs_val = get_rs_at_z(bg, 1089)
            ax1.annotate(f'{rs_val:.1f} Mpc', xy=(1089, rs_val),
                         xytext=(400, rs_val + 3), fontsize=8, color=colors[kc],
                         arrowprops=dict(arrowstyle='->', color=colors[kc], lw=0.8))

    ax1.axvline(1089, color='gray', ls='--', lw=0.8, alpha=0.7)
    ax1.text(1150, 5, r'$z_{\rm rec}$', fontsize=10, color='gray')
    ax1.set_xscale('log')
    ax1.set_xlabel(r'Redshift $z$')
    ax1.set_ylabel(r'Comoving sound horizon $r_s(z)$ [Mpc]')
    ax1.set_title(r'(a) Sound horizon accumulation')
    ax1.set_xlim(100, 1e5)
    ax1.legend(loc='upper left', fontsize=9)
    ax1.invert_xaxis()

    # Panel (b): r_s and D_A at z=1089 vs kappa
    kappas_with_data = []
    rs_vals = []
    DA_vals = []
    theta_vals = []

    for r in sorted(scan_results, key=lambda x: x['kappa_c']):
        if r.get('rs') and r.get('DA'):
            kappas_with_data.append(r['kappa_c'])
            rs_vals.append(r['rs'])
            DA_vals.append(r['DA'])
            theta_vals.append(r['rs'] / r['DA'] * 1000)  # mrad

    ax2.plot(kappas_with_data, rs_vals, 'o-', color='#1976D2', lw=2,
             markersize=5, label=r'$r_s(z_*)$ [Mpc]')
    ax2b = ax2.twinx()
    ax2b.plot(kappas_with_data, DA_vals, 's--', color='#D32F2F', lw=2,
              markersize=5, label=r'$D_A(z_*)$ [Mpc]')

    ax2.set_xlabel(r'$\kappa_c$')
    ax2.set_ylabel(r'$r_s(z_*)$ [Mpc]', color='#1976D2')
    ax2b.set_ylabel(r'$D_A(z_*)$ [Mpc]', color='#D32F2F')
    ax2.set_title(r'(b) Acoustic scale components at recombination')

    # Combined legend
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2b.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=9)
    ax2.grid(True, alpha=0.2)

    plt.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(FIGURES_DIR / f'fig4_background.{ext}')
    plt.close(fig)
    print("  Saved fig4_background")


# ============================================================
# SECTION 4: Tables
# ============================================================

def table_1_bestfit(scan_results):
    """Table 1: Best-fit parameters for each kappa."""

    # Plain text
    header = (f"{'kappa_c':>8} {'h':>7} {'omega_cdm':>10} {'omega_b':>8} "
              f"{'n_s':>7} {'A_s':>10} {'chi2/dof':>10} {'H0':>6} {'r_s':>8}")
    lines = [header, "-" * 80]
    for r in sorted(scan_results, key=lambda x: -x['kappa_c']):
        rs_str = f"{r['rs']:.2f}" if r.get('rs') else 'N/A'
        lines.append(
            f"{r['kappa_c']:>8.3f} {r['h']:>7.4f} {r['omega_cdm']:>10.5f} "
            f"{r['omega_b']:>8.5f} {r['n_s']:>7.4f} {r['A_s']:>10.2e} "
            f"{r['chi2_dof']:>10.4f} {r['H0']:>6.1f} {rs_str:>8}"
        )
    txt = '\n'.join(lines)
    (FIGURES_DIR / 'table1_bestfit.txt').write_text(txt)

    # LaTeX
    latex_lines = [
        r'\begin{table}',
        r'\centering',
        r'\caption{Best-fit cosmological parameters for each $\kappa_c$ value.}',
        r'\label{tab:bestfit}',
        r'\begin{tabular}{ccccccccr}',
        r'\hline',
        r'$\kappa_c$ & $h$ & $\omega_{\rm cdm}$ & $\omega_b$ & $n_s$ & '
        r'$A_s$ & $\chi^2$/dof & $H_0$ & $r_s$ \\',
        r' & & & & & $[\times 10^{-9}]$ & & [km/s/Mpc] & [Mpc] \\',
        r'\hline',
    ]
    for r in sorted(scan_results, key=lambda x: -x['kappa_c']):
        rs_str = f"{r['rs']:.1f}" if r.get('rs') else '--'
        As_scaled = r['A_s'] * 1e9
        latex_lines.append(
            f"{r['kappa_c']:.3f} & {r['h']:.4f} & {r['omega_cdm']:.4f} & "
            f"{r['omega_b']:.4f} & {r['n_s']:.4f} & {As_scaled:.2f} & "
            f"{r['chi2_dof']:.3f} & {r['H0']:.1f} & {rs_str} \\\\"
        )
    latex_lines += [r'\hline', r'\end{tabular}', r'\end{table}']
    (FIGURES_DIR / 'table1_bestfit.tex').write_text('\n'.join(latex_lines))
    print("  Saved table1_bestfit")
    print(txt)


def table_2_chi2_breakdown(scan_results, p_ell, p_dl, p_sigma):
    """Table 2: chi2 breakdown by l-range for key models."""
    key_kappas = [1.00, 0.99, 0.97, 0.95]
    ranges = [(2, 29), (30, 800), (801, 2500), (2, 2500)]
    range_labels = ['l=2-29', 'l=30-800', 'l=801-2500', 'Total']

    header = f"{'Model':<25}"
    for rl in range_labels:
        header += f"  {rl:>14}"
    lines = [header, "-" * 85]

    for kc in key_kappas:
        r = next((r for r in scan_results if abs(r['kappa_c'] - kc) < 0.001), None)
        if r and r.get('ell') is not None:
            ell = np.array(r['ell']) if not isinstance(r['ell'], np.ndarray) else r['ell']
            dl = np.array(r['dl']) if not isinstance(r['dl'], np.ndarray) else r['dl']
            name = 'LCDM' if kc == 1.0 else f'GD kappa={kc}'
            line = f"{name:<25}"
            for l_min, l_max in ranges:
                chi2, ndof = compute_chi2(ell, dl, p_ell, p_dl, p_sigma, l_min, l_max)
                line += f"  {chi2:.0f}/{ndof:>4} = {chi2/ndof:.2f}"
            lines.append(line)

    txt = '\n'.join(lines)
    (FIGURES_DIR / 'table2_chi2.txt').write_text(txt)
    print("  Saved table2_chi2")
    print(txt)


# ============================================================
# SECTION 5: Paper numbers
# ============================================================

def generate_paper_numbers(scan_results, p_ell, p_dl, p_sigma):
    """Generate all numbers needed for paper text."""
    lines = []

    def add(s):
        lines.append(s)
        print(s)

    add("=" * 60)
    add("PAPER NUMBERS")
    add("=" * 60)

    lcdm = next((r for r in scan_results if abs(r['kappa_c'] - 1.0) < 0.001), None)
    gd99 = next((r for r in scan_results if abs(r['kappa_c'] - 0.99) < 0.001), None)
    gd97 = next((r for r in scan_results if abs(r['kappa_c'] - 0.97) < 0.001), None)
    gd96 = next((r for r in scan_results if abs(r['kappa_c'] - 0.96) < 0.001), None)
    gd95 = next((r for r in scan_results if abs(r['kappa_c'] - 0.95) < 0.001), None)

    add("\n--- LCDM Reference ---")
    if lcdm:
        add(f"  r_s(z=1089) = {lcdm['rs']:.2f} Mpc")
        add(f"  D_A(z=1089) = {lcdm['DA']:.2f} Mpc")
        add(f"  H0 = {lcdm['H0']:.1f} km/s/Mpc")
        add(f"  chi2 = {lcdm['chi2']:.1f} / {lcdm['ndof']} dof")
        add(f"  chi2/dof = {lcdm['chi2_dof']:.4f}")

    for label, r in [('kappa=0.99', gd99), ('kappa=0.97', gd97),
                     ('kappa=0.96', gd96), ('kappa=0.95', gd95)]:
        if r and lcdm:
            add(f"\n--- GD {label} ---")
            add(f"  H0 = {r['H0']:.1f} km/s/Mpc")
            add(f"  r_s = {r['rs']:.2f} Mpc")
            delta_rs = (r['rs'] - lcdm['rs']) / lcdm['rs'] * 100
            add(f"  Delta_r_s = {delta_rs:.2f}%")
            add(f"  chi2/dof = {r['chi2_dof']:.4f}")
            delta_chi2 = r['chi2'] - lcdm['chi2']
            add(f"  Delta_chi2 vs LCDM = {delta_chi2:.1f}")
            # Hubble tension
            sigma_lcdm = abs(H0_SHOES - lcdm['H0']) / np.sqrt(H0_SHOES_ERR**2 + 0.54**2)
            sigma_gd = abs(H0_SHOES - r['H0']) / np.sqrt(H0_SHOES_ERR**2 + 0.54**2)
            add(f"  Hubble tension: {sigma_lcdm:.1f}sigma (LCDM) -> {sigma_gd:.1f}sigma (GD)")

    # 95% CL constraint
    add("\n--- 95% CL Constraint ---")
    if lcdm:
        chi2_threshold = lcdm['chi2'] + 3.84
        add(f"  LCDM chi2 = {lcdm['chi2']:.1f}")
        add(f"  Threshold (Delta chi2 = 3.84) = {chi2_threshold:.1f}")

        # Find bound by interpolation
        sorted_results = sorted(scan_results, key=lambda x: x['kappa_c'], reverse=True)
        kappa_bound = None
        for i, r in enumerate(sorted_results):
            if r['chi2'] > chi2_threshold and i > 0:
                r_prev = sorted_results[i-1]
                # Linear interpolation
                f = (chi2_threshold - r_prev['chi2']) / (r['chi2'] - r_prev['chi2'])
                kappa_bound = r_prev['kappa_c'] + f * (r['kappa_c'] - r_prev['kappa_c'])
                break

        if kappa_bound:
            add(f"  95% CL lower bound: kappa > {kappa_bound:.3f}")
            # Find corresponding H0
            h0_at_bound = np.interp(kappa_bound,
                                     [r['kappa_c'] for r in scan_results],
                                     [r['H0'] for r in scan_results])
            add(f"  Maximum H0 at 95% CL: {h0_at_bound:.1f} km/s/Mpc")

    add("\n--- Phase C Justification ---")
    add("  omega_BD = 50,000 (Cassini-compatible)")
    add("  Perturbation correction = 1/omega_BD = 0.002%")
    add("  Background-only is correct to leading order in 1/omega_BD")

    txt = '\n'.join(lines)
    (FIGURES_DIR / 'paper_numbers.txt').write_text(txt)
    return txt


# ============================================================
# SECTION 6: Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='GD Paper Analysis Pipeline')
    parser.add_argument('--figures', action='store_true',
                        help='Regenerate figures from cached data only')
    parser.add_argument('--scan-only', action='store_true',
                        help='Run kappa scan and cache, skip figures')
    parser.add_argument('--clean', action='store_true',
                        help='Clear cache and re-run everything')
    args = parser.parse_args()

    setup_pub_style()
    FIGURES_DIR.mkdir(exist_ok=True)

    if args.clean:
        import shutil
        if CACHE_DIR.exists():
            shutil.rmtree(CACHE_DIR)
        print("Cache cleared.")

    # Load Planck data
    p_ell, p_dl, p_sigma = load_planck()
    print(f"Planck data: {len(p_ell)} multipoles (l={p_ell[0]:.0f}-{p_ell[-1]:.0f})")

    # Run or load kappa scan
    if not args.figures:
        scan_results = run_kappa_scan(KAPPA_VALUES, p_ell, p_dl, p_sigma)
        save_master_results(scan_results)
    else:
        scan_results = load_master_results()
        print(f"Loaded {len(scan_results)} cached results")

    if args.scan_only:
        print("\nScan complete. Run without --scan-only to generate figures.")
        return

    # Generate all deliverables
    print("\n=== GENERATING FIGURES ===")
    figure_1_theory(scan_results)
    figure_2_cmb(scan_results, p_ell, p_dl, p_sigma)
    figure_3_constraint(scan_results)
    figure_4_background(scan_results)

    print("\n=== GENERATING TABLES ===")
    table_1_bestfit(scan_results)
    print()
    table_2_chi2_breakdown(scan_results, p_ell, p_dl, p_sigma)

    print("\n=== GENERATING PAPER NUMBERS ===")
    generate_paper_numbers(scan_results, p_ell, p_dl, p_sigma)

    print(f"\nAll outputs saved to {FIGURES_DIR}/")
    print("Done.")


if __name__ == '__main__':
    main()
