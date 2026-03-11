#!/usr/bin/env python3
"""
GD-CLASS MCMC: Bayesian parameter estimation for the Glassy Dynamics model.

Uses emcee (affine-invariant ensemble sampler) to sample 7 parameters:
  h, omega_b, omega_cdm, n_s, ln(10^10 A_s), tau_reio, gd_kappa_c

Compares to Planck 2018 TT data (l >= 30) via Gaussian likelihood.
Runs CLASS via subprocess for fork-safe parallelization.

Usage:
  python3 gd_mcmc.py --nwalkers 16 --nsteps 5 --workers 2    # smoke test
  python3 gd_mcmc.py --nwalkers 16 --nsteps 3000 --workers 8  # production
  python3 gd_mcmc.py --resume --nsteps 5000 --workers 8       # resume
  python3 gd_mcmc.py --plots-only --burn-in 500               # plots from chain
"""

import argparse
import logging
import os
import sys
import subprocess
import time
import tempfile
import uuid
from pathlib import Path
from multiprocessing import Pool

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Logging setup
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger('gd_mcmc')
import emcee
import h5py

GD_DIR = Path(__file__).parent
CLASS_BIN = GD_DIR / 'class'
PLANCK_FILE = GD_DIR / 'data' / 'planck_2018_TT.txt'
OUTPUT_DIR = GD_DIR / 'output'
CHAIN_FILE = OUTPUT_DIR / 'gd_mcmc_chains.h5'
T_CMB = 2.7255
TCMB2 = (T_CMB * 1e6)**2

# Parameter names and labels
PARAM_NAMES = ['h', 'omega_b', 'omega_cdm', 'n_s', 'ln10As', 'tau_reio', 'gd_kappa_c']
PARAM_LABELS = [r'$h$', r'$\omega_b$', r'$\omega_{cdm}$', r'$n_s$',
                r'$\ln(10^{10}A_s)$', r'$\tau_{reio}$', r'$\kappa_c$']

# Flat prior bounds
PRIOR_BOUNDS = {
    'h':          (0.60, 0.80),
    'omega_b':    (0.019, 0.025),
    'omega_cdm':  (0.08, 0.18),
    'n_s':        (0.90, 1.05),
    'ln10As':     (2.5, 3.5),
    'tau_reio':   (0.02, 0.12),
    'gd_kappa_c': (0.90, 1.00),
}

# Gaussian prior on tau (Planck low-ell EE constraint)
TAU_PRIOR_MEAN = 0.054
TAU_PRIOR_SIGMA = 0.007

# Fixed GD parameters
GD_FIXED = {
    'z_freeze': 1040,
    'z_onset': 1140,
    'beta': 0.5,
}

# Starting point (near grid-scan best-fit for kappa=0.98)
P0_CENTER = [0.71, 0.0220, 0.12, 0.97, 3.044, 0.054, 0.98]

# Planck data (loaded once, shared via module-level)
_planck_ell = None
_planck_dl = None
_planck_sigma = None


def load_planck():
    """Load Planck 2018 TT data. Returns ell, D_l, sigma."""
    data = np.loadtxt(PLANCK_FILE)
    ell = data[:, 0]
    dl = data[:, 1]
    sigma = (data[:, 2] + data[:, 3]) / 2.0
    return ell, dl, sigma


def ensure_planck_loaded():
    """Load Planck data into module-level variables if not already loaded."""
    global _planck_ell, _planck_dl, _planck_sigma
    if _planck_ell is None:
        _planck_ell, _planck_dl, _planck_sigma = load_planck()


def make_ini_content(h, omega_b, omega_cdm, A_s, n_s, tau_reio,
                     kappa_c, z_freeze, z_onset, beta, root):
    """Generate .ini content for a CLASS run."""
    return f"""root = {root}
write_background = no
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


def run_class_subprocess(params):
    """Run CLASS via subprocess. Returns (ell, dl) or None on failure.

    Uses unique temp file names to avoid collisions in parallel.
    """
    uid = f"{os.getpid()}_{uuid.uuid4().hex[:8]}"
    root = f'output/_mcmc_{uid}_'
    ini_content = make_ini_content(root=root, **params)
    ini_path = GD_DIR / f'_mcmc_{uid}.ini'

    try:
        with open(ini_path, 'w') as f:
            f.write(ini_content)

        result = subprocess.run(
            [str(CLASS_BIN), str(ini_path)],
            cwd=str(GD_DIR),
            capture_output=True, text=True, timeout=120
        )
        if result.returncode != 0:
            return None

        cl_file = GD_DIR / f'{root}cl.dat'
        if not cl_file.exists():
            # Try with numbering
            for i in range(5):
                cf = GD_DIR / f'{root}{i:02d}_cl.dat'
                if cf.exists():
                    cl_file = cf
                    break

        if not cl_file.exists():
            return None

        data = np.loadtxt(cl_file)
        ell = data[:, 0]
        cl_tt = data[:, 1]  # dimensionless l(l+1)/2pi C_l
        dl_tt = cl_tt * TCMB2
        return ell, dl_tt

    except (subprocess.TimeoutExpired, Exception):
        return None

    finally:
        # Clean up temp files
        ini_path.unlink(missing_ok=True)
        for f in (GD_DIR / 'output').glob(f'_mcmc_{uid}_*'):
            f.unlink(missing_ok=True)


def compute_chi2(model_ell, model_dl, p_ell, p_dl, p_sigma, l_min=30, l_max=2500):
    """Compute chi-squared between model and Planck."""
    mask = (p_ell >= l_min) & (p_ell <= l_max)
    model_interp = np.interp(p_ell[mask], model_ell, model_dl)
    chi2 = np.sum(((model_interp - p_dl[mask]) / p_sigma[mask]) ** 2)
    return chi2


def log_prior(theta):
    """Log prior: flat bounds on all params + Gaussian on tau."""
    h, omega_b, omega_cdm, n_s, ln10As, tau_reio, gd_kappa_c = theta

    # Flat bounds
    bounds = [PRIOR_BOUNDS[name] for name in PARAM_NAMES]
    for val, (lo, hi) in zip(theta, bounds):
        if not (lo < val < hi):
            return -np.inf

    # Gaussian prior on tau (Planck low-ell EE)
    lp = -0.5 * ((tau_reio - TAU_PRIOR_MEAN) / TAU_PRIOR_SIGMA) ** 2
    return lp


_eval_count = 0
_eval_start_time = None

def log_likelihood(theta):
    """Log likelihood: run CLASS, compare to Planck TT."""
    global _eval_count, _eval_start_time
    ensure_planck_loaded()

    h, omega_b, omega_cdm, n_s, ln10As, tau_reio, gd_kappa_c = theta
    A_s = np.exp(ln10As) * 1e-10

    params = {
        'h': h, 'omega_b': omega_b, 'omega_cdm': omega_cdm,
        'n_s': n_s, 'A_s': A_s, 'tau_reio': tau_reio,
        'kappa_c': gd_kappa_c, **GD_FIXED,
    }

    t0 = time.time()
    result = run_class_subprocess(params)
    dt = time.time() - t0

    if result is None:
        log.warning(f"CLASS FAILED | h={h:.4f} kc={gd_kappa_c:.4f} n_s={n_s:.4f} ({dt:.1f}s)")
        return -np.inf

    model_ell, model_dl = result
    chi2 = compute_chi2(model_ell, model_dl, _planck_ell, _planck_dl, _planck_sigma)
    ndof = np.sum((_planck_ell >= 30) & (_planck_ell <= 2500))
    chi2_red = chi2 / ndof

    _eval_count += 1
    if _eval_start_time is None:
        _eval_start_time = time.time()
    elapsed = time.time() - _eval_start_time
    rate = _eval_count / elapsed if elapsed > 0 else 0

    if _eval_count % 16 == 0:  # log every 16 evals (= 1 step for 16 walkers)
        log.info(f"eval #{_eval_count:>5d} | chi2/dof={chi2_red:.3f} H0={h*100:.1f} kc={gd_kappa_c:.3f} "
                 f"n_s={n_s:.4f} | {dt:.1f}s | {rate:.1f} eval/min")

    return -0.5 * chi2


def log_probability(theta):
    """Log posterior = log prior + log likelihood."""
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(theta)
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll


def gelman_rubin(chains):
    """Compute Gelman-Rubin R-hat statistic.

    Args:
        chains: array of shape (nwalkers, nsteps, ndim)
    Returns:
        R-hat for each parameter (ndim,)
    """
    m, n, d = chains.shape
    chain_means = np.mean(chains, axis=1)
    grand_mean = np.mean(chain_means, axis=0)
    B = n / (m - 1) * np.sum((chain_means - grand_mean) ** 2, axis=0)
    W = np.mean(np.var(chains, axis=1, ddof=1), axis=0)
    var_hat = (1 - 1.0 / n) * W + B / n
    R_hat = np.sqrt(var_hat / W)
    return R_hat


def run_mcmc(nwalkers, nsteps, workers, resume):
    """Run emcee MCMC sampler."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    ndim = len(PARAM_NAMES)

    # HDF5 backend for crash-safe saving
    backend = emcee.backends.HDFBackend(str(CHAIN_FILE))

    if resume and CHAIN_FILE.exists():
        log.info(f"Resuming from {CHAIN_FILE} ({backend.iteration} steps already done)")
        # Get last walker positions from saved chain
        p0 = backend.get_last_sample()
    else:
        backend.reset(nwalkers, ndim)
        # Initialize walkers in a ball around the best-fit
        p0 = np.array(P0_CENTER) + 1e-3 * np.random.randn(nwalkers, ndim)
        # Clip to prior bounds
        for i, name in enumerate(PARAM_NAMES):
            lo, hi = PRIOR_BOUNDS[name]
            p0[:, i] = np.clip(p0[:, i], lo + 1e-6, hi - 1e-6)
        log.info(f"Starting fresh: {nwalkers} walkers, {nsteps} steps, {workers} workers")
        log.info(f"Initial center: h={P0_CENTER[0]}, kc={P0_CENTER[6]}, n_s={P0_CENTER[3]}")

    log.info("Starting MCMC sampling...")
    step_start = time.time()

    if workers > 1:
        with Pool(workers) as pool:
            sampler = emcee.EnsembleSampler(
                nwalkers, ndim, log_probability, pool=pool, backend=backend
            )
            for step, state in enumerate(sampler.sample(p0, iterations=nsteps, progress=True)):
                if (step + 1) % 10 == 0:
                    elapsed = time.time() - step_start
                    rate = (step + 1) / elapsed * 3600
                    best_logp = np.max(state.log_prob)
                    median_params = np.median(state.coords, axis=0)
                    accept = np.mean(sampler.acceptance_fraction)
                    log.info(
                        f"Step {step+1:>4d}/{nsteps} | "
                        f"best logP={best_logp:.1f} | "
                        f"H0={median_params[0]*100:.1f} kc={median_params[6]:.3f} n_s={median_params[3]:.4f} | "
                        f"accept={accept:.2f} | "
                        f"{rate:.0f} steps/hr | "
                        f"ETA {(nsteps-step-1)/rate*60:.0f} min"
                    )
    else:
        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, log_probability, backend=backend
        )
        for step, state in enumerate(sampler.sample(p0, iterations=nsteps, progress=True)):
            if (step + 1) % 10 == 0:
                elapsed = time.time() - step_start
                rate = (step + 1) / elapsed * 3600
                best_logp = np.max(state.log_prob)
                median_params = np.median(state.coords, axis=0)
                accept = np.mean(sampler.acceptance_fraction)
                log.info(
                    f"Step {step+1:>4d}/{nsteps} | "
                    f"best logP={best_logp:.1f} | "
                    f"H0={median_params[0]*100:.1f} kc={median_params[6]:.3f} n_s={median_params[3]:.4f} | "
                    f"accept={accept:.2f} | "
                    f"{rate:.0f} steps/hr | "
                    f"ETA {(nsteps-step-1)/rate*60:.0f} min"
                )

    print(f"\nDone. Total steps: {backend.iteration}")
    print(f"Chain shape: {backend.get_chain().shape}")
    print(f"Acceptance fraction: {np.mean(sampler.acceptance_fraction):.3f}")

    try:
        tau = sampler.get_autocorr_time(quiet=True)
        print(f"Autocorrelation times: {tau}")
    except emcee.autocorr.AutocorrError:
        print("Warning: chain too short for reliable autocorrelation estimate")

    return sampler


def plot_corner(burn_in, thin=1):
    """Generate corner plot from saved chains."""
    backend = emcee.backends.HDFBackend(str(CHAIN_FILE))
    chain = backend.get_chain(discard=burn_in, thin=thin, flat=True)
    print(f"Flat chain shape after burn-in={burn_in}, thin={thin}: {chain.shape}")

    try:
        import corner as corner_pkg
    except ImportError:
        print("ERROR: corner package not installed. pip install corner")
        return

    fig = corner_pkg.corner(
        chain, labels=PARAM_LABELS,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_kwargs={"fontsize": 11},
        label_kwargs={"fontsize": 12},
    )
    fig.suptitle('GD-CLASS MCMC: 7-parameter posterior', fontsize=14, y=1.02)
    outfile = OUTPUT_DIR / 'gd_mcmc_corner.png'
    fig.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {outfile}")

    # Also make a corner plot with derived H0
    h_samples = chain[:, 0] * 100  # H0 = h * 100
    chain_with_h0 = np.column_stack([chain, h_samples])
    labels_h0 = PARAM_LABELS + [r'$H_0$']

    fig2 = corner_pkg.corner(
        chain_with_h0, labels=labels_h0,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_kwargs={"fontsize": 10},
        label_kwargs={"fontsize": 11},
    )
    fig2.suptitle('GD-CLASS MCMC: posterior with derived H0', fontsize=14, y=1.02)
    outfile2 = OUTPUT_DIR / 'gd_mcmc_corner_h0.png'
    fig2.savefig(outfile2, dpi=150, bbox_inches='tight')
    plt.close(fig2)
    print(f"Saved: {outfile2}")


def plot_traces(burn_in):
    """Generate trace plots for convergence inspection."""
    backend = emcee.backends.HDFBackend(str(CHAIN_FILE))
    chain = backend.get_chain()  # (nsteps, nwalkers, ndim)
    nsteps, nwalkers, ndim = chain.shape

    fig, axes = plt.subplots(ndim, 1, figsize=(12, 2.5 * ndim), sharex=True)
    for i in range(ndim):
        ax = axes[i]
        for w in range(nwalkers):
            ax.plot(chain[:, w, i], alpha=0.3, lw=0.5)
        ax.set_ylabel(PARAM_LABELS[i], fontsize=11)
        ax.axvline(burn_in, color='red', ls='--', lw=1, label='burn-in' if i == 0 else None)
    axes[0].legend(loc='upper right')
    axes[-1].set_xlabel('Step')
    fig.suptitle('MCMC Trace Plots', fontsize=14)
    plt.tight_layout()

    outfile = OUTPUT_DIR / 'gd_mcmc_traces.png'
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved: {outfile}")


def plot_bestfit(burn_in, thin=1):
    """Plot best-fit spectrum vs Planck."""
    ensure_planck_loaded()

    backend = emcee.backends.HDFBackend(str(CHAIN_FILE))
    chain = backend.get_chain(discard=burn_in, thin=thin, flat=True)
    log_probs = backend.get_log_prob(discard=burn_in, thin=thin, flat=True)

    # Find MAP (maximum a posteriori)
    best_idx = np.argmax(log_probs)
    best_theta = chain[best_idx]
    print(f"MAP parameters: {dict(zip(PARAM_NAMES, best_theta))}")

    # Run CLASS for best-fit
    h, omega_b, omega_cdm, n_s, ln10As, tau_reio, gd_kappa_c = best_theta
    A_s = np.exp(ln10As) * 1e-10
    params_gd = {
        'h': h, 'omega_b': omega_b, 'omega_cdm': omega_cdm,
        'n_s': n_s, 'A_s': A_s, 'tau_reio': tau_reio,
        'kappa_c': gd_kappa_c, **GD_FIXED,
    }
    gd_result = run_class_subprocess(params_gd)

    # Run LCDM baseline
    params_lcdm = {
        'h': 0.6736, 'omega_b': 0.02237, 'omega_cdm': 0.1200,
        'n_s': 0.9649, 'A_s': 2.1e-9, 'tau_reio': 0.0544,
        'kappa_c': 1.0, **GD_FIXED,
    }
    lcdm_result = run_class_subprocess(params_lcdm)

    if gd_result is None:
        print("ERROR: best-fit CLASS run failed")
        return

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), gridspec_kw={'height_ratios': [3, 1]})

    # Panel 1: TT spectrum
    ax = axes[0]
    ax.errorbar(_planck_ell, _planck_dl, yerr=_planck_sigma, fmt='.', color='gray',
                alpha=0.3, markersize=2, elinewidth=0.5, label='Planck 2018 TT', zorder=1)
    if lcdm_result:
        ax.plot(lcdm_result[0], lcdm_result[1], 'k-', lw=1.5, label=r'$\Lambda$CDM', zorder=2)
    ax.plot(gd_result[0], gd_result[1], 'r-', lw=1.5,
            label=rf'GD MAP ($\kappa_c$={gd_kappa_c:.3f}, h={h:.3f})', zorder=3)
    ax.set_xlim(2, 2500)
    ax.set_xscale('log')
    ax.set_xlabel(r'Multipole $\ell$')
    ax.set_ylabel(r'$\mathcal{D}_\ell^{TT}$ [$\mu K^2$]')
    ax.set_title('MCMC Best-fit vs Planck 2018 TT')
    ax.legend(loc='upper right')

    # Panel 2: Residuals
    ax = axes[1]
    gd_interp = np.interp(_planck_ell, gd_result[0], gd_result[1])
    gd_resid = (gd_interp - _planck_dl) / _planck_sigma
    mask = _planck_ell >= 2
    bin_size = 20
    n_bins = len(_planck_ell[mask]) // bin_size
    ell_binned, resid_binned = [], []
    for b in range(n_bins):
        s, e = b * bin_size, (b + 1) * bin_size
        ell_binned.append(np.mean(_planck_ell[mask][s:e]))
        resid_binned.append(np.mean(gd_resid[mask][s:e]))
    ax.plot(ell_binned, resid_binned, 'ro-', markersize=3, lw=1, alpha=0.7, label='GD MAP')

    if lcdm_result:
        lcdm_interp = np.interp(_planck_ell, lcdm_result[0], lcdm_result[1])
        lcdm_resid = (lcdm_interp - _planck_dl) / _planck_sigma
        ell_b2, resid_b2 = [], []
        for b in range(n_bins):
            s, e = b * bin_size, (b + 1) * bin_size
            ell_b2.append(np.mean(_planck_ell[mask][s:e]))
            resid_b2.append(np.mean(lcdm_resid[mask][s:e]))
        ax.plot(ell_b2, resid_b2, 'ko-', markersize=3, lw=1, alpha=0.7, label=r'$\Lambda$CDM')

    ax.axhline(0, color='gray', ls=':', lw=0.5)
    ax.set_xlim(2, 2500)
    ax.set_xscale('log')
    ax.set_xlabel(r'Multipole $\ell$')
    ax.set_ylabel(r'$(D_\ell^{model} - D_\ell^{Planck}) / \sigma$')
    ax.legend(loc='upper right')
    plt.tight_layout()

    outfile = OUTPUT_DIR / 'gd_mcmc_bestfit.png'
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved: {outfile}")


def print_summary(burn_in, thin=1):
    """Print parameter summary table with uncertainties."""
    backend = emcee.backends.HDFBackend(str(CHAIN_FILE))
    chain = backend.get_chain(discard=burn_in, thin=thin, flat=True)
    log_probs = backend.get_log_prob(discard=burn_in, thin=thin, flat=True)

    # Convergence
    chain_3d = backend.get_chain(discard=burn_in, thin=thin)  # (nsteps, nwalkers, ndim)
    # Transpose to (nwalkers, nsteps, ndim) for gelman_rubin
    R_hat = gelman_rubin(np.transpose(chain_3d, (1, 0, 2)))

    lines = []
    lines.append("=" * 80)
    lines.append("GD-CLASS MCMC PARAMETER SUMMARY")
    lines.append("=" * 80)
    lines.append(f"Total samples (post burn-in): {len(chain)}")
    lines.append(f"Best log-probability: {np.max(log_probs):.2f}")
    lines.append("")
    lines.append(f"{'Parameter':<15} {'Median':>10} {'68% CI':>20} {'95% CI':>20} {'R-hat':>8}")
    lines.append("-" * 80)

    for i, name in enumerate(PARAM_NAMES):
        samples = chain[:, i]
        median = np.median(samples)
        q16, q84 = np.percentile(samples, [16, 84])
        q2, q98 = np.percentile(samples, [2.5, 97.5])
        lines.append(f"{name:<15} {median:>10.5f} [{q16:.5f}, {q84:.5f}] [{q2:.5f}, {q98:.5f}] {R_hat[i]:>8.4f}")

    # Derived parameters
    lines.append("-" * 80)
    h_samples = chain[:, 0]
    H0_samples = h_samples * 100
    H0_med = np.median(H0_samples)
    H0_16, H0_84 = np.percentile(H0_samples, [16, 84])
    H0_2, H0_98 = np.percentile(H0_samples, [2.5, 97.5])
    lines.append(f"{'H0 (derived)':<15} {H0_med:>10.2f} [{H0_16:.2f}, {H0_84:.2f}] [{H0_2:.2f}, {H0_98:.2f}]")

    lines.append("")
    lines.append("Convergence: R-hat < 1.1 for all parameters = " +
                 ("PASS" if np.all(R_hat < 1.1) else "FAIL (need more steps)"))
    lines.append("=" * 80)

    summary = "\n".join(lines)
    print(summary)

    outfile = OUTPUT_DIR / 'gd_mcmc_summary.txt'
    with open(outfile, 'w') as f:
        f.write(summary + "\n")
    print(f"Saved: {outfile}")


def main():
    parser = argparse.ArgumentParser(description='GD-CLASS MCMC parameter estimation')
    parser.add_argument('--nwalkers', type=int, default=16, help='Number of walkers (default: 16)')
    parser.add_argument('--nsteps', type=int, default=2000, help='Number of steps (default: 2000)')
    parser.add_argument('--workers', type=int, default=1, help='Number of parallel workers (default: 1)')
    parser.add_argument('--burn-in', type=int, default=500, help='Burn-in steps to discard (default: 500)')
    parser.add_argument('--thin', type=int, default=1, help='Thinning factor (default: 1)')
    parser.add_argument('--resume', action='store_true', help='Resume from existing chain file')
    parser.add_argument('--plots-only', action='store_true', help='Skip MCMC, generate plots from chain')
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(exist_ok=True)

    if not args.plots_only:
        if not CLASS_BIN.exists():
            print(f"ERROR: CLASS binary not found at {CLASS_BIN}")
            print("Run 'make clean && make -j4' first")
            sys.exit(1)

        if not PLANCK_FILE.exists():
            print(f"ERROR: Planck data not found at {PLANCK_FILE}")
            sys.exit(1)

        print("=" * 70)
        print("GD-CLASS MCMC Parameter Estimation")
        print("=" * 70)
        print(f"Parameters: {PARAM_NAMES}")
        print(f"Walkers: {args.nwalkers}, Steps: {args.nsteps}, Workers: {args.workers}")
        print(f"Chain file: {CHAIN_FILE}")
        print()

        run_mcmc(args.nwalkers, args.nsteps, args.workers, args.resume)

    if not CHAIN_FILE.exists():
        print(f"ERROR: No chain file found at {CHAIN_FILE}")
        sys.exit(1)

    # Adjust burn-in if chain is too short
    backend = emcee.backends.HDFBackend(str(CHAIN_FILE))
    total_steps = backend.iteration
    burn_in = min(args.burn_in, max(0, total_steps - 2))
    if burn_in != args.burn_in:
        print(f"Warning: adjusted burn-in from {args.burn_in} to {burn_in} (chain has {total_steps} steps)")

    print("\nGenerating diagnostics and plots...")
    print_summary(burn_in, args.thin)
    plot_traces(burn_in)
    plot_corner(burn_in, args.thin)
    plot_bestfit(burn_in, args.thin)
    print("\nAll done.")


if __name__ == '__main__':
    main()
