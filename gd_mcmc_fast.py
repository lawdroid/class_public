#!/usr/bin/env python3
"""
Fast MCMC with Planck n_s prior — designed for ~30 min runs with immediate feedback.

Changes from gd_mcmc.py:
  - Gaussian prior on n_s (0.9649 ± 0.0042) from Planck TT+TE+EE+lensing
  - Initialize walkers at grid-scan best-fit (h=0.674, not 0.71)
  - Print every step (not every 10)
  - Default: 100 steps, 4 workers

Usage:
  python3 gd_mcmc_fast.py                        # 100 steps, ~25 min
  python3 gd_mcmc_fast.py --nsteps 200            # 200 steps, ~50 min
  python3 gd_mcmc_fast.py --plots-only --burn-in 20  # plots from chain
"""

import argparse
import logging
import os
import sys
import subprocess
import time
import uuid
from pathlib import Path
from multiprocessing import Pool

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger('gd_mcmc_fast')
import emcee
import h5py

GD_DIR = Path(__file__).parent
CLASS_BIN = GD_DIR / 'class'
PLANCK_FILE = GD_DIR / 'data' / 'planck_2018_TT.txt'
OUTPUT_DIR = GD_DIR / 'output'
CHAIN_FILE = OUTPUT_DIR / 'gd_mcmc_fast_chains.h5'
T_CMB = 2.7255
TCMB2 = (T_CMB * 1e6)**2

PARAM_NAMES = ['h', 'omega_b', 'omega_cdm', 'n_s', 'ln10As', 'tau_reio', 'gd_kappa_c']
PARAM_LABELS = [r'$h$', r'$\omega_b$', r'$\omega_{cdm}$', r'$n_s$',
                r'$\ln(10^{10}A_s)$', r'$\tau_{reio}$', r'$\kappa_c$']

# Prior bounds
PRIOR_BOUNDS = {
    'h':          (0.60, 0.80),
    'omega_b':    (0.019, 0.025),
    'omega_cdm':  (0.08, 0.18),
    'n_s':        (0.90, 1.05),
    'ln10As':     (2.5, 3.5),
    'tau_reio':   (0.02, 0.12),
    'gd_kappa_c': (0.90, 1.00),
}

# Gaussian priors
TAU_PRIOR_MEAN = 0.054
TAU_PRIOR_SIGMA = 0.007
NS_PRIOR_MEAN = 0.9649    # Planck 2018 TT+TE+EE+lensing
NS_PRIOR_SIGMA = 0.0042

# Fixed GD parameters
GD_FIXED = {
    'z_freeze': 1040,
    'z_onset': 1140,
    'beta': 0.5,
}

# Start near LCDM best-fit (grid scan showed this is where the posterior is)
P0_CENTER = [0.6736, 0.02237, 0.1200, 0.9649, 3.044, 0.054, 0.995]

_planck_ell = None
_planck_dl = None
_planck_sigma = None


def load_planck():
    data = np.loadtxt(PLANCK_FILE)
    return data[:, 0], data[:, 1], (data[:, 2] + data[:, 3]) / 2.0


def ensure_planck_loaded():
    global _planck_ell, _planck_dl, _planck_sigma
    if _planck_ell is None:
        _planck_ell, _planck_dl, _planck_sigma = load_planck()


def make_ini_content(h, omega_b, omega_cdm, A_s, n_s, tau_reio,
                     kappa_c, z_freeze, z_onset, beta, root):
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
    uid = f"{os.getpid()}_{uuid.uuid4().hex[:8]}"
    root = f'output/_mcmc_fast_{uid}_'
    ini_content = make_ini_content(root=root, **params)
    ini_path = GD_DIR / f'_mcmc_fast_{uid}.ini'

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
            for i in range(5):
                cf = GD_DIR / f'{root}{i:02d}_cl.dat'
                if cf.exists():
                    cl_file = cf
                    break

        if not cl_file.exists():
            return None

        data = np.loadtxt(cl_file)
        ell = data[:, 0]
        cl_tt = data[:, 1]
        dl_tt = cl_tt * TCMB2
        return ell, dl_tt

    except (subprocess.TimeoutExpired, Exception):
        return None

    finally:
        ini_path.unlink(missing_ok=True)
        for f in (GD_DIR / 'output').glob(f'_mcmc_fast_{uid}_*'):
            f.unlink(missing_ok=True)


def compute_chi2(model_ell, model_dl, p_ell, p_dl, p_sigma, l_min=30, l_max=2500):
    mask = (p_ell >= l_min) & (p_ell <= l_max)
    model_interp = np.interp(p_ell[mask], model_ell, model_dl)
    return np.sum(((model_interp - p_dl[mask]) / p_sigma[mask]) ** 2)


def log_prior(theta):
    h, omega_b, omega_cdm, n_s, ln10As, tau_reio, gd_kappa_c = theta

    bounds = [PRIOR_BOUNDS[name] for name in PARAM_NAMES]
    for val, (lo, hi) in zip(theta, bounds):
        if not (lo < val < hi):
            return -np.inf

    # Gaussian prior on tau (Planck low-ell EE)
    lp = -0.5 * ((tau_reio - TAU_PRIOR_MEAN) / TAU_PRIOR_SIGMA) ** 2
    # Gaussian prior on n_s (Planck TT+TE+EE+lensing)
    lp += -0.5 * ((n_s - NS_PRIOR_MEAN) / NS_PRIOR_SIGMA) ** 2
    return lp


_eval_count = 0
_eval_start_time = None

def log_likelihood(theta):
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
        return -np.inf

    model_ell, model_dl = result
    chi2 = compute_chi2(model_ell, model_dl, _planck_ell, _planck_dl, _planck_sigma)
    ndof = np.sum((_planck_ell >= 30) & (_planck_ell <= 2500))
    chi2_red = chi2 / ndof

    _eval_count += 1
    if _eval_start_time is None:
        _eval_start_time = time.time()

    return -0.5 * chi2


def log_probability(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(theta)
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll


def gelman_rubin(chains):
    m, n, d = chains.shape
    chain_means = np.mean(chains, axis=1)
    grand_mean = np.mean(chain_means, axis=0)
    B = n / (m - 1) * np.sum((chain_means - grand_mean) ** 2, axis=0)
    W = np.mean(np.var(chains, axis=1, ddof=1), axis=0)
    var_hat = (1 - 1.0 / n) * W + B / n
    R_hat = np.sqrt(var_hat / W)
    return R_hat


def run_mcmc(nwalkers, nsteps, workers, resume):
    OUTPUT_DIR.mkdir(exist_ok=True)
    ndim = len(PARAM_NAMES)

    backend = emcee.backends.HDFBackend(str(CHAIN_FILE))

    if resume and CHAIN_FILE.exists():
        log.info(f"Resuming from {CHAIN_FILE} ({backend.iteration} steps)")
        p0 = backend.get_last_sample()
    else:
        backend.reset(nwalkers, ndim)
        p0 = np.array(P0_CENTER) + 1e-3 * np.random.randn(nwalkers, ndim)
        for i, name in enumerate(PARAM_NAMES):
            lo, hi = PRIOR_BOUNDS[name]
            p0[:, i] = np.clip(p0[:, i], lo + 1e-6, hi - 1e-6)
        log.info(f"Starting: {nwalkers} walkers, {nsteps} steps, {workers} workers")
        log.info(f"Center: h={P0_CENTER[0]}, kc={P0_CENTER[6]}, n_s={P0_CENTER[3]}")
        log.info(f"Priors: n_s = {NS_PRIOR_MEAN} ± {NS_PRIOR_SIGMA}, tau = {TAU_PRIOR_MEAN} ± {TAU_PRIOR_SIGMA}")

    step_start = time.time()

    pool = Pool(workers) if workers > 1 else None
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_probability,
        pool=pool, backend=backend
    )

    try:
        for step, state in enumerate(sampler.sample(p0, iterations=nsteps, progress=False)):
            elapsed = time.time() - step_start
            rate = (step + 1) / elapsed * 3600 if elapsed > 0 else 0
            best_logp = np.max(state.log_prob)
            med = np.median(state.coords, axis=0)
            accept = np.mean(sampler.acceptance_fraction)

            # Print EVERY step for immediate feedback
            chi2_approx = -2 * best_logp / 2471
            print(
                f"Step {step+1:>4d}/{nsteps} | "
                f"H0={med[0]*100:.1f} kc={med[6]:.4f} ns={med[3]:.4f} | "
                f"chi2/dof≈{chi2_approx:.3f} | "
                f"accept={accept:.2f} | "
                f"{elapsed:.0f}s | "
                f"{rate:.0f} steps/hr",
                flush=True
            )
    finally:
        if pool:
            pool.close()
            pool.join()

    print(f"\nDone. {backend.iteration} steps. Accept: {np.mean(sampler.acceptance_fraction):.3f}")
    try:
        tau = sampler.get_autocorr_time(quiet=True)
        print(f"Autocorrelation: {tau}")
    except:
        print("Chain too short for autocorrelation estimate (expected for fast run)")

    return sampler


def print_summary(burn_in, thin=1):
    backend = emcee.backends.HDFBackend(str(CHAIN_FILE))
    chain = backend.get_chain(discard=burn_in, thin=thin, flat=True)
    log_probs = backend.get_log_prob(discard=burn_in, thin=thin, flat=True)
    chain_3d = backend.get_chain(discard=burn_in, thin=thin)

    if chain_3d.shape[0] < 3:
        print("Too few steps for summary after burn-in")
        return

    R_hat = gelman_rubin(np.transpose(chain_3d, (1, 0, 2)))

    print("=" * 70)
    print("GD-CLASS FAST MCMC — WITH n_s PRIOR (0.9649 ± 0.0042)")
    print("=" * 70)
    print(f"Samples (post burn-in): {len(chain)}")
    print(f"Best log-probability: {np.max(log_probs):.2f}")
    print()
    print(f"{'Parameter':<15} {'Median':>10} {'68% CI':>20} {'95% CI':>20} {'R-hat':>8}")
    print("-" * 70)

    for i, name in enumerate(PARAM_NAMES):
        s = chain[:, i]
        med = np.median(s)
        q16, q84 = np.percentile(s, [16, 84])
        q025, q975 = np.percentile(s, [2.5, 97.5])
        print(f"{name:<15} {med:>10.5f} [{q16:.5f}, {q84:.5f}] [{q025:.5f}, {q975:.5f}] {R_hat[i]:>8.4f}")

    print("-" * 70)
    H0 = chain[:, 0] * 100
    med = np.median(H0)
    q16, q84 = np.percentile(H0, [16, 84])
    q025, q975 = np.percentile(H0, [2.5, 97.5])
    print(f"{'H0 (derived)':<15} {med:>10.2f} [{q16:.2f}, {q84:.2f}] [{q025:.2f}, {q975:.2f}]")
    print()
    print("R-hat < 1.1: " + ("PASS" if np.all(R_hat < 1.1) else "FAIL (need more steps)"))
    print("=" * 70)

    outfile = OUTPUT_DIR / 'gd_mcmc_fast_summary.txt'
    # Also save to file
    with open(outfile, 'w') as f:
        f.write(f"Samples: {len(chain)}, Best logP: {np.max(log_probs):.2f}\n")
        for i, name in enumerate(PARAM_NAMES):
            s = chain[:, i]
            f.write(f"{name}: {np.median(s):.5f} [{np.percentile(s,2.5):.5f}, {np.percentile(s,97.5):.5f}]\n")
        f.write(f"H0: {np.median(H0):.2f} [{np.percentile(H0,2.5):.2f}, {np.percentile(H0,97.5):.2f}]\n")
    print(f"Saved: {outfile}")


def main():
    parser = argparse.ArgumentParser(description='Fast GD-CLASS MCMC with n_s prior')
    parser.add_argument('--nwalkers', type=int, default=16)
    parser.add_argument('--nsteps', type=int, default=100)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--burn-in', type=int, default=20)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--plots-only', action='store_true')
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(exist_ok=True)

    if not args.plots_only:
        if not CLASS_BIN.exists():
            print(f"ERROR: CLASS binary not found at {CLASS_BIN}")
            sys.exit(1)

        print("=" * 70)
        print("GD-CLASS FAST MCMC — with Planck n_s prior")
        print(f"  n_s = {NS_PRIOR_MEAN} ± {NS_PRIOR_SIGMA}")
        print(f"  tau = {TAU_PRIOR_MEAN} ± {TAU_PRIOR_SIGMA}")
        print(f"  {args.nwalkers} walkers, {args.nsteps} steps, {args.workers} workers")
        print("=" * 70)

        run_mcmc(args.nwalkers, args.nsteps, args.workers, args.resume)

    if CHAIN_FILE.exists():
        backend = emcee.backends.HDFBackend(str(CHAIN_FILE))
        burn_in = min(args.burn_in, max(0, backend.iteration - 2))
        print_summary(burn_in)
    else:
        print("No chain file found")


if __name__ == '__main__':
    main()
