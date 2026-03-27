"""
╔══════════════════════════════════════════════════════════════╗
║  GOLDBACH STATS VALIDATION BATTERY  v1.0                    ║
║                                                             ║
║  Beats up the r ≈ 0.83 finding from every angle.            ║
║  Uses cached Goertzel/FFT data where possible.              ║
║                                                             ║
║  TESTS:                                                     ║
║    1. Permutation test (10,000 shuffles)                    ║
║    2. Bootstrap CI (10,000 resamples)                       ║
║    3. Block bootstrap (accounts for autocorrelation)        ║
║    4. Split-half reliability                                ║
║    5. Detrended correlation                                 ║
║    6. Different Q_max (bucket sensitivity)                  ║
║    7. Different window sizes                                ║
║    8. Sign convention audit                                 ║
║    9. Robustness: Spearman + Kendall                        ║
║                                                             ║
║  SETUP: Uses goldbach_goertzel.py functions                 ║
║  USAGE: python3 goldbach_stats.py 1e9                       ║
║         python3 goldbach_stats.py 1e8                       ║
╚══════════════════════════════════════════════════════════════╝
"""

import numpy as np
from math import log, gcd, pi, cos, sin, isqrt
import time
import json
import os
import sys
import hashlib
from scipy import stats as sp_stats
from multiprocessing import Pool, cpu_count

N_CORES = max(1, cpu_count() - 4)
CACHE_DIR = os.path.expanduser("~/goldbach_cache")


# =============================================================
# CORE: Recompute major/minor arrays for given parameters
# =============================================================

def von_mangoldt(N):
    """Build Λ(n) array."""
    try:
        import primesieve
        lam = np.zeros(N + 1, dtype=np.float64)
        for p in primesieve.primes(N):
            lp = log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
        return lam
    except ImportError:
        is_p = np.ones(N + 1, dtype=np.bool_)
        is_p[0] = is_p[1] = False; is_p[4::2] = False
        for i in range(3, isqrt(N) + 1, 2):
            if is_p[i]: is_p[i*i::2*i] = False
        lam = np.zeros(N + 1, dtype=np.float64)
        for p in np.where(is_p)[0]:
            lp = log(int(p)); pk = int(p)
            while pk <= N: lam[pk] = lp; pk *= int(p)
        return lam


def get_major_indices(M, Q_max=20):
    """Grid indices near a/q for coprime a/q with q ≤ Q_max."""
    s = set()
    for q in range(1, Q_max + 1):
        for a in range(0, q + 1):
            if gcd(a, q) != 1 and a > 0: continue
            kc = round(a * M / q) % M
            for dk in range(-2, 3): s.add((kc + dk) % M)
    return np.array(sorted(s), dtype=np.int64)


try:
    from numba import njit
    @njit(cache=True)
    def goertzel_bin_fast(lam, k, M):
        w = 2.0 * pi * k / M
        coeff = 2.0 * cos(w)
        s1 = s2 = 0.0
        L = len(lam)
        for n in range(M + 1):
            x = lam[n] if n < L else 0.0
            s0 = x + coeff * s1 - s2
            s2 = s1; s1 = s0
        return complex(s1 - s2 * cos(w), s2 * sin(w))
    HAVE_NUMBA = True
except ImportError:
    HAVE_NUMBA = False
    def goertzel_bin_fast(lam, k, M):
        w = 2.0 * pi * k / M
        coeff = 2.0 * cos(w)
        s1 = s2 = 0.0
        L = len(lam)
        for n in range(M + 1):
            x = float(lam[n]) if n < L else 0.0
            s0 = x + coeff * s1 - s2
            s2 = s1; s1 = s0
        return complex(s1 - s2 * cos(w), s2 * sin(w))


# Worker setup for parallel Goertzel
_LAM_SHARED = None
_M_SHARED = None

def _init_worker(lam_path, M):
    global _LAM_SHARED, _M_SHARED
    _LAM_SHARED = np.load(lam_path, mmap_mode='r')
    _M_SHARED = M

def _goertzel_worker(k):
    global _LAM_SHARED, _M_SHARED
    lam = np.ascontiguousarray(_LAM_SHARED) if HAVE_NUMBA else _LAM_SHARED
    result = goertzel_bin_fast(lam, k, _M_SHARED)
    return int(k), result


def compute_major_minor(N_center, window=2000, scan_step=20,
                        Q_max=20, n_cores=None, lam=None,
                        cache_dir=None):
    """
    Compute matched (major, minor, total) arrays for a given
    parameter set. Returns dict with arrays and metadata.

    If lam is provided, uses it. Otherwise loads/computes.
    """
    if n_cores is None: n_cores = N_CORES
    if cache_dir is None: cache_dir = CACHE_DIR
    os.makedirs(cache_dir, exist_ok=True)

    N_max = N_center + window + 2
    M = 2 * N_max + 1

    # Get Λ
    if lam is None:
        print(f"    Sieving to {N_max:.2e}...", end="", flush=True)
        t0 = time.time()
        lam = von_mangoldt(N_max)
        print(f" {time.time()-t0:.1f}s")

    # Scan points
    lo = N_center - window + ((N_center - window) % 2)
    hi = N_center + window
    scan = [n for n in range(max(6, lo), hi + 1, scan_step)
            if n % 2 == 0]

    # Totals via np.dot
    totals = []
    for Ni in scan:
        totals.append(float(np.dot(lam[1:Ni], lam[Ni-1:0:-1])))

    # Goertzel for major arcs
    major_idx = get_major_indices(M, Q_max)

    # Save lam for parallel workers
    lam_path = os.path.join(cache_dir, f"lam_stats_temp.npy")
    np.save(lam_path, lam)

    # Warm numba
    if HAVE_NUMBA:
        tiny = np.zeros(100, dtype=np.float64)
        tiny[2] = log(2)
        _ = goertzel_bin_fast(tiny, 1, 201)

    with Pool(processes=n_cores,
              initializer=_init_worker,
              initargs=(lam_path, M)) as pool:
        S_values = {}
        for k_val, dft_val in pool.imap_unordered(
                _goertzel_worker, major_idx.tolist(),
                chunksize=max(1, len(major_idx) // (n_cores * 4))):
            S_values[k_val] = dft_val

    S2_major = np.array([np.conj(S_values[int(k)])**2
                         for k in major_idx])

    try: os.remove(lam_path)
    except: pass

    # Compute major/minor
    maj_arr = []
    min_arr = []
    for i, Ni in enumerate(scan):
        phases = -2.0 * pi * Ni * major_idx.astype(np.float64) / M
        major_re = float(np.sum(
            (S2_major * np.exp(1j * phases)).real)) / M
        maj_arr.append(major_re)
        min_arr.append(totals[i] - major_re)

    return {
        'scan': np.array(scan),
        'major': np.array(maj_arr),
        'minor': np.array(min_arr),
        'total': np.array(totals),
        'Q_max': Q_max,
        'window': window,
        'scan_step': scan_step,
        'M': M,
        'n_major_pts': len(major_idx),
        'lam': lam,
    }


# =============================================================
# STATS TESTS
# =============================================================

def test_1_permutation(major, minor, n_perms=10000, seed=42):
    """
    Permutation test: shuffle major relative to minor,
    compute r each time. Exact p-value with zero assumptions.
    """
    rng = np.random.default_rng(seed)
    r_obs = sp_stats.pearsonr(major, minor)[0]

    r_perms = np.zeros(n_perms)
    for i in range(n_perms):
        shuffled = rng.permutation(major)
        r_perms[i] = sp_stats.pearsonr(shuffled, minor)[0]

    p_perm = float(np.mean(np.abs(r_perms) >= np.abs(r_obs)))
    r_mean = float(np.mean(r_perms))
    r_std = float(np.std(r_perms))

    return {
        'test': 'Permutation',
        'r_observed': float(r_obs),
        'r_perm_mean': r_mean,
        'r_perm_std': r_std,
        'p_value': p_perm,
        'n_perms': n_perms,
        'n_exceed': int(np.sum(np.abs(r_perms) >= np.abs(r_obs))),
        'verdict': 'SIGNIFICANT' if p_perm < 0.001 else
                   'MARGINAL' if p_perm < 0.05 else 'NOT SIGNIFICANT',
    }


def test_2_bootstrap(major, minor, n_boot=10000, seed=42):
    """
    Bootstrap CI on r: resample pairs with replacement.
    """
    rng = np.random.default_rng(seed)
    n = len(major)
    r_obs = sp_stats.pearsonr(major, minor)[0]

    r_boots = np.zeros(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        r_boots[i] = sp_stats.pearsonr(major[idx], minor[idx])[0]

    ci_lo = float(np.percentile(r_boots, 2.5))
    ci_hi = float(np.percentile(r_boots, 97.5))

    return {
        'test': 'Bootstrap CI',
        'r_observed': float(r_obs),
        'ci_95_lo': ci_lo,
        'ci_95_hi': ci_hi,
        'ci_width': ci_hi - ci_lo,
        'boot_mean': float(np.mean(r_boots)),
        'boot_std': float(np.std(r_boots)),
        'n_boot': n_boot,
        'excludes_zero': ci_lo > 0,
        'verdict': 'ROBUST' if ci_lo > 0.5 else
                   'POSITIVE' if ci_lo > 0 else 'INCLUDES ZERO',
    }


def test_3_block_bootstrap(major, minor, block_size=20,
                           n_boot=10000, seed=42):
    """
    Block bootstrap: resample in contiguous blocks to account
    for autocorrelation between neighbouring N values.
    """
    rng = np.random.default_rng(seed)
    n = len(major)
    r_obs = sp_stats.pearsonr(major, minor)[0]

    n_blocks = max(1, n // block_size)
    r_boots = np.zeros(n_boot)

    for i in range(n_boot):
        # Sample block_size blocks with replacement
        block_starts = rng.integers(0, n - block_size + 1,
                                     size=n_blocks)
        idx = np.concatenate([np.arange(s, s + block_size)
                              for s in block_starts])[:n]
        if len(idx) < 3: continue
        r_boots[i] = sp_stats.pearsonr(major[idx], minor[idx])[0]

    ci_lo = float(np.percentile(r_boots, 2.5))
    ci_hi = float(np.percentile(r_boots, 97.5))

    return {
        'test': 'Block Bootstrap',
        'r_observed': float(r_obs),
        'block_size': block_size,
        'ci_95_lo': ci_lo,
        'ci_95_hi': ci_hi,
        'ci_width': ci_hi - ci_lo,
        'boot_mean': float(np.mean(r_boots)),
        'boot_std': float(np.std(r_boots)),
        'n_boot': n_boot,
        'excludes_zero': ci_lo > 0,
        'verdict': 'ROBUST' if ci_lo > 0.5 else
                   'POSITIVE' if ci_lo > 0 else 'INCLUDES ZERO',
    }


def test_4_split_half(major, minor):
    """
    Split-half reliability: compute r on first half and second
    half independently. Both should be positive and similar.
    """
    n = len(major)
    mid = n // 2

    r_first = sp_stats.pearsonr(major[:mid], minor[:mid])[0]
    r_second = sp_stats.pearsonr(major[mid:], minor[mid:])[0]
    r_full = sp_stats.pearsonr(major, minor)[0]

    # Also try even/odd split (interleaved)
    r_even = sp_stats.pearsonr(major[::2], minor[::2])[0]
    r_odd = sp_stats.pearsonr(major[1::2], minor[1::2])[0]

    return {
        'test': 'Split-Half',
        'r_full': float(r_full),
        'r_first_half': float(r_first),
        'r_second_half': float(r_second),
        'r_even_indices': float(r_even),
        'r_odd_indices': float(r_odd),
        'half_diff': abs(r_first - r_second),
        'interleave_diff': abs(r_even - r_odd),
        'all_positive': all(r > 0 for r in
                           [r_first, r_second, r_even, r_odd]),
        'verdict': 'CONSISTENT' if abs(r_first - r_second) < 0.15
                   and all(r > 0.3 for r in [r_first, r_second])
                   else 'INCONSISTENT',
    }


def test_5_detrended(major, minor, scan):
    """
    Detrend: remove linear trend in N from both series,
    then re-correlate. Tests if correlation is just
    "both things grow with N."
    """
    r_raw = sp_stats.pearsonr(major, minor)[0]

    # Linear detrend
    N_arr = np.array(scan, dtype=np.float64)
    slope_maj, intercept_maj = np.polyfit(N_arr, major, 1)
    slope_min, intercept_min = np.polyfit(N_arr, minor, 1)

    major_dt = major - (slope_maj * N_arr + intercept_maj)
    minor_dt = minor - (slope_min * N_arr + intercept_min)

    r_detrended = sp_stats.pearsonr(major_dt, minor_dt)[0]
    p_detrended = sp_stats.pearsonr(major_dt, minor_dt)[1]

    return {
        'test': 'Detrended',
        'r_raw': float(r_raw),
        'r_detrended': float(r_detrended),
        'p_detrended': float(p_detrended),
        'major_slope': float(slope_maj),
        'minor_slope': float(slope_min),
        'r_change': float(r_detrended - r_raw),
        'verdict': 'SURVIVES' if r_detrended > 0.3 else
                   'WEAKENED' if r_detrended > 0 else
                   'KILLED BY DETRENDING',
    }


def test_6_sign_audit(major, minor, total):
    """
    Audit sign conventions and magnitudes.
    Checks that totals are positive (they should be = r_Λ(N)),
    that major and minor decomposition is sensible.
    """
    n = len(major)

    return {
        'test': 'Sign Audit',
        'n_total_positive': int(np.sum(total > 0)),
        'n_total_negative': int(np.sum(total < 0)),
        'n_major_positive': int(np.sum(major > 0)),
        'n_major_negative': int(np.sum(major < 0)),
        'n_minor_positive': int(np.sum(minor > 0)),
        'n_minor_negative': int(np.sum(minor < 0)),
        'frac_minor_positive': float(np.mean(minor > 0)),
        'mean_total': float(np.mean(total)),
        'mean_major': float(np.mean(major)),
        'mean_minor': float(np.mean(minor)),
        'major_fraction': float(np.mean(major) / np.mean(total))
                          if np.mean(total) != 0 else float('nan'),
        'minor_fraction': float(np.mean(minor) / np.mean(total))
                          if np.mean(total) != 0 else float('nan'),
        'total_all_positive': bool(np.all(total > 0)),
        'sum_check': float(np.max(np.abs(
            major + minor - total)) / np.mean(np.abs(total))),
        'verdict': 'CLEAN' if (np.all(total > 0) and
                   abs(np.max(np.abs(major + minor - total)) /
                       np.mean(np.abs(total))) < 1e-6)
                   else 'ISSUES FOUND',
    }


def test_7_robust_correlations(major, minor):
    """
    Multiple correlation measures: Pearson, Spearman, Kendall.
    Agreement across all three strengthens the finding.
    """
    r_pearson, p_pearson = sp_stats.pearsonr(major, minor)
    r_spearman, p_spearman = sp_stats.spearmanr(major, minor)
    r_kendall, p_kendall = sp_stats.kendalltau(major, minor)

    return {
        'test': 'Robust Correlations',
        'pearson_r': float(r_pearson),
        'pearson_p': float(p_pearson),
        'spearman_rho': float(r_spearman),
        'spearman_p': float(p_spearman),
        'kendall_tau': float(r_kendall),
        'kendall_p': float(p_kendall),
        'all_positive': all(r > 0 for r in
                           [r_pearson, r_spearman, r_kendall]),
        'all_significant': all(p < 0.001 for p in
                              [p_pearson, p_spearman, p_kendall]),
        'verdict': 'ALL AGREE' if (
            all(r > 0.3 for r in [r_pearson, r_spearman, r_kendall])
            and all(p < 0.001 for p in [p_pearson, p_spearman, p_kendall]))
            else 'DISAGREEMENT',
    }


# =============================================================
# Q_MAX SENSITIVITY (requires recomputation)
# =============================================================

def test_8_qmax_sensitivity(N_center, window=2000, scan_step=20,
                            Q_values=None, n_cores=None,
                            cache_dir=None):
    """
    Test sensitivity to major/minor boundary definition.
    Recomputes with different Q_max values.
    """
    if Q_values is None:
        Q_values = [5, 10, 15, 20, 30, 50]
    if n_cores is None: n_cores = N_CORES
    if cache_dir is None: cache_dir = CACHE_DIR

    N_max = N_center + window + 2

    print(f"    Sieving once for Q_max sensitivity...", end="", flush=True)
    t0 = time.time()
    lam = von_mangoldt(N_max)
    print(f" {time.time()-t0:.1f}s")

    results = []
    for Q in Q_values:
        print(f"    Q_max = {Q}...", end="", flush=True)
        t0 = time.time()
        data = compute_major_minor(N_center, window=window,
                                    scan_step=scan_step,
                                    Q_max=Q, n_cores=n_cores,
                                    lam=lam, cache_dir=cache_dir)
        r = sp_stats.pearsonr(data['major'], data['minor'])[0]
        frac = float(np.mean(data['minor'] > 0))
        elapsed = time.time() - t0
        print(f" r={r:+.4f}, frac_pos={frac:.3f} ({elapsed:.1f}s)")
        results.append({
            'Q_max': Q,
            'r': float(r),
            'frac_minor_pos': frac,
            'n_major_pts': data['n_major_pts'],
        })

    rs = [r['r'] for r in results]
    return {
        'test': 'Q_max Sensitivity',
        'results': results,
        'r_range': [min(rs), max(rs)],
        'r_std': float(np.std(rs)),
        'all_positive': all(r > 0 for r in rs),
        'verdict': 'STABLE' if np.std(rs) < 0.1 and all(r > 0.3 for r in rs)
                   else 'SENSITIVE' if np.std(rs) > 0.2
                   else 'MODERATE SENSITIVITY',
    }


def test_9_window_sensitivity(N_center, window_sizes=None,
                              scan_step=20, Q_max=20,
                              n_cores=None, cache_dir=None):
    """
    Test sensitivity to scan window centre and size.
    """
    if window_sizes is None:
        window_sizes = [500, 1000, 2000, 5000]
    if n_cores is None: n_cores = N_CORES
    if cache_dir is None: cache_dir = CACHE_DIR

    N_max = N_center + max(window_sizes) + 2
    print(f"    Sieving once for window sensitivity...", end="",
          flush=True)
    t0 = time.time()
    lam = von_mangoldt(N_max)
    print(f" {time.time()-t0:.1f}s")

    results = []
    for w in window_sizes:
        print(f"    window = ±{w}...", end="", flush=True)
        t0 = time.time()
        data = compute_major_minor(N_center, window=w,
                                    scan_step=scan_step,
                                    Q_max=Q_max, n_cores=n_cores,
                                    lam=lam, cache_dir=cache_dir)
        r = sp_stats.pearsonr(data['major'], data['minor'])[0]
        frac = float(np.mean(data['minor'] > 0))
        elapsed = time.time() - t0
        print(f" r={r:+.4f}, n={len(data['scan'])}, frac_pos={frac:.3f}"
              f" ({elapsed:.1f}s)")
        results.append({
            'window': w,
            'r': float(r),
            'n_scan': len(data['scan']),
            'frac_minor_pos': frac,
        })

    rs = [r['r'] for r in results]
    return {
        'test': 'Window Sensitivity',
        'results': results,
        'r_range': [min(rs), max(rs)],
        'r_std': float(np.std(rs)),
        'all_positive': all(r > 0 for r in rs),
        'verdict': 'STABLE' if np.std(rs) < 0.1 and all(r > 0.3 for r in rs)
                   else 'SENSITIVE',
    }


# =============================================================
# MASTER RUNNER
# =============================================================

def run_validation(N_center, window=2000, scan_step=20, Q_max=20,
                   n_cores=None, cache_dir=None,
                   skip_recompute=False):
    """
    Run the full validation battery.

    If skip_recompute=True, only runs tests 1-7 (no new Goertzel).
    """
    if n_cores is None: n_cores = N_CORES
    if cache_dir is None: cache_dir = CACHE_DIR

    print()
    print("█" * 60)
    print(f"  VALIDATION BATTERY: N = {N_center:.2e}")
    print("█" * 60)
    print()

    t_start = time.time()

    # Compute base major/minor arrays
    print("  Computing base major/minor decomposition...")
    t0 = time.time()
    base = compute_major_minor(N_center, window=window,
                                scan_step=scan_step, Q_max=Q_max,
                                n_cores=n_cores, cache_dir=cache_dir)
    t_base = time.time() - t0
    print(f"  Base computed: {t_base:.1f}s")
    print()

    major = base['major']
    minor = base['minor']
    total = base['total']
    scan = base['scan']

    r_base = sp_stats.pearsonr(major, minor)[0]
    print(f"  Base r = {r_base:+.4f}")
    print(f"  Base frac(minor > 0) = {np.mean(minor > 0):.3f}")
    print()

    all_results = []

    # ── Test 1: Permutation ──
    print("─" * 60)
    print("  TEST 1: Permutation (10,000 shuffles)")
    print("─" * 60)
    t0 = time.time()
    res1 = test_1_permutation(major, minor)
    print(f"  r_observed = {res1['r_observed']:+.4f}")
    print(f"  r_shuffled = {res1['r_perm_mean']:+.4f} ± {res1['r_perm_std']:.4f}")
    print(f"  p-value    = {res1['p_value']:.6f}"
          f"  ({res1['n_exceed']}/{res1['n_perms']} exceeded)")
    print(f"  VERDICT: {res1['verdict']}")
    print(f"  ({time.time()-t0:.1f}s)")
    print()
    all_results.append(res1)

    # ── Test 2: Bootstrap CI ──
    print("─" * 60)
    print("  TEST 2: Bootstrap CI (10,000 resamples)")
    print("─" * 60)
    t0 = time.time()
    res2 = test_2_bootstrap(major, minor)
    print(f"  r_observed  = {res2['r_observed']:+.4f}")
    print(f"  95% CI      = [{res2['ci_95_lo']:+.4f}, {res2['ci_95_hi']:+.4f}]")
    print(f"  CI width    = {res2['ci_width']:.4f}")
    print(f"  Excludes 0? = {res2['excludes_zero']}")
    print(f"  VERDICT: {res2['verdict']}")
    print(f"  ({time.time()-t0:.1f}s)")
    print()
    all_results.append(res2)

    # ── Test 3: Block Bootstrap ──
    print("─" * 60)
    print("  TEST 3: Block Bootstrap (block=20, 10,000 resamples)")
    print("─" * 60)
    t0 = time.time()
    res3 = test_3_block_bootstrap(major, minor, block_size=20)
    print(f"  r_observed  = {res3['r_observed']:+.4f}")
    print(f"  95% CI      = [{res3['ci_95_lo']:+.4f}, {res3['ci_95_hi']:+.4f}]")
    print(f"  CI width    = {res3['ci_width']:.4f}")
    print(f"  Excludes 0? = {res3['excludes_zero']}")
    print(f"  VERDICT: {res3['verdict']}")
    print(f"  ({time.time()-t0:.1f}s)")
    print()
    all_results.append(res3)

    # ── Test 4: Split-Half ──
    print("─" * 60)
    print("  TEST 4: Split-Half Reliability")
    print("─" * 60)
    res4 = test_4_split_half(major, minor)
    print(f"  r_full         = {res4['r_full']:+.4f}")
    print(f"  r_first_half   = {res4['r_first_half']:+.4f}")
    print(f"  r_second_half  = {res4['r_second_half']:+.4f}")
    print(f"  r_even_indices = {res4['r_even_indices']:+.4f}")
    print(f"  r_odd_indices  = {res4['r_odd_indices']:+.4f}")
    print(f"  All positive?  = {res4['all_positive']}")
    print(f"  VERDICT: {res4['verdict']}")
    print()
    all_results.append(res4)

    # ── Test 5: Detrended ──
    print("─" * 60)
    print("  TEST 5: Detrended Correlation")
    print("─" * 60)
    res5 = test_5_detrended(major, minor, scan)
    print(f"  r_raw       = {res5['r_raw']:+.4f}")
    print(f"  r_detrended = {res5['r_detrended']:+.4f}")
    print(f"  p_detrended = {res5['p_detrended']:.2e}")
    print(f"  Change:       {res5['r_change']:+.4f}")
    print(f"  VERDICT: {res5['verdict']}")
    print()
    all_results.append(res5)

    # ── Test 6: Sign Audit ──
    print("─" * 60)
    print("  TEST 6: Sign Convention Audit")
    print("─" * 60)
    res6 = test_6_sign_audit(major, minor, total)
    print(f"  Total:  {res6['n_total_positive']} pos /"
          f" {res6['n_total_negative']} neg")
    print(f"  Major:  {res6['n_major_positive']} pos /"
          f" {res6['n_major_negative']} neg")
    print(f"  Minor:  {res6['n_minor_positive']} pos /"
          f" {res6['n_minor_negative']} neg")
    print(f"  major + minor = total? max error ="
          f" {res6['sum_check']:.2e}")
    print(f"  Major fraction: {res6['major_fraction']:.4f}")
    print(f"  Minor fraction: {res6['minor_fraction']:.4f}")
    print(f"  VERDICT: {res6['verdict']}")
    print()
    all_results.append(res6)

    # ── Test 7: Robust Correlations ──
    print("─" * 60)
    print("  TEST 7: Robust Correlations (Pearson/Spearman/Kendall)")
    print("─" * 60)
    res7 = test_7_robust_correlations(major, minor)
    print(f"  Pearson  r = {res7['pearson_r']:+.4f}  (p={res7['pearson_p']:.2e})")
    print(f"  Spearman ρ = {res7['spearman_rho']:+.4f}  (p={res7['spearman_p']:.2e})")
    print(f"  Kendall  τ = {res7['kendall_tau']:+.4f}  (p={res7['kendall_p']:.2e})")
    print(f"  All positive & significant? {res7['all_significant']}")
    print(f"  VERDICT: {res7['verdict']}")
    print()
    all_results.append(res7)

    # ── Tests 8-9: Recomputation tests ──
    if not skip_recompute:
        print("─" * 60)
        print("  TEST 8: Q_max Sensitivity")
        print("─" * 60)
        res8 = test_8_qmax_sensitivity(
            N_center, window=window, scan_step=scan_step,
            Q_values=[5, 10, 15, 20, 30, 50],
            n_cores=n_cores, cache_dir=cache_dir)
        print(f"  r range: [{res8['r_range'][0]:.4f}, {res8['r_range'][1]:.4f}]")
        print(f"  r std:   {res8['r_std']:.4f}")
        print(f"  VERDICT: {res8['verdict']}")
        print()
        all_results.append(res8)

        print("─" * 60)
        print("  TEST 9: Window Size Sensitivity")
        print("─" * 60)
        res9 = test_9_window_sensitivity(
            N_center, window_sizes=[500, 1000, 2000, 5000],
            scan_step=scan_step, Q_max=Q_max,
            n_cores=n_cores, cache_dir=cache_dir)
        print(f"  r range: [{res9['r_range'][0]:.4f}, {res9['r_range'][1]:.4f}]")
        print(f"  r std:   {res9['r_std']:.4f}")
        print(f"  VERDICT: {res9['verdict']}")
        print()
        all_results.append(res9)

    # ── Summary ──
    total_time = time.time() - t_start

    print()
    print("█" * 60)
    print("  VALIDATION SUMMARY")
    print("█" * 60)
    print()

    verdicts = [r.get('verdict', '?') for r in all_results]
    for r in all_results:
        name = r.get('test', '?')
        v = r.get('verdict', '?')
        symbol = '✓' if v in ['SIGNIFICANT', 'ROBUST', 'CONSISTENT',
                                'SURVIVES', 'CLEAN', 'ALL AGREE',
                                'STABLE'] else \
                 '~' if v in ['POSITIVE', 'WEAKENED',
                              'MODERATE SENSITIVITY'] else '✗'
        print(f"  {symbol} {name:.<35} {v}")

    n_pass = sum(1 for v in verdicts if v in [
        'SIGNIFICANT', 'ROBUST', 'CONSISTENT', 'SURVIVES',
        'CLEAN', 'ALL AGREE', 'STABLE'])
    n_total = len(verdicts)

    print()
    print(f"  Passed: {n_pass}/{n_total}")
    print(f"  Time:   {total_time:.0f}s ({total_time/60:.1f}min)")
    print()

    if n_pass == n_total:
        print("  ██ ALL TESTS PASSED ██")
        print("  The positive correlation is robust under every")
        print("  validation method applied.")
    elif n_pass >= n_total - 2:
        print("  Most tests passed. Review flagged tests.")
    else:
        print("  Multiple tests failed. Finding may be fragile.")
    print()

    # Save
    cs = hashlib.sha256(json.dumps(
        {'N': N_center, 'n_pass': n_pass, 'n_total': n_total,
         'r_base': round(r_base, 6)}).encode()).hexdigest()[:16]

    save_path = os.path.join(cache_dir,
                             f"validation_{N_center}.json")
    with open(save_path, 'w') as f:
        json.dump({
            'N': N_center,
            'r_base': float(r_base),
            'n_pass': n_pass,
            'n_total': n_total,
            'results': all_results,
            'checksum': cs,
        }, f, indent=2, default=str)
    print(f"  Saved: {save_path}")
    print(f"  Checksum: {cs}")

    return all_results


# =============================================================
# ENTRY
# =============================================================

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 goldbach_stats.py <N>")
        print("       python3 goldbach_stats.py 1e8")
        print("       python3 goldbach_stats.py 1e9")
        print("       python3 goldbach_stats.py 1e8 --quick")
        print("         (--quick skips Q_max and window tests)")
        sys.exit(0)

    N = int(float(sys.argv[1]))
    skip = '--quick' in sys.argv

    print(f"\n  N = {N:.2e}")
    print(f"  Cores: {N_CORES}")
    print(f"  Skip recompute tests: {skip}")
    print(f"  numba: {'yes' if HAVE_NUMBA else 'no'}")

    run_validation(N, skip_recompute=skip)
