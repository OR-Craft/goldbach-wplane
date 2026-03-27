"""
╔══════════════════════════════════════════════════════════════╗
║  GOLDBACH LARGE-N FFT CORRELATION  v3.0                     ║
║  Hybrid: compute once, save to disk, rescan freely          ║
║                                                             ║
║  WORKFLOW:                                                  ║
║    1. compute_and_save(N) — FFT + IFFT, saves to disk       ║
║       (expensive, do once per N)                            ║
║    2. scan_from_cache(N) — reload + correlation scan         ║
║       (cheap, rerun with different params)                  ║
║    3. quick_run(N) — does both in one call                  ║
║                                                             ║
║  MEMORY:                                                    ║
║    10^8 → ~8 GB peak, then frees to ~100 MB for scanning   ║
║    10^9 → ~40 GB peak, then frees to ~100 MB               ║
║                                                             ║
║  DISK (per N):                                              ║
║    ifft array:  ~16 bytes × 2N (10^9 → 16 GB file)         ║
║    S2_major:    ~10 KB                                      ║
║    metadata:    ~1 KB                                       ║
║                                                             ║
║  USAGE:                                                     ║
║    CACHE_DIR = "/path/to/ssd"   # set this first!           ║
║                                                             ║
║    # Full pipeline:                                         ║
║    quick_run(10**8)                                         ║
║    quick_run(10**9)                                         ║
║                                                             ║
║    # Or split compute/scan:                                 ║
║    compute_and_save(10**9)     # 30 min, 40 GB peak         ║
║    scan_from_cache(10**9)      # seconds, ~100 MB           ║
║    scan_from_cache(10**9, window=5000)   # different params ║
║    scan_from_cache(10**9, Q_max=10)      # different Q      ║
║                                                             ║
║  DEPENDENCIES: numpy, scipy                                 ║
║  OPTIONAL: pip install primesieve (10x faster sieve)        ║
╚══════════════════════════════════════════════════════════════╝
"""

import numpy as np
from math import log, gcd, isqrt
import time
import hashlib
import json
import os
import gc as GC
import warnings
from scipy import stats

warnings.filterwarnings('ignore')


# ╔════════════════════════════════════════════════╗
# ║  SET THIS TO YOUR SSD PATH                    ║
# ╚════════════════════════════════════════════════╝
CACHE_DIR = os.path.expanduser("~/goldbach_cache")


# =============================================================
# SIEVE
# =============================================================

def von_mangoldt(N):
    """Λ(n) for n=0..N. Tries primesieve, falls back to numpy."""
    try:
        import primesieve
        print(f"  Sieve (primesieve) to {N:.2e}...", end="", flush=True)
        t0 = time.time()
        lam = np.zeros(N + 1, dtype=np.float64)
        for p in primesieve.primes(N):
            lp = log(p)
            pk = p
            while pk <= N:
                lam[pk] = lp
                pk *= p
        elapsed = time.time() - t0
        n_primes = len(primesieve.primes(N))
        print(f" {n_primes} primes, {elapsed:.1f}s")
        return lam
    except ImportError:
        print(f"  Sieve (numpy) to {N:.2e}...", end="", flush=True)
        print(" (pip install primesieve for 10x speedup)")
        t0 = time.time()
        is_p = np.ones(N + 1, dtype=np.bool_)
        is_p[0] = is_p[1] = False
        is_p[4::2] = False
        for i in range(3, isqrt(N) + 1, 2):
            if is_p[i]:
                is_p[i*i::2*i] = False
        lam = np.zeros(N + 1, dtype=np.float64)
        primes_found = 0
        for p in np.where(is_p)[0]:
            primes_found += 1
            lp = log(int(p))
            pk = int(p)
            while pk <= N:
                lam[pk] = lp
                pk *= int(p)
        del is_p
        GC.collect()
        elapsed = time.time() - t0
        print(f" {primes_found} primes, {elapsed:.1f}s")
        return lam


# =============================================================
# MAJOR-ARC CLASSIFICATION
# =============================================================

def get_major_indices(M, Q_max=20):
    """
    Find grid indices near a/q for all q ≤ Q_max.
    Returns sorted numpy array of indices.
    """
    major_set = set()
    for q in range(1, Q_max + 1):
        for a in range(0, q + 1):
            if gcd(a, q) != 1 and a > 0:
                continue
            kc = round(a * M / q) % M
            for dk in range(-2, 3):
                major_set.add((kc + dk) % M)
    return np.array(sorted(major_set), dtype=np.int64)


# =============================================================
# STEP 1: COMPUTE AND SAVE
# =============================================================

def compute_and_save(N_center, window=2000, Q_max=20,
                     cache_dir=None):
    """
    Compute FFT, IFFT, extract major-arc data, save to disk.

    This is the expensive step. Do it once per N.
    Peak memory: ~16 bytes × 4N (two complex arrays during FFT).

    Saves:
      {cache_dir}/ifft_{N}.npy     — IFFT(conj(S²)) array
      {cache_dir}/s2_major_{N}.npy — S² at major-arc points
      {cache_dir}/meta_{N}.json    — metadata (M, Q_max, etc.)
    """
    if cache_dir is None:
        cache_dir = CACHE_DIR
    os.makedirs(cache_dir, exist_ok=True)

    t_start = time.time()
    N_max = N_center + window + 2
    M = 2 * N_max + 1

    mem_peak_gb = M * 16 / 1e9  # two complex128 arrays
    disk_gb = M * 16 / 1e9      # one complex128 array on disk

    print()
    print("█" * 60)
    print(f"  COMPUTE & SAVE: N = {N_center:.2e}")
    print("█" * 60)
    print(f"  M = {M}")
    print(f"  Peak RAM:  ~{mem_peak_gb:.1f} GB")
    print(f"  Disk:      ~{disk_gb:.1f} GB")
    print(f"  Cache dir: {cache_dir}")
    print()

    # Check disk space
    try:
        st = os.statvfs(cache_dir)
        free_gb = st.f_bavail * st.f_frsize / 1e9
        print(f"  Disk free: {free_gb:.1f} GB")
        if free_gb < disk_gb * 1.2:
            print(f"  *** WARNING: Tight on disk! Need ~{disk_gb:.1f} GB ***")
    except:
        pass
    print()

    # Step 1: Build Λ(n)
    print("  [1/5] Building von Mangoldt array...")
    lam = von_mangoldt(N_max)
    psi = float(np.sum(lam))
    print(f"    ψ({N_max}) = {psi:.0f} (ratio: {psi/N_max:.4f})")
    print()

    # Step 2: FFT
    print(f"  [2/5] FFT (M = {M})...", end="", flush=True)
    t0 = time.time()
    F = np.fft.fft(lam, n=M)
    del lam
    GC.collect()
    t_fft = time.time() - t0
    print(f" {t_fft:.1f}s")

    # Step 3: S² and major-arc extraction
    print(f"  [3/5] Computing S² and extracting major arcs...",
          end="", flush=True)
    t0 = time.time()
    S = np.conj(F)
    S2 = S * S
    del S, F
    GC.collect()

    major_idx = get_major_indices(M, Q_max)
    S2_major = S2[major_idx].copy()
    print(f" {time.time()-t0:.1f}s ({len(major_idx)} major pts)")

    # Step 4: IFFT(conj(S²))
    print(f"  [4/5] IFFT for total integrals...", end="", flush=True)
    t0 = time.time()
    conj_S2 = np.conj(S2)
    del S2
    GC.collect()
    ifft_cs2 = np.fft.ifft(conj_S2)
    del conj_S2
    GC.collect()
    t_ifft = time.time() - t0
    print(f" {t_ifft:.1f}s")

    # Verify: total at N_center
    r_check = float(M * ifft_cs2[N_center].real)
    print(f"    TRUST: r_Λ({N_center}) = {r_check/M:.2f}")

    # Step 5: Save to disk
    print(f"  [5/5] Saving to disk...", end="", flush=True)
    t0 = time.time()

    ifft_path = os.path.join(cache_dir, f"ifft_{N_center}.npy")
    s2maj_path = os.path.join(cache_dir, f"s2_major_{N_center}.npy")
    meta_path = os.path.join(cache_dir, f"meta_{N_center}.json")
    majidx_path = os.path.join(cache_dir, f"major_idx_{N_center}.npy")

    np.save(ifft_path, ifft_cs2)
    np.save(s2maj_path, S2_major)
    np.save(majidx_path, major_idx)

    meta = {
        'N_center': N_center,
        'N_max': N_max,
        'M': M,
        'Q_max': Q_max,
        'window': window,
        'n_major': len(major_idx),
        'psi_ratio': psi / N_max,
        'r_check': r_check / M,
        'checksum': hashlib.sha256(
            f"{N_center}_{M}_{Q_max}".encode()
        ).hexdigest()[:16],
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S UTC',
                                    time.gmtime()),
    }
    with open(meta_path, 'w') as f:
        json.dump(meta, f, indent=2)

    del ifft_cs2, S2_major, major_idx
    GC.collect()

    t_save = time.time() - t0
    total_time = time.time() - t_start

    # File sizes
    sizes = {}
    for path in [ifft_path, s2maj_path, majidx_path, meta_path]:
        if os.path.exists(path):
            sizes[os.path.basename(path)] = os.path.getsize(path) / 1e9

    print(f" {t_save:.1f}s")
    print()
    print(f"  Files saved:")
    for name, gb in sizes.items():
        if gb > 0.001:
            print(f"    {name}: {gb:.2f} GB")
        else:
            print(f"    {name}: {gb*1000:.1f} MB")
    print()
    print(f"  Total time: {total_time:.1f}s")
    print(f"  RAM now: ~free (all big arrays deleted)")
    print()

    return meta


# =============================================================
# STEP 2: SCAN FROM CACHE
# =============================================================

def scan_from_cache(N_center, window=2000, scan_step=20,
                    Q_max=None, cache_dir=None):
    """
    Load cached IFFT + major-arc data, run correlation scan.

    This is cheap: ~100 MB RAM, runs in seconds.
    Can be called repeatedly with different parameters.

    Args:
        N_center: the N value (must have been compute_and_save'd)
        window: half-width of scan (can differ from compute step)
        scan_step: spacing between N values in scan
        Q_max: if set, recompute major indices (otherwise use cached)
        cache_dir: where the cache files live
    """
    if cache_dir is None:
        cache_dir = CACHE_DIR

    t_start = time.time()

    print()
    print("─" * 60)
    print(f"  SCAN FROM CACHE: N = {N_center:.2e}")
    print("─" * 60)

    # Load metadata
    meta_path = os.path.join(cache_dir, f"meta_{N_center}.json")
    if not os.path.exists(meta_path):
        print(f"  *** Cache not found for N={N_center}.")
        print(f"  *** Run compute_and_save({N_center}) first.")
        return None

    with open(meta_path) as f:
        meta = json.load(f)

    M = meta['M']
    cached_window = meta['window']
    cached_Q = meta['Q_max']

    print(f"  M = {M}, cached window = ±{cached_window},"
          f" cached Q = {cached_Q}")

    # Check window fits
    if window > cached_window:
        print(f"  *** Requested window {window} > cached {cached_window}")
        print(f"  *** Using window = {cached_window}")
        window = cached_window

    # Load IFFT array
    ifft_path = os.path.join(cache_dir, f"ifft_{N_center}.npy")
    print(f"  Loading IFFT ({os.path.getsize(ifft_path)/1e9:.2f} GB)...",
          end="", flush=True)
    t0 = time.time()
    ifft_cs2 = np.load(ifft_path)
    print(f" {time.time()-t0:.1f}s")

    # Load or recompute major indices
    if Q_max is not None and Q_max != cached_Q:
        print(f"  Recomputing major indices for Q_max={Q_max}...")
        major_idx = get_major_indices(M, Q_max)
        # Need to recompute S2_major from ifft... 
        # Actually we can't directly. We'd need the original S2.
        # So we must use cached Q_max.
        print(f"  *** Cannot change Q_max without recomputing FFT.")
        print(f"  *** Using cached Q_max = {cached_Q}")
        Q_max = cached_Q

    majidx_path = os.path.join(cache_dir, f"major_idx_{N_center}.npy")
    s2maj_path = os.path.join(cache_dir, f"s2_major_{N_center}.npy")
    major_idx = np.load(majidx_path)
    S2_major = np.load(s2maj_path)

    print(f"  Major: {len(major_idx)} pts, Q_max = {cached_Q}")

    # Build scan range
    lo = N_center - window + ((N_center - window) % 2)
    hi = N_center + window
    scan = [n for n in range(max(6, lo), hi + 1, scan_step)
            if n % 2 == 0]

    print(f"  Scanning {len(scan)} even N (step={scan_step})...",
          end="", flush=True)
    t0 = time.time()

    maj_c, min_c, tot_c = [], [], []

    for Ni in scan:
        # Total from IFFT: O(1) lookup
        total_re = float(M * ifft_cs2[Ni].real)

        # Major: sum over ~640 points
        phases = -2.0 * np.pi * Ni * major_idx / M
        major_re = float(np.sum((S2_major * np.exp(1j * phases)).real))

        # Minor = total - major
        minor_re = total_re - major_re

        maj_c.append(major_re)
        min_c.append(minor_re)
        tot_c.append(total_re)

    t_scan = time.time() - t0
    print(f" {t_scan:.1f}s")

    # Correlation
    ma, mi = np.array(maj_c), np.array(min_c)
    ta = np.array(tot_c)

    if np.std(ma) > 1e-10 and np.std(mi) > 1e-10:
        r, p_val = stats.pearsonr(ma, mi)
        if np.isnan(r):
            r, p_val = 0.0, 1.0
    else:
        r, p_val = 0.0, 1.0

    # Additional diagnostics
    # Fraction of N where minor is positive
    frac_min_pos = float(np.mean(np.array(min_c) > 0))

    # Mean contributions
    mean_maj = float(np.mean(ma))
    mean_min = float(np.mean(mi))
    mean_tot = float(np.mean(ta))

    # Spearman rank correlation (robust alternative)
    r_spearman, p_spearman = stats.spearmanr(ma, mi)

    del ifft_cs2
    GC.collect()

    total_time = time.time() - t_start

    # Report
    tag = '✓ REINFORCING' if r > 0.3 else \
          '~ MODERATE' if r > 0.1 else \
          '○ WEAK' if r > -0.1 else '✗ OPPOSING'

    print()
    print(f"  ┌──────────────────────────────────────────────┐")
    print(f"  │  N = {N_center:>14}                           │")
    print(f"  │  Pearson  r = {r:>+8.4f}       {tag:>14}  │")
    print(f"  │  Spearman ρ = {r_spearman:>+8.4f}                       │")
    print(f"  │  p-value    = {p_val:>10.2e}                   │")
    print(f"  │  n_scan     = {len(scan):>6}                       │")
    print(f"  │  mean(major) / mean(total) = "
          f"{mean_maj/mean_tot:.4f}              │"
          if abs(mean_tot) > 0 else
          f"  │  mean(total) ≈ 0                              │")
    print(f"  │  frac(minor > 0) = {frac_min_pos:.3f}                    │")
    print(f"  │  Scan time: {t_scan:.1f}s  Total: {total_time:.1f}s"
          f"                │")
    print(f"  └──────────────────────────────────────────────┘")

    cs = hashlib.sha256(json.dumps(
        {'N': N_center, 'r': round(r, 6), 'window': window,
         'step': scan_step}, sort_keys=True
    ).encode()).hexdigest()[:16]
    print(f"  Checksum: {cs}")
    print()

    return {
        'N': N_center, 'r': float(r), 'p_val': float(p_val),
        'r_spearman': float(r_spearman),
        'frac_min_pos': frac_min_pos,
        'mean_maj': mean_maj, 'mean_min': mean_min,
        'mean_tot': mean_tot,
        'n_scan': len(scan), 'window': window,
        'scan_step': scan_step, 'Q_max': cached_Q,
        'time_scan': t_scan, 'time_total': total_time,
        'checksum': cs,
    }


# =============================================================
# CONVENIENCE: DO BOTH IN ONE CALL
# =============================================================

def quick_run(N_center, window=2000, scan_step=20, Q_max=20,
              cache_dir=None, force_recompute=False):
    """
    Compute (if needed) and scan in one call.

    If cache exists and force_recompute is False, skips FFT.
    """
    if cache_dir is None:
        cache_dir = CACHE_DIR

    meta_path = os.path.join(cache_dir, f"meta_{N_center}.json")

    if os.path.exists(meta_path) and not force_recompute:
        print(f"  Cache found for N={N_center}. Skipping FFT.")
    else:
        compute_and_save(N_center, window=window, Q_max=Q_max,
                         cache_dir=cache_dir)

    return scan_from_cache(N_center, window=window,
                           scan_step=scan_step, cache_dir=cache_dir)


# =============================================================
# MULTI-N SWEEP
# =============================================================

def run_all(targets, window=2000, scan_step=20, Q_max=20,
            cache_dir=None):
    """Run at multiple N values and summarise."""
    results = []
    for N in targets:
        GC.collect()
        res = quick_run(N, window=window, scan_step=scan_step,
                        Q_max=Q_max, cache_dir=cache_dir)
        if res:
            results.append(res)

    if not results:
        print("No results.")
        return results

    print()
    print("█" * 60)
    print("  SUMMARY")
    print("█" * 60)
    print()
    print(f"  {'N':>14} {'Pearson r':>10} {'Spearman ρ':>11}"
          f" {'p-value':>12} {'min>0':>7}")
    print("  " + "─" * 56)
    for res in results:
        print(f"  {res['N']:>14} {res['r']:>+10.4f}"
              f" {res['r_spearman']:>+11.4f}"
              f" {res['p_val']:>12.2e}"
              f" {res['frac_min_pos']:>7.3f}")

    rs = [r['r'] for r in results]
    print()
    print(f"  Pearson range: [{min(rs):.4f}, {max(rs):.4f}]"
          f"  Mean: {np.mean(rs):.4f}")

    if len(results) >= 3:
        logNs = [log(r['N']) for r in results]
        rt, pt = stats.pearsonr(logNs, rs)
        print(f"  Trend vs log(N): {rt:+.4f} (p={pt:.2e})")
        if pt < 0.05:
            if rt > 0:
                print("  → STRENGTHENING with scale")
            else:
                print("  → WEAKENING with scale")
        else:
            print("  → STABLE (no significant trend)")

    print()

    if all(r > 0.3 for r in rs):
        print("  ✓ POSITIVE CORRELATION HOLDS at all tested scales.")
    elif all(r > 0 for r in rs):
        print("  ~ Positive but weakening.")
    else:
        print("  ✗ Correlation breaks at some scales.")
    print()

    return results


# =============================================================
# CACHE MANAGEMENT
# =============================================================

def list_cache(cache_dir=None):
    """Show what's been computed and cached."""
    if cache_dir is None:
        cache_dir = CACHE_DIR
    if not os.path.exists(cache_dir):
        print(f"  Cache dir {cache_dir} does not exist.")
        return

    print(f"  Cache directory: {cache_dir}")
    print()

    metas = sorted(f for f in os.listdir(cache_dir)
                   if f.startswith("meta_"))
    if not metas:
        print("  (empty)")
        return

    print(f"  {'N':>14} {'M':>12} {'Q':>4} {'r_Λ(N)':>12}"
          f" {'timestamp':>22}")
    print("  " + "─" * 68)

    total_size = 0
    for mf in metas:
        with open(os.path.join(cache_dir, mf)) as f:
            meta = json.load(f)
        N = meta['N_center']
        ifft_path = os.path.join(cache_dir, f"ifft_{N}.npy")
        size = os.path.getsize(ifft_path) if os.path.exists(ifft_path) else 0
        total_size += size
        print(f"  {N:>14} {meta['M']:>12} {meta['Q_max']:>4}"
              f" {meta.get('r_check', '?'):>12}"
              f" {meta.get('timestamp', '?'):>22}"
              f"  [{size/1e9:.1f} GB]")

    print()
    print(f"  Total cache size: {total_size/1e9:.1f} GB")


def clear_cache(N_center=None, cache_dir=None):
    """Delete cached files. If N_center is None, clears all."""
    if cache_dir is None:
        cache_dir = CACHE_DIR
    if not os.path.exists(cache_dir):
        return

    if N_center is not None:
        # Delete specific N
        for prefix in ['ifft_', 's2_major_', 'major_idx_', 'meta_']:
            path = os.path.join(cache_dir, f"{prefix}{N_center}.npy")
            if path.endswith('.npy') or path.endswith('.json'):
                if os.path.exists(path):
                    os.remove(path)
                    print(f"  Deleted {path}")
            path_json = os.path.join(cache_dir, f"{prefix}{N_center}.json")
            if os.path.exists(path_json):
                os.remove(path_json)
                print(f"  Deleted {path_json}")
    else:
        # Clear everything
        for f in os.listdir(cache_dir):
            os.remove(os.path.join(cache_dir, f))
            print(f"  Deleted {f}")


# =============================================================
# SELF-TEST
# =============================================================

def self_test(cache_dir=None):
    """Verify everything works at small scale."""
    if cache_dir is None:
        cache_dir = CACHE_DIR

    print("=" * 60)
    print("  SELF-TESTS")
    print("=" * 60)

    # Test 1: sieve
    lam = von_mangoldt(100)
    ok1 = abs(lam[7] - log(7)) < 0.01 and abs(lam[6]) < 0.01
    print(f"  Λ(7)=log7, Λ(6)=0:  {'✓' if ok1 else '✗'}")

    # Test 2: FFT cross-check
    lam_t = von_mangoldt(1000)
    M_t = 2 * 1000 + 1
    F_t = np.fft.fft(lam_t, n=M_t)
    S2_t = np.conj(F_t) ** 2
    ifft_t = np.fft.ifft(np.conj(S2_t))
    # ifft already has 1/M factor, so ifft[N] = r_Λ(N)/M... no:
    # ifft(F²)[N] = (1/M) Σ F²[k] e(2πiNk/M) = Σ_{a+b=N} Λ(a)Λ(b) = r_Λ(N)
    r_fft = float(ifft_t[1000].real)
    r_direct = sum(float(lam_t[m]) * float(lam_t[1000 - m])
                   for m in range(1, 1000))
    ok2 = abs(r_fft - r_direct) / max(abs(r_direct), 1) < 0.01
    print(f"  FFT vs direct:       {'✓' if ok2 else '✗'}"
          f"  (fft={r_fft:.2f}, direct={r_direct:.2f},"
          f" err={abs(r_fft-r_direct)/max(abs(r_direct),1):.2e})")

    # Test 3: full pipeline at N=10000
    print("  Full pipeline at N=10000...")
    res = quick_run(10000, window=500, scan_step=10, Q_max=20,
                    cache_dir=cache_dir)

    ok3 = res is not None and res['r'] > 0.5
    print(f"  r(10000) > 0.5:      {'✓' if ok3 else '✗'}"
          f"  (got {res['r']:.4f})" if res else "  (failed)")

    # Test 4: rescan from cache (should be instant)
    print("  Rescan from cache...")
    t0 = time.time()
    res2 = scan_from_cache(10000, window=500, scan_step=10,
                           cache_dir=cache_dir)
    t_rescan = time.time() - t0
    ok4 = res2 is not None and abs(res2['r'] - res['r']) < 1e-10
    print(f"  Cache rescan:        {'✓' if ok4 else '✗'}"
          f"  ({t_rescan:.2f}s, r matches)")

    all_ok = ok1 and ok2 and ok3 and ok4
    print()
    print(f"  {'ALL PASSED' if all_ok else '*** FAILURES ***'}")
    print()

    # Clean up test cache
    clear_cache(10000, cache_dir=cache_dir)

    return all_ok


# =============================================================
# NO AUTO-RUN
# =============================================================

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        self_test()
    elif len(sys.argv) > 1:
        try:
            targets = [int(float(x)) for x in sys.argv[1:]]
            run_all(targets)
        except Exception as e:
            print(f"Error: {e}")
            print("Usage: python script.py test")
            print("       python script.py 1e7 1e8 1e9")
    else:
        print("Loaded. Functions available:")
        print()
        print("  self_test()                # verify everything works")
        print("  quick_run(10**8)           # compute + scan in one go")
        print("  compute_and_save(10**9)    # expensive FFT, save to disk")
        print("  scan_from_cache(10**9)     # cheap rescan from disk")
        print("  run_all([1e7, 1e8, 1e9])   # full sweep")
        print("  list_cache()               # show cached computations")
        print("  clear_cache()              # delete cache files")
        print()
        print(f"  Cache dir: {CACHE_DIR}")
        print(f"  (change CACHE_DIR at top of script if needed)")
