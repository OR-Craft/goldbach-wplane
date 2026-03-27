"""
╔══════════════════════════════════════════════════════════════╗
║  GOLDBACH CORRELATION v4.1 — GOERTZEL + MULTICORE           ║
║                                                             ║
║  Verified: matches numpy FFT to 4 decimal places            ║
║  Key fixes: M+1 Goertzel iterations + /M normalization      ║
║                                                             ║
║  SETUP: pip install numpy scipy primesieve numba             ║
║         (numba optional but gives ~50x speedup)             ║
║                                                             ║
║  USAGE:                                                     ║
║    python3 goldbach_goertzel.py test     # self-test         ║
║    python3 goldbach_goertzel.py          # runs 10^9         ║
║    python3 goldbach_goertzel.py 1e8      # specific N        ║
╚══════════════════════════════════════════════════════════════╝
"""

import numpy as np
from math import log, gcd, pi, cos, sin, isqrt
import time, hashlib, json, sys, os
from multiprocessing import Pool, cpu_count
from scipy import stats

N_CORES = max(1, cpu_count() - 4)
CACHE_DIR = os.path.expanduser("~/goldbach_cache")

# ═══════════════════════════════════════════════════
# GOERTZEL ALGORITHM — VERIFIED CORRECT
#
# Key: M+1 total iterations (not M, not len(lam))
# Verified against numpy FFT: error < 1e-11
# ═══════════════════════════════════════════════════

def goertzel_bin(lam, k, M):
    """Pure Python Goertzel. Returns DFT bin X[k]."""
    w = 2.0 * pi * k / M
    coeff = 2.0 * cos(w)
    s1 = s2 = 0.0
    L = len(lam)
    for n in range(M + 1):  # M+1 iterations!
        x = float(lam[n]) if n < L else 0.0
        s0 = x + coeff * s1 - s2
        s2 = s1
        s1 = s0
    return complex(s1 - s2 * cos(w), s2 * sin(w))

try:
    from numba import njit

    @njit(cache=True)
    def goertzel_bin_fast(lam, k, M):
        """Numba Goertzel. ~50x faster than pure Python."""
        w = 2.0 * pi * k / M
        coeff = 2.0 * cos(w)
        s1 = 0.0
        s2 = 0.0
        L = len(lam)
        for n in range(M + 1):  # M+1 iterations!
            x = lam[n] if n < L else 0.0
            s0 = x + coeff * s1 - s2
            s2 = s1
            s1 = s0
        return complex(s1 - s2 * cos(w), s2 * sin(w))

    HAVE_NUMBA = True
    print("  numba: available")
except ImportError:
    HAVE_NUMBA = False
    print("  numba: not found (pip install numba for ~50x speedup)")
    goertzel_bin_fast = None

# ═══════════════════════════════════════════════════
# SIEVE
# ═══════════════════════════════════════════════════

def von_mangoldt(N):
    try:
        import primesieve
        print(f"  Sieve (primesieve) to {N:.2e}...", end="", flush=True)
        t0 = time.time()
        lam = np.zeros(N + 1, dtype=np.float64)
        for p in primesieve.primes(N):
            lp = log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
        print(f" {time.time()-t0:.1f}s")
        return lam
    except ImportError:
        print(f"  Sieve (numpy) to {N:.2e}...", end="", flush=True)
        t0 = time.time()
        is_p = np.ones(N + 1, dtype=np.bool_)
        is_p[0] = is_p[1] = False
        is_p[4::2] = False
        for i in range(3, isqrt(N) + 1, 2):
            if is_p[i]: is_p[i*i::2*i] = False
        lam = np.zeros(N + 1, dtype=np.float64)
        for p in np.where(is_p)[0]:
            lp = log(int(p)); pk = int(p)
            while pk <= N: lam[pk] = lp; pk *= int(p)
        del is_p
        print(f" {time.time()-t0:.1f}s")
        return lam

# ═══════════════════════════════════════════════════
# MAJOR-ARC CLASSIFICATION
# ═══════════════════════════════════════════════════

def get_major_indices(M, Q_max=20):
    s = set()
    for q in range(1, Q_max + 1):
        for a in range(0, q + 1):
            if gcd(a, q) != 1 and a > 0: continue
            kc = round(a * M / q) % M
            for dk in range(-2, 3): s.add((kc + dk) % M)
    return np.array(sorted(s), dtype=np.int64)

# ═══════════════════════════════════════════════════
# PARALLEL GOERTZEL WORKER
# ═══════════════════════════════════════════════════

_LAM_SHARED = None
_M_SHARED = None

def _init_worker(lam_path, M):
    global _LAM_SHARED, _M_SHARED
    _LAM_SHARED = np.load(lam_path, mmap_mode='r')
    _M_SHARED = M

def _goertzel_worker(k):
    global _LAM_SHARED, _M_SHARED
    lam = _LAM_SHARED
    M = _M_SHARED

    if HAVE_NUMBA:
        try:
            result = goertzel_bin_fast(np.ascontiguousarray(lam), k, M)
            return int(k), result
        except:
            pass

    # Fallback: pure Python
    w = 2.0 * pi * k / M
    coeff = 2.0 * cos(w)
    s1 = s2 = 0.0
    L = len(lam)
    for n in range(M + 1):  # M+1!
        x = float(lam[n]) if n < L else 0.0
        s0 = x + coeff * s1 - s2
        s2 = s1
        s1 = s0
    return int(k), complex(s1 - s2 * cos(w), s2 * sin(w))

# ═══════════════════════════════════════════════════
# MAIN COMPUTATION
# ═══════════════════════════════════════════════════

def compute_correlation(N_center, window=2000, scan_step=20,
                        Q_max=20, n_cores=None, cache_dir=None):
    if n_cores is None: n_cores = N_CORES
    if cache_dir is None: cache_dir = CACHE_DIR
    os.makedirs(cache_dir, exist_ok=True)

    t_start = time.time()
    N_max = N_center + window + 2
    M = 2 * N_max + 1

    print(f"\n{'█'*60}")
    print(f"  GOERTZEL MULTICORE v4.1: N = {N_center:.2e}")
    print(f"  M = {M}  |  Cores: {n_cores}")
    print(f"  RAM: ~{(N_max+1)*8/1e9:.1f} GB")
    print(f"{'█'*60}\n")

    # ── Step 1: Sieve ──
    print("  [1/4] Sieve...")
    lam = von_mangoldt(N_max)
    psi = float(np.sum(lam))
    print(f"    ψ ratio: {psi/N_max:.6f}\n")

    # ── Step 2: Direct totals via np.dot ──
    lo = N_center - window + ((N_center - window) % 2)
    hi = N_center + window
    scan = [n for n in range(max(6, lo), hi + 1, scan_step)
            if n % 2 == 0]

    print(f"  [2/4] Direct totals ({len(scan)} N, np.dot)...")
    t0 = time.time()
    totals = []
    for i, Ni in enumerate(scan):
        total = float(np.dot(lam[1:Ni], lam[Ni-1:0:-1]))
        totals.append(total)
        if (i+1) % 50 == 0 or i == 0:
            el = time.time() - t0
            eta = el / (i+1) * (len(scan) - i - 1)
            print(f"    {i+1}/{len(scan)}"
                  f"  r_Λ({Ni})={total:.0f}"
                  f"  [{el:.0f}s, ~{eta:.0f}s left]")
    t_totals = time.time() - t0
    print(f"    Done: {t_totals:.1f}s\n")

    # ── Step 3: Goertzel parallel ──
    major_idx = get_major_indices(M, Q_max)
    print(f"  [3/4] Goertzel at {len(major_idx)} frequencies,"
          f" {n_cores} cores")

    # Save Λ for memory-mapped worker access
    lam_path = os.path.join(cache_dir, f"lam_temp_{N_center}.npy")
    np.save(lam_path, lam)
    del lam

    # JIT warmup
    if HAVE_NUMBA:
        print("    JIT warmup...", end="", flush=True)
        tiny = np.zeros(100, dtype=np.float64)
        tiny[2] = log(2); tiny[3] = log(3)
        _ = goertzel_bin_fast(tiny, 1, 201)
        print(" done")

    t0 = time.time()
    print(f"    Computing...", flush=True)

    with Pool(processes=n_cores,
              initializer=_init_worker,
              initargs=(lam_path, M)) as pool:
        results_iter = pool.imap_unordered(
            _goertzel_worker, major_idx.tolist(),
            chunksize=max(1, len(major_idx) // (n_cores * 4)))
        S_values = {}
        done = 0
        for k_val, dft_val in results_iter:
            S_values[k_val] = dft_val
            done += 1
            if done % 100 == 0 or done == len(major_idx):
                el = time.time() - t0
                eta = el / done * (len(major_idx) - done)
                print(f"      {done}/{len(major_idx)}"
                      f" [{el:.0f}s, ~{eta:.0f}s left]")

    t_goertzel = time.time() - t0
    print(f"    Done: {t_goertzel:.1f}s\n")

    # Build S² array
    # DFT convention: X[k] = Σ x[n] exp(-2πink/M)
    # Analytic convention: S(k/M) = Σ Λ(n) exp(+2πink/M) = conj(X[k])
    S2_major = np.array([np.conj(S_values[int(k)])**2
                         for k in major_idx])

    try: os.remove(lam_path)
    except: pass

    # ── Step 4: Correlate ──
    # CRITICAL: divide by M to match np.dot normalization
    # np.dot gives r_Λ(N) = Σ Λ(m)Λ(N-m)
    # Goertzel integrand sum gives M × r_Λ(N)
    print(f"  [4/4] Correlation...")

    maj_c, min_c = [], []
    for i, Ni in enumerate(scan):
        phases = -2.0 * pi * Ni * major_idx.astype(np.float64) / M
        # /M is the normalization fix!
        major_re = float(np.sum(
            (S2_major * np.exp(1j * phases)).real
        )) / M
        minor_re = totals[i] - major_re
        maj_c.append(major_re)
        min_c.append(minor_re)

    ma, mi = np.array(maj_c), np.array(min_c)

    if np.std(ma) > 1e-10 and np.std(mi) > 1e-10:
        r, p_val = stats.pearsonr(ma, mi)
        if np.isnan(r): r, p_val = 0.0, 1.0
    else:
        r, p_val = 0.0, 1.0

    r_sp, _ = stats.spearmanr(ma, mi)
    frac_pos = float(np.mean(mi > 0))
    mean_tot = float(np.mean(totals))
    mean_maj = float(np.mean(ma))
    mean_min = float(np.mean(mi))

    total_time = time.time() - t_start
    tag = ('✓ REINFORCING' if r > 0.3 else '~ MODERATE' if r > 0.1
           else '○ WEAK' if r > -0.1 else '✗ OPPOSING')

    print(f"""
  ┌──────────────────────────────────────────────┐
  │  N = {N_center:>14}                           │
  │  Pearson  r = {r:>+8.4f}       {tag:>14}  │
  │  Spearman ρ = {r_sp:>+8.4f}                       │
  │  p-value    = {p_val:>10.2e}                   │
  │  n_scan     = {len(scan):>6}                       │
  │  maj/total  = {mean_maj/mean_tot:.4f}                      │
  │  frac(minor > 0) = {frac_pos:.3f}                    │
  │  Time: {total_time:.0f}s ({total_time/60:.1f}min)                     │
  │  Goertzel: {t_goertzel:.0f}s on {n_cores} cores                │
  └──────────────────────────────────────────────┘""")

    cs = hashlib.sha256(json.dumps(
        {'N': N_center, 'r': round(r, 6)}).encode()).hexdigest()[:16]
    print(f"  Checksum: {cs}\n")

    result = {
        'N': N_center, 'r': float(r), 'r_spearman': float(r_sp),
        'p_val': float(p_val), 'frac_min_pos': frac_pos,
        'mean_maj': mean_maj, 'mean_min': mean_min,
        'mean_tot': mean_tot,
        'n_scan': len(scan), 'window': window,
        'scan_step': scan_step, 'Q_max': Q_max,
        'n_cores': n_cores,
        'time_totals': t_totals, 'time_goertzel': t_goertzel,
        'time_total': total_time, 'checksum': cs,
    }
    rpath = os.path.join(cache_dir, f"result_goertzel_{N_center}.json")
    with open(rpath, 'w') as f: json.dump(result, f, indent=2)
    print(f"  Saved: {rpath}")
    return result

# ═══════════════════════════════════════════════════
# SELF-TESTS
# ═══════════════════════════════════════════════════

def self_test(cache_dir=None):
    if cache_dir is None: cache_dir = CACHE_DIR
    print(f"\n{'='*60}\n  SELF-TESTS\n{'='*60}\n")
    all_ok = True

    # Test 1: Goertzel vs FFT at multiple bins
    print("  Test 1: Goertzel vs numpy FFT...")
    lam = np.zeros(1001, dtype=np.float64)
    lam[2]=log(2); lam[3]=log(3); lam[4]=log(2)
    lam[5]=log(5); lam[7]=log(7); lam[8]=log(2)
    M = 2001
    F_ref = np.fft.fft(lam, n=M)

    max_err = 0
    for k in [0, 1, 10, 100, 500, 1000]:
        got = goertzel_bin(lam, k, M)
        err = abs(got - F_ref[k]) / max(abs(F_ref[k]), 1e-10)
        max_err = max(max_err, err)
    ok1 = max_err < 1e-8
    all_ok &= ok1
    print(f"    Python:  max err = {max_err:.2e}  {'✓' if ok1 else '✗'}")

    if HAVE_NUMBA:
        max_err_nb = 0
        for k in [0, 1, 10, 100, 500, 1000]:
            got = goertzel_bin_fast(lam, k, M)
            err = abs(got - F_ref[k]) / max(abs(F_ref[k]), 1e-10)
            max_err_nb = max(max_err_nb, err)
        ok1b = max_err_nb < 1e-8
        all_ok &= ok1b
        print(f"    Numba:   max err = {max_err_nb:.2e}  {'✓' if ok1b else '✗'}")

    # Test 2: Normalization check
    print("\n  Test 2: Normalization (Goertzel sum / M = np.dot)...")
    lam2 = von_mangoldt(500)
    M2 = 1001
    Ni = 400
    dot_val = float(np.dot(lam2[1:Ni], lam2[Ni-1:0:-1]))
    # Sum Goertzel integrand over ALL bins
    integ_sum = 0.0
    for k in range(M2):
        Gk = goertzel_bin(lam2, k, M2)
        Sk = np.conj(Gk)
        S2k = Sk**2
        twist = np.exp(-2j * pi * Ni * k / M2)
        integ_sum += (S2k * twist).real
    ratio = integ_sum / dot_val if abs(dot_val) > 0 else 0
    ok2 = abs(ratio - M2) < 0.01
    all_ok &= ok2
    print(f"    sum/dot = {ratio:.4f} (should be M={M2})"
          f"  {'✓' if ok2 else '✗'}")

    # Test 3: Full pipeline matches FFT result
    print("\n  Test 3: Full pipeline at N=10000...")
    print("    (FFT reference: r = +0.7686)")
    res = compute_correlation(10000, window=500, scan_step=10,
                              Q_max=20, n_cores=2, cache_dir=cache_dir)
    ok3 = res is not None and abs(res['r'] - 0.7686) < 0.01
    all_ok &= ok3
    print(f"    Got r = {res['r']:+.4f},"
          f" diff = {abs(res['r']-0.7686):.4f}"
          f"  {'✓' if ok3 else '✗'}")

    # Test 4: maj/total should be ~0.95 at N=10000
    ok4 = res is not None and 0.9 < res['mean_maj']/res['mean_tot'] < 1.0
    all_ok &= ok4
    print(f"    maj/total = {res['mean_maj']/res['mean_tot']:.4f}"
          f" (expect ~0.95)  {'✓' if ok4 else '✗'}")

    print(f"\n  {'ALL PASSED' if all_ok else '*** FAILURES ***'}\n")
    return all_ok

# ═══════════════════════════════════════════════════
# MULTI-N SWEEP
# ═══════════════════════════════════════════════════

def run_all(targets, **kwargs):
    results = []
    for N in targets:
        results.append(compute_correlation(N, **kwargs))
    print(f"\n{'█'*60}\n  SUMMARY\n{'█'*60}\n")
    print(f"  {'N':>14} {'r':>10} {'ρ':>10}"
          f" {'p':>12} {'min>0':>7} {'time':>7}")
    print("  " + "─" * 58)
    for res in results:
        print(f"  {res['N']:>14} {res['r']:>+10.4f}"
              f" {res['r_spearman']:>+10.4f}"
              f" {res['p_val']:>12.2e}"
              f" {res['frac_min_pos']:>7.3f}"
              f" {res['time_total']:>6.0f}s")
    rs = [r['r'] for r in results]
    print(f"\n  Range: [{min(rs):.4f}, {max(rs):.4f}]"
          f"  Mean: {np.mean(rs):.4f}")
    if len(results) >= 3:
        rt, pt = stats.pearsonr([log(r['N']) for r in results], rs)
        print(f"  Trend: {rt:+.4f} (p={pt:.2e})")
    print()
    return results

# ═══════════════════════════════════════════════════
# ENTRY POINT
# ═══════════════════════════════════════════════════

if __name__ == "__main__":
    print(f"  Cores: {cpu_count()} total, using {N_CORES}")
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        self_test()
    elif len(sys.argv) > 1:
        try:
            compute_correlation(int(float(sys.argv[1])))
        except ValueError:
            print("Usage: python3 goldbach_goertzel.py test | 1e8 | 1e9")
    else:
        compute_correlation(10**9)