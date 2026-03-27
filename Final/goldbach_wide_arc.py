"""
╔══════════════════════════════════════════════════════════════╗
║  GOLDBACH WIDE ARC v5.0                                     ║
║                                                             ║
║  Fixes the grid neighbourhood issue from v4.                ║
║  Uses ±20 grid points (instead of ±2) to capture the full   ║
║  major arc peak structure.                                  ║
║                                                             ║
║  DIAGNOSIS: At ±2, we capture 12.5% of each major arc peak. ║
║  The phase twist e(-Nα) makes this subsample oscillate       ║
║  between positive and negative across N values.             ║
║  At ±20, the fraction converges to ~1.04 (stable, positive). ║
║                                                             ║
║  COST: ~5,200 Goertzel bins instead of 640.                 ║
║  Time at 10^9: ~50 min (12 cores) instead of ~6 min.        ║
║  But the absolute magnitudes are now trustworthy.            ║
║                                                             ║
║  ALSO COMPUTES:                                             ║
║  - Multiple half-widths (2, 5, 10, 20) for convergence test ║
║  - w-plane cloud: w(N) = major + i*minor                    ║
║  - Beige Master follow-ups with correct normalisation       ║
║                                                             ║
║  REQUIRES: goldbach_stats.py (for compute_major_minor)      ║
║  USAGE: python3 goldbach_wide_arc.py 1e9                    ║
║         python3 goldbach_wide_arc.py 1e8                    ║
║         python3 goldbach_wide_arc.py 1e7 --quick            ║
╚══════════════════════════════════════════════════════════════╝
"""

import numpy as np
from math import log, gcd, pi, cos, sin, isqrt
import time, hashlib, json, sys, os
from multiprocessing import Pool, cpu_count
from scipy import stats

N_CORES = max(1, cpu_count() - 4)
CACHE_DIR = os.path.expanduser("~/goldbach_cache")

# ── Goertzel ──

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
    print("  numba: available")
except ImportError:
    HAVE_NUMBA = False
    print("  numba: not found")

# ── Sieve ──

def von_mangoldt(N):
    try:
        import primesieve
        print(f"  Sieve to {N:.2e}...", end="", flush=True)
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
        is_p[0] = is_p[1] = False; is_p[4::2] = False
        for i in range(3, isqrt(N) + 1, 2):
            if is_p[i]: is_p[i*i::2*i] = False
        lam = np.zeros(N + 1, dtype=np.float64)
        for p in np.where(is_p)[0]:
            lp = log(int(p)); pk = int(p)
            while pk <= N: lam[pk] = lp; pk *= int(p)
        del is_p
        print(f" {time.time()-t0:.1f}s")
        return lam

# ── Major arc classification with variable width ──

def get_major_indices(M, Q_max=20, half_width=20):
    """Major arc grid points with specified neighbourhood width."""
    s = set()
    for q in range(1, Q_max + 1):
        for a in range(0, q + 1):
            if gcd(a, q) != 1 and a > 0: continue
            kc = round(a * M / q) % M
            for dk in range(-half_width, half_width + 1):
                s.add((kc + dk) % M)
    return np.array(sorted(s), dtype=np.int64)

# ── Parallel Goertzel ──

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
        except: pass
    # Fallback
    w = 2.0 * pi * k / M
    coeff = 2.0 * cos(w)
    s1 = s2 = 0.0
    L = len(lam)
    for n in range(M + 1):
        x = float(lam[n]) if n < L else 0.0
        s0 = x + coeff * s1 - s2
        s2 = s1; s1 = s0
    return int(k), complex(s1 - s2 * cos(w), s2 * sin(w))

# ── Core computation ──

def compute_wide_arc(N_center, window=2000, scan_step=20,
                     Q_max=20, half_widths=None,
                     n_cores=None, cache_dir=None):
    """
    Compute major/minor decomposition at multiple neighbourhood widths.
    Returns results for each half_width.
    """
    if half_widths is None:
        half_widths = [2, 5, 10, 20]
    if n_cores is None: n_cores = N_CORES
    if cache_dir is None: cache_dir = CACHE_DIR
    os.makedirs(cache_dir, exist_ok=True)

    t_start = time.time()
    N_max = N_center + window + 2
    M = 2 * N_max + 1

    print(f"\n{'█'*60}")
    print(f"  WIDE ARC ANALYSIS: N = {N_center:.2e}")
    print(f"  M = {M}  |  Cores: {n_cores}")
    print(f"  Half-widths to test: {half_widths}")
    print(f"{'█'*60}\n")

    # Step 1: Sieve
    print("  [1/4] Sieve...")
    lam = von_mangoldt(N_max)

    # Step 2: Totals
    lo = N_center - window + ((N_center - window) % 2)
    hi = N_center + window
    scan = [n for n in range(max(6, lo), hi + 1, scan_step) if n % 2 == 0]

    print(f"  [2/4] Direct totals ({len(scan)} N, np.dot)...")
    t0 = time.time()
    totals = []
    for i, Ni in enumerate(scan):
        totals.append(float(np.dot(lam[1:Ni], lam[Ni-1:0:-1])))
        if (i+1) % 50 == 0 or i == 0:
            el = time.time()-t0; eta = el/(i+1)*(len(scan)-i-1)
            print(f"    {i+1}/{len(scan)} [{el:.0f}s, ~{eta:.0f}s left]")
    print(f"    Done: {time.time()-t0:.1f}s\n")

    # Step 3: Goertzel for ALL unique indices across all half_widths
    # Compute the union of all major indices, then subset per width
    all_indices_by_hw = {}
    for hw in half_widths:
        all_indices_by_hw[hw] = get_major_indices(M, Q_max, hw)

    # Union of all indices
    union_set = set()
    for hw in half_widths:
        union_set.update(all_indices_by_hw[hw].tolist())
    union_idx = np.array(sorted(union_set), dtype=np.int64)

    print(f"  [3/4] Goertzel at {len(union_idx)} unique frequencies")
    print(f"    (union across half-widths {half_widths})")
    print(f"    Using {n_cores} cores...")

    # Save lam for workers
    lam_path = os.path.join(cache_dir, f"lam_wide_temp.npy")
    np.save(lam_path, lam)
    del lam

    if HAVE_NUMBA:
        print("    JIT warmup...", end="", flush=True)
        tiny = np.zeros(100, dtype=np.float64)
        tiny[2] = log(2)
        _ = goertzel_bin_fast(tiny, 1, 201)
        print(" done")

    t0 = time.time()
    with Pool(processes=n_cores,
              initializer=_init_worker,
              initargs=(lam_path, M)) as pool:
        S_values = {}
        done = 0
        total_bins = len(union_idx)
        for k_val, dft_val in pool.imap_unordered(
                _goertzel_worker, union_idx.tolist(),
                chunksize=max(1, total_bins // (n_cores * 4))):
            S_values[k_val] = dft_val
            done += 1
            if done % 500 == 0 or done == total_bins:
                el = time.time()-t0
                eta = el/done*(total_bins-done)
                print(f"      {done}/{total_bins} [{el:.0f}s, ~{eta:.0f}s left]")

    t_goertzel = time.time() - t0
    print(f"    Done: {t_goertzel:.1f}s\n")

    try: os.remove(lam_path)
    except: pass

    # Step 4: Compute decomposition for each half_width
    print(f"  [4/4] Computing decompositions...")

    all_results = {}

    for hw in half_widths:
        major_idx = all_indices_by_hw[hw]
        S2_major = np.array([np.conj(S_values[int(k)])**2
                             for k in major_idx])

        maj_c, min_c = [], []
        for i, Ni in enumerate(scan):
            phases = -2.0 * pi * Ni * major_idx.astype(np.float64) / M
            major_re = float(np.sum(
                (S2_major * np.exp(1j * phases)).real)) / M
            maj_c.append(major_re)
            min_c.append(totals[i] - major_re)

        ma = np.array(maj_c)
        mi = np.array(min_c)
        ta = np.array(totals)

        r_val, p_val = stats.pearsonr(ma, mi)
        r_sp, _ = stats.spearmanr(ma, mi)
        frac_pos = float(np.mean(mi > 0))
        mean_frac = float(np.mean(ma / ta))
        std_frac = float(np.std(ma / ta))

        # Magnitude correlation (sign-aware)
        # If major is consistently positive, use |major| directly
        # If major is consistently negative, flip sign for magnitude test
        maj_sign = 'POSITIVE' if np.mean(ma) > 0 else 'NEGATIVE'
        if np.mean(ma) < 0:
            r_mag, _ = stats.pearsonr(-ma, mi)  # compare magnitudes
        else:
            r_mag, _ = stats.pearsonr(ma, mi)

        all_results[hw] = {
            'half_width': hw,
            'n_major_pts': len(major_idx),
            'r_pearson': float(r_val),
            'r_spearman': float(r_sp),
            'p_value': float(p_val),
            'frac_minor_pos': frac_pos,
            'mean_major_frac': mean_frac,
            'std_major_frac': std_frac,
            'major_sign': maj_sign,
            'r_magnitude': float(r_mag),
            'major': ma,
            'minor': mi,
            'total': ta,
            'scan': np.array(scan),
        }

        tag = '✓' if r_val > 0.3 else '~' if r_val > 0 else '✗'
        print(f"    hw={hw:>3}: n_pts={len(major_idx):>6}"
              f"  r={r_val:>+.4f}  maj/tot={mean_frac:>+.4f}"
              f"  min>0={frac_pos:.3f}  {tag}")

    # ═══════════════════════════════════════════════════════
    # CONVERGENCE TEST
    # ═══════════════════════════════════════════════════════

    print(f"\n{'━'*60}")
    print(f"  CONVERGENCE TEST")
    print(f"{'━'*60}\n")

    print(f"  {'hw':>4} {'n_pts':>6} {'r':>8} {'maj/tot':>10}"
          f" {'std(frac)':>10} {'maj_sign':>10} {'min>0':>7}")
    print(f"  {'─'*58}")
    for hw in half_widths:
        res = all_results[hw]
        print(f"  {hw:>4} {res['n_major_pts']:>6} {res['r_pearson']:>+8.4f}"
              f" {res['mean_major_frac']:>+10.4f}"
              f" {res['std_major_frac']:>10.4f}"
              f" {res['major_sign']:>10}"
              f" {res['frac_minor_pos']:>7.3f}")

    # Check: does r converge?
    rs = [all_results[hw]['r_pearson'] for hw in half_widths]
    fracs = [all_results[hw]['mean_major_frac'] for hw in half_widths]
    print()
    print(f"  r range: [{min(rs):.4f}, {max(rs):.4f}]"
          f"  std: {np.std(rs):.4f}")
    print(f"  Fraction range: [{min(fracs):.4f}, {max(fracs):.4f}]")

    # Does the fraction stabilise (converge)?
    if len(half_widths) >= 3:
        last_two = [all_results[hw]['mean_major_frac'] for hw in half_widths[-2:]]
        frac_converged = abs(last_two[1] - last_two[0]) < 0.05
        print(f"  Fraction converged? {frac_converged}"
              f" (|Δ| = {abs(last_two[1]-last_two[0]):.4f})")
    print()

    # Does the major fraction stay positive at the widest hw?
    widest = half_widths[-1]
    widest_frac = all_results[widest]['mean_major_frac']
    widest_sign = all_results[widest]['major_sign']
    if widest_frac > 0.5:
        print(f"  ✓ At hw={widest}: major/total = {widest_frac:+.4f} (POSITIVE)")
        print(f"    The full major arc is positive as expected.")
        print(f"    Narrower hw values may flip sign due to subsampling.")
    elif widest_frac > 0:
        print(f"  ~ At hw={widest}: major/total = {widest_frac:+.4f} (weakly positive)")
    else:
        print(f"  ! At hw={widest}: major/total = {widest_frac:+.4f} (NEGATIVE)")
        print(f"    Even the wide arc is negative. May need wider hw.")

    # ═══════════════════════════════════════════════════════
    # W-PLANE CLOUD (using widest hw)
    # ═══════════════════════════════════════════════════════

    print(f"\n{'━'*60}")
    print(f"  W-PLANE ANALYSIS (hw={widest})")
    print(f"{'━'*60}\n")

    res_w = all_results[widest]
    w_real = res_w['major']  # Re(w) = major
    w_imag = res_w['minor']  # Im(w) = minor
    w_arg = np.arctan2(w_imag, w_real)  # argument
    w_abs = np.sqrt(w_real**2 + w_imag**2)  # modulus

    print(f"  w(N) = major(N) + i·minor(N)")
    print(f"  Using half_width = {widest}")
    print()
    print(f"  Re(w) = major:  mean={np.mean(w_real):.2e}  std={np.std(w_real):.2e}")
    print(f"  Im(w) = minor:  mean={np.mean(w_imag):.2e}  std={np.std(w_imag):.2e}")
    print(f"  |w|:            mean={np.mean(w_abs):.2e}")
    print(f"  arg(w):         mean={np.mean(w_arg):.4f} rad ({np.degrees(np.mean(w_arg)):.1f}°)")
    print(f"                  std ={np.std(w_arg):.4f} rad ({np.degrees(np.std(w_arg)):.1f}°)")
    print()

    # Quadrant distribution
    q1 = np.sum((w_real > 0) & (w_imag > 0))
    q2 = np.sum((w_real < 0) & (w_imag > 0))
    q3 = np.sum((w_real < 0) & (w_imag < 0))
    q4 = np.sum((w_real > 0) & (w_imag < 0))
    print(f"  Quadrant distribution (n={len(w_real)}):")
    print(f"    Q1 (+,+): {q1:>4} ({q1/len(w_real):.1%})")
    print(f"    Q2 (-,+): {q2:>4} ({q2/len(w_real):.1%})")
    print(f"    Q3 (-,-): {q3:>4} ({q3/len(w_real):.1%})")
    print(f"    Q4 (+,-): {q4:>4} ({q4/len(w_real):.1%})")
    print()

    # Distance from counterexample line Re+Im=0
    # For counterexample Z: major(Z) + minor(Z) = R(Z) = 0
    # So Re(w) + Im(w) = 0, i.e., the line y = -x
    dist_to_line = (w_real + w_imag) / np.sqrt(2)  # signed distance to y=-x
    print(f"  Distance to counterexample line (Re+Im=0):")
    print(f"    min:  {np.min(dist_to_line):.2e}")
    print(f"    mean: {np.mean(dist_to_line):.2e}")
    print(f"    All positive? {np.all(dist_to_line > 0)}")
    print()

    # arg(w) relative to π (counterexample needs arg = π + n·2π for some interpretation)
    # Actually: counterexample needs Re(w)+Im(w) = 0, not arg = π
    # The "danger angle" is arg = 3π/4 (the direction toward the line Re+Im=0)
    danger_angle = 3 * pi / 4  # 135°
    angular_dist = np.abs(w_arg - danger_angle)
    print(f"  Angular distance from 'danger direction' (135°):")
    print(f"    min:  {np.min(angular_dist):.4f} rad ({np.degrees(np.min(angular_dist)):.1f}°)")
    print(f"    mean: {np.mean(angular_dist):.4f} rad ({np.degrees(np.mean(angular_dist)):.1f}°)")
    print()

    # Modular family check in w-plane
    scan_arr = np.array(scan)
    print(f"  Modular family structure in w-plane:")
    for mod_val in [6]:
        residues = scan_arr % mod_val
        unique_res = sorted(set(residues))
        for res in unique_res:
            mask = residues == res
            if np.sum(mask) >= 5:
                mean_arg_res = np.mean(w_arg[mask])
                std_arg_res = np.std(w_arg[mask])
                mean_abs_res = np.mean(w_abs[mask])
                print(f"    N ≡ {res} mod {mod_val}: "
                      f"mean(arg)={np.degrees(mean_arg_res):>+7.1f}°"
                      f"  std={np.degrees(std_arg_res):>5.1f}°"
                      f"  mean(|w|)={mean_abs_res:.2e}")
    print()

    # ═══════════════════════════════════════════════════════
    # BEIGE MASTER FOLLOW-UPS (with corrected normalisation)
    # ═══════════════════════════════════════════════════════

    print(f"{'━'*60}")
    print(f"  FOLLOW-UPS (corrected, hw={widest})")
    print(f"{'━'*60}\n")

    ma_w = res_w['major']
    mi_w = res_w['minor']
    ta_w = res_w['total']

    # 1. Magnitude correlation (now with consistent positive major)
    if np.mean(ma_w) > 0:
        r_mag_corr, _ = stats.pearsonr(ma_w, mi_w)
        r_abs_mag, _ = stats.pearsonr(np.abs(ma_w), np.abs(mi_w))
        print(f"  Magnitude correlation:")
        print(f"    r(major, minor) = {r_mag_corr:+.4f}")
        print(f"    r(|major|, |minor|) = {r_abs_mag:+.4f}")
        print(f"    Major is positive → both should agree")
    else:
        print(f"  NOTE: Major still negative at hw={widest}")
        print(f"    Using sign-flipped magnitude test")
        r_abs_mag, _ = stats.pearsonr(-ma_w, mi_w)
        print(f"    r(-major, minor) = {r_abs_mag:+.4f}")
    print()

    # 2. A/B ratio stability (Beige Master request)
    A = np.abs(ma_w)
    B = mi_w if np.all(mi_w > 0) else np.abs(mi_w)
    ratio_AB = A / B
    print(f"  A(N)/B(N) ratio:")
    print(f"    mean = {np.mean(ratio_AB):.4f}")
    print(f"    std  = {np.std(ratio_AB):.4f}")
    print(f"    CV   = {np.std(ratio_AB)/np.mean(ratio_AB):.4f}")
    print(f"    → {'STABLE' if np.std(ratio_AB)/np.mean(ratio_AB) < 0.1 else 'VARIABLE'}")
    print()

    # 3. Shared variance
    r_AB, _ = stats.pearsonr(A, B)
    cov_AB = np.cov(A, B)[0, 1]
    shared = max(0, 2 * cov_AB / (np.var(A) + np.var(B)))
    print(f"  Shared variance:")
    print(f"    r(A, B) = {r_AB:+.4f}")
    print(f"    Shared fraction = {shared:.1%}")
    print()

    # ═══════════════════════════════════════════════════════
    # SUMMARY
    # ═══════════════════════════════════════════════════════

    total_time = time.time() - t_start

    print(f"\n{'█'*60}")
    print(f"  WIDE ARC SUMMARY")
    print(f"{'█'*60}\n")

    print(f"  N = {N_center:.2e}")
    print(f"  Widest half-width tested: {widest}")
    print(f"  Major arc fraction (hw={widest}): {all_results[widest]['mean_major_frac']:+.4f}")
    print(f"  r(major, minor) at hw={widest}: {all_results[widest]['r_pearson']:+.4f}")
    print()

    print(f"  CONVERGENCE:")
    for hw in half_widths:
        res = all_results[hw]
        print(f"    hw={hw:>3}: r={res['r_pearson']:>+.4f}  frac={res['mean_major_frac']:>+.4f}")
    print()

    print(f"  W-PLANE:")
    print(f"    Mean angle: {np.degrees(np.mean(w_arg)):.1f}°")
    print(f"    Q1 fraction: {q1/len(w_real):.1%}")
    print(f"    All R(N)>0: {np.all(dist_to_line > 0)}")
    print()

    print(f"  Time: {total_time:.0f}s ({total_time/60:.1f}min)")
    print(f"  Goertzel: {t_goertzel:.0f}s on {n_cores} cores")
    print()

    # Save
    save_data = {
        'N': N_center,
        'half_widths': half_widths,
        'convergence': {hw: {
            'r': all_results[hw]['r_pearson'],
            'frac': all_results[hw]['mean_major_frac'],
            'n_pts': all_results[hw]['n_major_pts'],
        } for hw in half_widths},
        'w_plane': {
            'mean_arg_deg': float(np.degrees(np.mean(w_arg))),
            'std_arg_deg': float(np.degrees(np.std(w_arg))),
            'q1_frac': float(q1/len(w_real)),
            'all_positive': bool(np.all(dist_to_line > 0)),
        },
        'time': total_time,
    }

    cs = hashlib.sha256(json.dumps(
        {'N': N_center, 'r': round(all_results[widest]['r_pearson'], 6)},
        sort_keys=True).encode()).hexdigest()[:16]
    save_data['checksum'] = cs

    save_path = os.path.join(cache_dir, f"wide_arc_{N_center}.json")
    with open(save_path, 'w') as f:
        json.dump(save_data, f, indent=2, default=str)
    print(f"  Saved: {save_path}")
    print(f"  Checksum: {cs}")

    return all_results


# ── Entry ──

if __name__ == "__main__":
    N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**9
    quick = '--quick' in sys.argv

    if quick:
        hws = [2, 10, 20]
    else:
        hws = [2, 5, 10, 20]

    print(f"\n  Wide Arc Analysis for N = {N:.2e}")
    print(f"  Half-widths: {hws}")
    print(f"  Cores: {N_CORES}")
    print(f"  Cache: {CACHE_DIR}\n")

    compute_wide_arc(N, half_widths=hws)
