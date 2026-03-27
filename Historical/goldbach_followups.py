"""
╔══════════════════════════════════════════════════════════════╗
║  BEIGE MASTER FOLLOW-UPS                                    ║
║                                                             ║
║  1. Magnitude correlation: corr(|major|, |minor|)           ║
║     Tests "shared driver" hypothesis directly               ║
║                                                             ║
║  2. Fractional scaling: |major|/total and |minor|/total     ║
║     vs N — do relative contributions stabilise or drift?    ║
║                                                             ║
║  3. BONUS: Spectral coherence — corr(|major|, total)        ║
║     and corr(|minor|, total) separately                     ║
║                                                             ║
║  REQUIRES: goldbach_stats.py in same directory              ║
║  USAGE: python3 goldbach_followups.py 1e9                   ║
║         python3 goldbach_followups.py 1e8                   ║
║                                                             ║
║  Also runs at multiple N if you have cached data:           ║
║    python3 goldbach_followups.py multi                      ║
╚══════════════════════════════════════════════════════════════╝
"""

import numpy as np
import json
import os
import sys
import time
from scipy import stats

CACHE_DIR = os.path.expanduser("~/goldbach_cache")


def run_followups(major, minor, total, scan, N_center):
    """
    Run both Beige Master follow-up analyses plus bonus tests.
    """
    n = len(major)

    print(f"\n{'█'*60}")
    print(f"  BEIGE MASTER FOLLOW-UPS: N = {N_center:.2e}")
    print(f"{'█'*60}")
    print()
    print(f"  n = {n}")
    print(f"  r(major, minor) = {stats.pearsonr(major, minor)[0]:+.4f}")
    print()

    results = {}

    # ═══════════════════════════════════════════════════════
    # FOLLOW-UP 1: Magnitude correlation
    #
    # "If corr(|major|, |minor|) is also strongly positive,
    # it confirms the shared driver picture."
    #   — Beige Master
    # ═══════════════════════════════════════════════════════

    print("━" * 60)
    print("  FOLLOW-UP 1: Magnitude Correlation")
    print("  corr(|major|, |minor|)")
    print("━" * 60)
    print()

    abs_major = np.abs(major)
    abs_minor = np.abs(minor)

    r_mag, p_mag = stats.pearsonr(abs_major, abs_minor)
    rho_mag, p_rho_mag = stats.spearmanr(abs_major, abs_minor)

    print(f"  Pearson  r(|major|, |minor|) = {r_mag:+.4f}  (p={p_mag:.2e})")
    print(f"  Spearman ρ(|major|, |minor|) = {rho_mag:+.4f}  (p={p_rho_mag:.2e})")
    print()

    # Also check individual correlations with total
    r_maj_tot, p_maj_tot = stats.pearsonr(abs_major, total)
    r_min_tot, p_min_tot = stats.pearsonr(abs_minor, total)

    print(f"  r(|major|, total) = {r_maj_tot:+.4f}  (p={p_maj_tot:.2e})")
    print(f"  r(|minor|, total) = {r_min_tot:+.4f}  (p={p_min_tot:.2e})")
    print()

    # Interpretation
    if r_mag > 0.7:
        mag_interp = ("STRONG shared driver: when one component is large "
                       "in magnitude, the other is too. Both respond to "
                       "the same underlying spectral amplitude.")
    elif r_mag > 0.3:
        mag_interp = ("MODERATE shared driver: magnitudes co-move but "
                       "with substantial independent variation.")
    elif r_mag > 0:
        mag_interp = ("WEAK shared driver: slight positive tendency "
                       "but largely independent magnitudes.")
    else:
        mag_interp = ("NO shared driver: magnitudes are independent "
                       "or anti-correlated. The co-movement in the "
                       "signed correlation may be a sign effect only.")

    print(f"  INTERPRETATION: {mag_interp}")
    print()

    results['magnitude_correlation'] = {
        'r_abs_major_minor': float(r_mag),
        'p_abs_major_minor': float(p_mag),
        'rho_abs_major_minor': float(rho_mag),
        'r_abs_major_total': float(r_maj_tot),
        'r_abs_minor_total': float(r_min_tot),
        'interpretation': mag_interp,
    }

    # ═══════════════════════════════════════════════════════
    # FOLLOW-UP 2: Fractional Scaling
    #
    # "Plot |major|/total and |minor|/total versus N.
    # That shows whether the relative contributions
    # stabilise or drift."
    #   — Beige Master
    # ═══════════════════════════════════════════════════════

    print("━" * 60)
    print("  FOLLOW-UP 2: Fractional Scaling")
    print("  |major|/total and |minor|/total vs N")
    print("━" * 60)
    print()

    frac_abs_major = abs_major / total
    frac_abs_minor = abs_minor / total
    frac_signed_major = major / total
    frac_signed_minor = minor / total

    scan_arr = np.array(scan, dtype=np.float64)

    # Statistics of fractions
    print(f"  |major|/total:  mean={np.mean(frac_abs_major):.4f}"
          f"  std={np.std(frac_abs_major):.4f}"
          f"  CV={np.std(frac_abs_major)/np.mean(frac_abs_major):.4f}")
    print(f"  |minor|/total:  mean={np.mean(frac_abs_minor):.4f}"
          f"  std={np.std(frac_abs_minor):.4f}"
          f"  CV={np.std(frac_abs_minor)/np.mean(frac_abs_minor):.4f}")
    print(f"  major/total:    mean={np.mean(frac_signed_major):.4f}"
          f"  std={np.std(frac_signed_major):.4f}")
    print(f"  minor/total:    mean={np.mean(frac_signed_minor):.4f}"
          f"  std={np.std(frac_signed_minor):.4f}")
    print()

    # Trend: does the fraction drift with N?
    slope_maj, intercept_maj, r_trend_maj, p_trend_maj, _ = \
        stats.linregress(scan_arr, frac_abs_major)
    slope_min, intercept_min, r_trend_min, p_trend_min, _ = \
        stats.linregress(scan_arr, frac_abs_minor)

    print(f"  |major|/total trend with N:")
    print(f"    slope = {slope_maj:.2e}  r = {r_trend_maj:+.4f}"
          f"  p = {p_trend_maj:.2e}")
    print(f"  |minor|/total trend with N:")
    print(f"    slope = {slope_min:.2e}  r = {r_trend_min:+.4f}"
          f"  p = {p_trend_min:.2e}")
    print()

    # Are fractions stable?
    cv_major = np.std(frac_abs_major) / np.mean(frac_abs_major)
    cv_minor = np.std(frac_abs_minor) / np.mean(frac_abs_minor)

    if cv_major < 0.1 and cv_minor < 0.1:
        frac_interp = ("STABLE: Both fractions have CV < 10%. "
                        "The relative contributions are approximately "
                        "constant across the scan window.")
    elif cv_major < 0.3 and cv_minor < 0.3:
        frac_interp = ("MODERATELY STABLE: CV between 10-30%. "
                        "Some variation but no dramatic drift.")
    else:
        frac_interp = ("UNSTABLE: Large variation in fractional "
                        "contributions. The decomposition is not "
                        "settled at this scale.")

    # Is there a significant trend?
    if p_trend_maj < 0.01 or p_trend_min < 0.01:
        frac_interp += (" DRIFT DETECTED: significant linear trend "
                         "in at least one fraction.")
    else:
        frac_interp += (" NO DRIFT: no significant linear trend.")

    print(f"  INTERPRETATION: {frac_interp}")
    print()

    # Print a mini ASCII scatter for visual inspection
    print("  Mini profile (10 equally spaced samples):")
    step = max(1, n // 10)
    print(f"  {'N':>14}  {'|maj|/tot':>10}  {'|min|/tot':>10}  {'total':>14}")
    print(f"  {'─'*52}")
    for i in range(0, n, step):
        print(f"  {scan[i]:>14}  {frac_abs_major[i]:>10.4f}"
              f"  {frac_abs_minor[i]:>10.4f}  {total[i]:>14.0f}")
    print()

    results['fractional_scaling'] = {
        'mean_abs_major_frac': float(np.mean(frac_abs_major)),
        'std_abs_major_frac': float(np.std(frac_abs_major)),
        'cv_abs_major': float(cv_major),
        'mean_abs_minor_frac': float(np.mean(frac_abs_minor)),
        'std_abs_minor_frac': float(np.std(frac_abs_minor)),
        'cv_abs_minor': float(cv_minor),
        'slope_major': float(slope_maj),
        'r_trend_major': float(r_trend_maj),
        'p_trend_major': float(p_trend_maj),
        'slope_minor': float(slope_min),
        'r_trend_minor': float(r_trend_min),
        'p_trend_minor': float(p_trend_min),
        'interpretation': frac_interp,
    }

    # ═══════════════════════════════════════════════════════
    # FOLLOW-UP 3 (BONUS): Component-level analysis
    #
    # Break down exactly how major and minor relate to total
    # separately, to understand the shared driver mechanism.
    # ═══════════════════════════════════════════════════════

    print("━" * 60)
    print("  FOLLOW-UP 3 (BONUS): Component-Level Analysis")
    print("━" * 60)
    print()

    # Signed correlations with total
    r_maj_signed, p_ms = stats.pearsonr(major, total)
    r_min_signed, p_ns = stats.pearsonr(minor, total)

    print(f"  r(major, total) = {r_maj_signed:+.4f}  (p={p_ms:.2e})")
    print(f"  r(minor, total) = {r_min_signed:+.4f}  (p={p_ns:.2e})")
    print()

    # If major is negative and correlates negatively with total,
    # that means |major| correlates POSITIVELY with total
    # (larger total → more negative major → larger |major|)
    if major[0] < 0 and r_maj_signed < 0:
        print(f"  major is negative and anti-correlates with total.")
        print(f"  → |major| POSITIVELY correlates with total")
        print(f"     (confirmed: r(|major|, total) = {r_maj_tot:+.4f})")
    elif major[0] < 0 and r_maj_signed > 0:
        print(f"  major is negative but positively correlates with total.")
        print(f"  → |major| NEGATIVELY correlates with total")
        print(f"     (confirmed: r(|major|, total) = {r_maj_tot:+.4f})")
    print()

    print(f"  minor positively correlates with total:")
    print(f"    r(minor, total) = {r_min_signed:+.4f}")
    print()

    # The Beige Master's "A(N) and B(N)" decomposition
    # major ≈ -A(N), minor ≈ +B(N), total ≈ B(N) - A(N)
    # Check if A and B (magnitudes) correlate
    A = np.abs(major)  # A(N)
    B = minor.copy()   # B(N), already positive

    r_AB, p_AB = stats.pearsonr(A, B)
    print(f"  Beige Master decomposition:")
    print(f"    major = -A(N),  minor = +B(N),  total = B(N) - A(N)")
    print(f"    r(A, B) = {r_AB:+.4f}  (p={p_AB:.2e})")
    print()

    # What fraction of total variance is explained by A vs B?
    # total = B - A
    # Var(total) = Var(B) + Var(A) - 2Cov(A,B)
    var_A = np.var(A)
    var_B = np.var(B)
    cov_AB = np.cov(A, B)[0, 1]
    var_total = np.var(total)

    print(f"    Var(A)     = {var_A:.2e}")
    print(f"    Var(B)     = {var_B:.2e}")
    print(f"    Cov(A,B)   = {cov_AB:.2e}")
    print(f"    Var(total) = {var_total:.2e}")
    print(f"    Var(B) + Var(A) - 2Cov(A,B) = "
          f"{var_B + var_A - 2*cov_AB:.2e}")
    print(f"    (should equal Var(total): "
          f"ratio = {(var_B + var_A - 2*cov_AB)/var_total:.6f})")
    print()

    # Relative magnitudes
    mean_A = np.mean(A)
    mean_B = np.mean(B)
    mean_T = np.mean(total)
    print(f"    mean(A) = {mean_A:.2e}  ({mean_A/mean_T:.2f}× total)")
    print(f"    mean(B) = {mean_B:.2e}  ({mean_B/mean_T:.2f}× total)")
    print(f"    mean(T) = {mean_T:.2e}")
    print(f"    B/A ratio = {mean_B/mean_A:.4f}")
    print()

    # The key question: is the covariance between A and B
    # larger than what random co-movement would produce?
    # If Cov(A,B) > 0 and large, they share a driver.
    if cov_AB > 0:
        # Correlation of A and B relative to their total influence
        frac_shared = 2 * cov_AB / (var_A + var_B)
        print(f"    Shared variance fraction: "
              f"2*Cov(A,B)/(Var(A)+Var(B)) = {frac_shared:.4f}")
        if frac_shared > 0.5:
            comp_interp = ("STRONGLY COUPLED: A and B share more than "
                            "half their variance. They are driven by "
                            "the same underlying signal.")
        elif frac_shared > 0.2:
            comp_interp = ("MODERATELY COUPLED: Significant shared "
                            "variance between A and B.")
        else:
            comp_interp = ("WEAKLY COUPLED: Some shared variance "
                            "but mostly independent.")
    else:
        frac_shared = 0.0
        comp_interp = ("ANTI-COUPLED or INDEPENDENT: A and B do not "
                        "share positive covariance.")

    print(f"    INTERPRETATION: {comp_interp}")
    print()

    results['component_analysis'] = {
        'r_major_total_signed': float(r_maj_signed),
        'r_minor_total_signed': float(r_min_signed),
        'r_A_B': float(r_AB),
        'p_A_B': float(p_AB),
        'var_A': float(var_A),
        'var_B': float(var_B),
        'cov_AB': float(cov_AB),
        'var_total': float(var_total),
        'variance_check_ratio': float((var_B + var_A - 2*cov_AB)/var_total),
        'mean_A': float(mean_A),
        'mean_B': float(mean_B),
        'B_over_A': float(mean_B / mean_A),
        'shared_variance_fraction': float(frac_shared),
        'interpretation': comp_interp,
    }

    # ═══════════════════════════════════════════════════════
    # FOLLOW-UP 4 (BONUS): Modular structure check
    #
    # Do the correlations vary by N mod 6 or N mod 30?
    # This connects back to the earlier family analysis.
    # ═══════════════════════════════════════════════════════

    print("━" * 60)
    print("  FOLLOW-UP 4 (BONUS): Modular Family Check")
    print("━" * 60)
    print()

    scan_arr_int = np.array(scan)
    for mod_val in [6, 30]:
        print(f"  N mod {mod_val} analysis:")
        residues = scan_arr_int % mod_val
        unique_res = sorted(set(residues))
        rs_by_res = {}
        for res in unique_res:
            mask = residues == res
            if np.sum(mask) >= 5:
                r_res = stats.pearsonr(major[mask], minor[mask])[0]
                rs_by_res[res] = r_res
                print(f"    N ≡ {res:>2} mod {mod_val}: "
                      f"r = {r_res:+.4f}  (n={np.sum(mask)})")

        if rs_by_res:
            rs_vals = list(rs_by_res.values())
            print(f"    Range: [{min(rs_vals):.4f}, {max(rs_vals):.4f}]"
                  f"  Std: {np.std(rs_vals):.4f}")
            if np.std(rs_vals) < 0.1:
                print(f"    → UNIVERSAL (no family dependence)")
            else:
                print(f"    → FAMILY-DEPENDENT (some variation)")
        print()

    results['modular_families'] = {}
    for mod_val in [6, 30]:
        residues = scan_arr_int % mod_val
        unique_res = sorted(set(residues))
        family_rs = {}
        for res in unique_res:
            mask = residues == res
            if np.sum(mask) >= 5:
                family_rs[int(res)] = float(
                    stats.pearsonr(major[mask], minor[mask])[0])
        results['modular_families'][f'mod_{mod_val}'] = family_rs

    # ═══════════════════════════════════════════════════════
    # SUMMARY
    # ═══════════════════════════════════════════════════════

    print()
    print("█" * 60)
    print("  FOLLOW-UP SUMMARY")
    print("█" * 60)
    print()
    print(f"  1. Magnitude correlation:")
    print(f"     r(|major|, |minor|) = {r_mag:+.4f}")
    print(f"     → {'Shared driver CONFIRMED' if r_mag > 0.5 else 'Shared driver WEAK' if r_mag > 0 else 'No shared driver'}")
    print()
    print(f"  2. Fractional scaling:")
    print(f"     |major|/total CV = {cv_major:.4f}")
    print(f"     |minor|/total CV = {cv_minor:.4f}")
    print(f"     → {'STABLE' if cv_major < 0.1 and cv_minor < 0.1 else 'DRIFTING'}")
    print()
    print(f"  3. Component analysis:")
    print(f"     r(A, B) = {r_AB:+.4f}")
    print(f"     Shared variance = {frac_shared:.1%}")
    print(f"     → {'STRONGLY COUPLED' if frac_shared > 0.5 else 'MODERATELY COUPLED' if frac_shared > 0.2 else 'WEAKLY COUPLED'}")
    print()
    print(f"  4. Modular families:")
    for mod_val in [6, 30]:
        fam = results['modular_families'].get(f'mod_{mod_val}', {})
        if fam:
            vals = list(fam.values())
            print(f"     mod {mod_val}: range [{min(vals):.4f}, {max(vals):.4f}]"
                  f" std={np.std(vals):.4f}"
                  f" → {'UNIVERSAL' if np.std(vals) < 0.1 else 'VARIABLE'}")
    print()

    # Overall
    all_confirmed = (r_mag > 0.5 and
                     cv_major < 0.3 and cv_minor < 0.3 and
                     frac_shared > 0.2)

    if all_confirmed:
        print("  ╔══════════════════════════════════════════════╗")
        print("  ║  ALL FOLLOW-UPS CONFIRM SHARED DRIVER       ║")
        print("  ╚══════════════════════════════════════════════╝")
        print()
        print("  Both major and minor arc magnitudes co-move,")
        print("  their relative contributions are stable,")
        print("  and they share substantial variance.")
        print("  The prime Fourier spectrum has global coherence.")
    else:
        print("  Some follow-ups show mixed results.")
        print("  Further investigation needed.")
    print()

    # Save
    save_path = os.path.join(CACHE_DIR,
                             f"followups_{N_center}.json")
    with open(save_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"  Saved: {save_path}")

    return results


# ═══════════════════════════════════════════════════════
# MULTI-SCALE: run at all cached N values
# ═══════════════════════════════════════════════════════

def run_multi_scale(N_values, cache_dir=None):
    """Run follow-ups at multiple N and compare."""
    if cache_dir is None:
        cache_dir = CACHE_DIR

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from goldbach_stats import compute_major_minor

    all_results = []
    for N in N_values:
        print(f"\n{'═'*60}")
        print(f"  Computing for N = {N:.2e}...")
        data = compute_major_minor(N, cache_dir=cache_dir)
        res = run_followups(data['major'], data['minor'],
                            data['total'], data['scan'].tolist(), N)
        res['N'] = N
        all_results.append(res)

    # Cross-scale summary
    print(f"\n{'█'*60}")
    print(f"  CROSS-SCALE COMPARISON")
    print(f"{'█'*60}")
    print()
    print(f"  {'N':>14} {'r(|m|,|n|)':>11} {'CV_maj':>8} {'CV_min':>8}"
          f" {'r(A,B)':>8} {'shared%':>8}")
    print(f"  {'─'*60}")
    for res in all_results:
        mc = res['magnitude_correlation']
        fs = res['fractional_scaling']
        ca = res['component_analysis']
        print(f"  {res['N']:>14}"
              f" {mc['r_abs_major_minor']:>+11.4f}"
              f" {fs['cv_abs_major']:>8.4f}"
              f" {fs['cv_abs_minor']:>8.4f}"
              f" {ca['r_A_B']:>+8.4f}"
              f" {ca['shared_variance_fraction']:>7.1%}")
    print()

    return all_results


# ═══════════════════════════════════════════════════════
# ENTRY
# ═══════════════════════════════════════════════════════

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "multi":
        # Run at multiple scales
        targets = [10**7, 10**8, 10**9]
        run_multi_scale(targets)
    else:
        N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**9

        print(f"\n  Follow-ups for N = {N:.2e}")
        print(f"  Cache dir: {CACHE_DIR}")
        print(f"  Requires goldbach_stats.py in same directory.\n")

        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from goldbach_stats import compute_major_minor

        print(f"  Computing major/minor decomposition...\n")
        t0 = time.time()
        data = compute_major_minor(N, cache_dir=CACHE_DIR)
        print(f"\n  Computed in {time.time()-t0:.1f}s\n")

        run_followups(data['major'], data['minor'],
                      data['total'], data['scan'].tolist(), N)
