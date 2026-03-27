"""
╔══════════════════════════════════════════════════════════════╗
║  CAVEAT 1 TEST v1.1: Mechanical Coupling vs Structure        ║
║                                                             ║
║  Fixed: Model C correctly identified as POSITIVE CONTROL    ║
║  Verdict based on NEGATIVE CONTROLS only (A, B, D, E)       ║
║                                                             ║
║  REQUIRES: goldbach_stats.py in the same directory          ║
║  USAGE: python3 goldbach_caveat1.py 1e9                     ║
╚══════════════════════════════════════════════════════════════╝
"""

import numpy as np
import json
import os
import sys
import time
from scipy import stats

CACHE_DIR = os.path.expanduser("~/goldbach_cache")


def run_caveat1_test(major, minor, total, n_synthetic=10000,
                     seed=42):
    rng = np.random.default_rng(seed)
    n = len(total)

    r_observed = stats.pearsonr(major, minor)[0]
    frac_arr = major / total
    observed_major_frac = np.mean(frac_arr)
    observed_major_std = np.std(frac_arr)

    print(f"\n{'█'*60}")
    print(f"  CAVEAT 1: MECHANICAL COUPLING TEST v1.1")
    print(f"{'█'*60}")
    print()
    print(f"  Observed r(major, minor) = {r_observed:+.4f}")
    print(f"  n = {n}")
    print(f"  mean(total) = {np.mean(total):.2f}")
    print(f"  mean(major/total) = {observed_major_frac:.4f}")
    print(f"  std(major/total)  = {observed_major_std:.6f}")
    print(f"  major sign: {'ALL NEGATIVE' if np.all(major < 0) else 'ALL POSITIVE' if np.all(major > 0) else 'MIXED'}")
    print(f"  minor sign: {'ALL POSITIVE' if np.all(minor > 0) else 'ALL NEGATIVE' if np.all(minor < 0) else 'MIXED'}")
    print()

    all_results = {}

    print("━" * 60)
    print("  NEGATIVE CONTROLS (destroy structure)")
    print("━" * 60)
    print()

    # ── Model A: Fixed fraction + noise ──
    print("─" * 60)
    print("  Model A: Fixed fraction + Gaussian noise")
    print("  (major = constant × total + noise)")
    print("─" * 60)

    f_obs = observed_major_frac
    noise_std = observed_major_std * np.mean(np.abs(total))

    r_synth_a = np.zeros(n_synthetic)
    for i in range(n_synthetic):
        noise = rng.normal(0, noise_std, size=n)
        major_synth = f_obs * total + noise
        minor_synth = total - major_synth
        r_synth_a[i] = stats.pearsonr(major_synth, minor_synth)[0]

    p_a = float(np.mean(r_synth_a >= r_observed))
    excess_a = r_observed - np.mean(r_synth_a)
    print(f"  Baseline r: {np.mean(r_synth_a):+.4f} ± {np.std(r_synth_a):.4f}")
    print(f"  Observed r: {r_observed:+.4f}")
    print(f"  Excess:     {excess_a:+.4f}")
    print(f"  p(synth >= obs) = {p_a:.6f}")
    print()

    all_results['model_a'] = {
        'name': 'Fixed fraction + noise', 'role': 'NEGATIVE CONTROL',
        'r_synth_mean': float(np.mean(r_synth_a)),
        'r_synth_std': float(np.std(r_synth_a)),
        'r_observed': float(r_observed),
        'p_value': p_a, 'excess': float(excess_a),
    }

    # ── Model B: Random fraction ──
    print("─" * 60)
    print("  Model B: Random fraction per N (uniform)")
    print("  (each N gets an independent random split)")
    print("─" * 60)

    frac_lo = np.min(frac_arr)
    frac_hi = np.max(frac_arr)
    frac_range = frac_hi - frac_lo
    frac_lo_gen = frac_lo - 0.1 * frac_range
    frac_hi_gen = frac_hi + 0.1 * frac_range

    r_synth_b = np.zeros(n_synthetic)
    for i in range(n_synthetic):
        f_random = rng.uniform(frac_lo_gen, frac_hi_gen, size=n)
        major_synth = f_random * total
        minor_synth = total - major_synth
        r_synth_b[i] = stats.pearsonr(major_synth, minor_synth)[0]

    p_b = float(np.mean(r_synth_b >= r_observed))
    excess_b = r_observed - np.mean(r_synth_b)
    print(f"  Frac range: [{frac_lo:.6f}, {frac_hi:.6f}]")
    print(f"  Baseline r: {np.mean(r_synth_b):+.4f} ± {np.std(r_synth_b):.4f}")
    print(f"  Observed r: {r_observed:+.4f}")
    print(f"  Excess:     {excess_b:+.4f}")
    print(f"  p(synth >= obs) = {p_b:.6f}")
    print()

    all_results['model_b'] = {
        'name': 'Random fraction (uniform)', 'role': 'NEGATIVE CONTROL',
        'r_synth_mean': float(np.mean(r_synth_b)),
        'r_synth_std': float(np.std(r_synth_b)),
        'r_observed': float(r_observed),
        'p_value': p_b, 'excess': float(excess_b),
    }

    # ── Model D: Constant fraction ──
    print("─" * 60)
    print("  Model D: Constant fraction (analytical)")
    print("  (if f were identical for all N)")
    print("─" * 60)

    mean_f = np.mean(frac_arr)
    if mean_f < 0 or mean_f > 1:
        r_constant_f = -1.0
        explanation = (f"f = {mean_f:.4f} is outside [0,1], "
                       "so major and minor have opposite signs. "
                       "Constant fraction -> r = -1.0")
    else:
        r_constant_f = +1.0
        explanation = (f"f = {mean_f:.4f} is in [0,1]. "
                       "Constant fraction -> r = +1.0")

    excess_d = r_observed - r_constant_f
    print(f"  mean(f) = {mean_f:.4f}")
    print(f"  {explanation}")
    print(f"  Analytical r: {r_constant_f:+.1f}")
    print(f"  Excess:       {excess_d:+.4f}")
    print()

    all_results['model_d'] = {
        'name': 'Constant fraction (analytical)', 'role': 'NEGATIVE CONTROL',
        'mean_fraction': float(mean_f),
        'r_analytical': float(r_constant_f),
        'excess': float(excess_d),
    }

    # ── Model E: Shuffled fractions ──
    print("─" * 60)
    print("  Model E: Shuffled fraction assignment")
    print("  (real fractions, randomly reassigned to N values)")
    print("─" * 60)

    r_frac_total, p_frac_total = stats.pearsonr(frac_arr, total)
    r_auto = stats.pearsonr(frac_arr[:-1], frac_arr[1:])[0]

    print(f"  r(fraction, total) = {r_frac_total:+.4f} (p={p_frac_total:.2e})")
    print(f"  Fraction autocorr  = {r_auto:+.4f}")

    r_synth_e = np.zeros(n_synthetic)
    for i in range(n_synthetic):
        f_shuf = rng.permutation(frac_arr)
        maj_shuf = f_shuf * total
        min_shuf = total - maj_shuf
        r_synth_e[i] = stats.pearsonr(maj_shuf, min_shuf)[0]

    p_e = float(np.mean(r_synth_e >= r_observed))
    excess_e = r_observed - np.mean(r_synth_e)
    print(f"  Baseline r: {np.mean(r_synth_e):+.4f} ± {np.std(r_synth_e):.4f}")
    print(f"  Observed r: {r_observed:+.4f}")
    print(f"  Excess:     {excess_e:+.4f}")
    print(f"  p(shuffled >= obs) = {p_e:.6f}")
    print()

    all_results['model_e'] = {
        'name': 'Shuffled fractions', 'role': 'NEGATIVE CONTROL',
        'r_frac_total': float(r_frac_total),
        'p_frac_total': float(p_frac_total),
        'frac_autocorr': float(r_auto),
        'r_synth_mean': float(np.mean(r_synth_e)),
        'r_synth_std': float(np.std(r_synth_e)),
        'p_value': float(p_e), 'excess': float(excess_e),
    }

    # ═══════════════════════════════════════════════════════
    # POSITIVE CONTROL
    # ═══════════════════════════════════════════════════════

    print("━" * 60)
    print("  POSITIVE CONTROL (preserve structure)")
    print("━" * 60)
    print()

    print("─" * 60)
    print("  Model C: Shuffled residuals from linear fit")
    print("  (preserves major-total relationship, shuffles deviations)")
    print("  This SHOULD reproduce r ~ observed. It's a calibration check.")
    print("─" * 60)

    slope, intercept = np.polyfit(total, major, 1)
    predicted_major = slope * total + intercept
    residuals = major - predicted_major

    r_synth_c = np.zeros(n_synthetic)
    for i in range(n_synthetic):
        shuffled_res = rng.permutation(residuals)
        major_synth = predicted_major + shuffled_res
        minor_synth = total - major_synth
        r_synth_c[i] = stats.pearsonr(major_synth, minor_synth)[0]

    p_c = float(np.mean(r_synth_c >= r_observed))
    diff_c = r_observed - np.mean(r_synth_c)
    calibrated = abs(diff_c) < 0.02

    print(f"  Linear fit: major = {slope:.6f} x total + {intercept:.2f}")
    print(f"  Baseline r: {np.mean(r_synth_c):+.4f} ± {np.std(r_synth_c):.4f}")
    print(f"  Observed r: {r_observed:+.4f}")
    print(f"  Difference: {diff_c:+.4f}")
    print(f"  p-value:    {p_c:.4f}")
    c_verdict = "CALIBRATED (code measures what it should)" if calibrated \
                else "MISCALIBRATED (investigate)"
    print(f"  {c_verdict}")
    print()

    all_results['model_c'] = {
        'name': 'Shuffled residuals (POSITIVE CONTROL)',
        'role': 'POSITIVE CONTROL',
        'r_synth_mean': float(np.mean(r_synth_c)),
        'r_synth_std': float(np.std(r_synth_c)),
        'r_observed': float(r_observed),
        'p_value': p_c, 'calibrated': calibrated,
    }

    # ═══════════════════════════════════════════════════════
    # SUMMARY
    # ═══════════════════════════════════════════════════════

    print()
    print("█" * 60)
    print("  CAVEAT 1 SUMMARY")
    print("█" * 60)
    print()
    print(f"  Observed r(major, minor) = {r_observed:+.4f}")
    print()
    print(f"  NEGATIVE CONTROLS (mechanical baselines):")
    print(f"    A (fixed frac + noise):  {all_results['model_a']['r_synth_mean']:+.4f}"
          f"  excess = {all_results['model_a']['excess']:+.4f}")
    print(f"    B (random frac):         {all_results['model_b']['r_synth_mean']:+.4f}"
          f"  excess = {all_results['model_b']['excess']:+.4f}")
    print(f"    D (constant f):          {all_results['model_d']['r_analytical']:+.1f}"
          f"     excess = {all_results['model_d']['excess']:+.4f}")
    print(f"    E (shuffled fracs):      {all_results['model_e']['r_synth_mean']:+.4f}"
          f"  excess = {all_results['model_e']['excess']:+.4f}")
    print()
    print(f"  POSITIVE CONTROL (calibration):")
    print(f"    C (shuffled residuals):  {all_results['model_c']['r_synth_mean']:+.4f}"
          f"  {'CALIBRATED' if calibrated else 'MISCALIBRATED'}")
    print()
    print(f"  STRUCTURAL INDICATOR:")
    print(f"    r(fraction, total) = {all_results['model_e']['r_frac_total']:+.4f}")
    print(f"    The split depends on the magnitude of the total.")
    print()

    # Verdict based on NEGATIVE CONTROLS ONLY
    neg_baselines = {
        'A': all_results['model_a']['r_synth_mean'],
        'B': all_results['model_b']['r_synth_mean'],
        'D': all_results['model_d']['r_analytical'],
        'E': all_results['model_e']['r_synth_mean'],
    }
    max_neg = max(neg_baselines.values())
    max_neg_name = max(neg_baselines, key=neg_baselines.get)
    excess = r_observed - max_neg

    neg_p = {k: all_results[f'model_{k.lower()}']['p_value']
             for k in ['A', 'B', 'E']}
    all_p_zero = all(p < 0.001 for p in neg_p.values())
    sign_reversal = max_neg < 0 and r_observed > 0

    print(f"  Highest negative-control baseline: {max_neg:+.4f}"
          f" (Model {max_neg_name})")
    print(f"  Excess over negative controls:     {excess:+.4f}")
    if sign_reversal:
        print(f"  *** SIGN REVERSAL: mechanics predicts negative,"
              f" observed positive ***")
    print(f"  All neg-control p < 0.001: {all_p_zero}")
    print()

    if sign_reversal and excess > 1.0 and all_p_zero:
        verdict = "STRUCTURAL"
        detail = (f"Mechanical coupling predicts r ~ {max_neg:+.4f}. "
                   f"Observed r = {r_observed:+.4f}. "
                   f"Excess of {excess:+.4f} with complete sign reversal. "
                   "The structural component OVERCOMES mechanical coupling "
                   "and reverses its direction. No mechanical model "
                   "reproduces the observed result. "
                   f"r(fraction, total) = {all_results['model_e']['r_frac_total']:+.4f} "
                   "confirms the split is governed by prime structure. "
                   "Positive control calibrated: "
                   f"{'yes' if calibrated else 'no'}.")
    elif excess > 0.5 and all_p_zero:
        verdict = "STRONGLY STRUCTURAL"
        detail = (f"Excess of {excess:+.4f} over all mechanical baselines. "
                   "Strong evidence for structural coupling.")
    elif excess > 0.2:
        verdict = "LIKELY STRUCTURAL"
        detail = (f"Excess of {excess:+.4f}. Evidence for structural "
                   "coupling beyond the constraint.")
    elif excess > 0.05:
        verdict = "MIXED"
        detail = ("Modest excess. Some structural component likely.")
    elif excess > -0.05:
        verdict = "MOSTLY MECHANICAL"
        detail = ("Observed r close to baselines. Mostly mechanical.")
    else:
        verdict = "MECHANICAL"
        detail = ("Observed r at or below baselines. Fully mechanical.")

    print(f"  ╔══════════════════════════════════════════════╗")
    print(f"  ║  VERDICT: {verdict:<36} ║")
    print(f"  ╚══════════════════════════════════════════════╝")
    print()
    print(f"  {detail}")
    print()

    all_results['verdict'] = verdict
    all_results['detail'] = detail
    all_results['r_observed'] = float(r_observed)
    all_results['excess_over_negative_controls'] = float(excess)
    all_results['sign_reversal'] = sign_reversal
    all_results['n'] = n
    all_results['positive_control_calibrated'] = calibrated

    save_path = os.path.join(CACHE_DIR, "caveat1_v1.1_results.json")
    with open(save_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"  Saved: {save_path}")

    return all_results


if __name__ == "__main__":
    N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**9

    print(f"\n  CAVEAT 1 TEST v1.1 for N = {N:.2e}")
    print(f"  Cache dir: {CACHE_DIR}")
    print(f"  Requires goldbach_stats.py in same directory.\n")

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from goldbach_stats import compute_major_minor

    print(f"  Computing major/minor decomposition (~12 min at 10^9)...\n")

    t0 = time.time()
    data = compute_major_minor(N, cache_dir=CACHE_DIR)
    print(f"\n  Computed in {time.time()-t0:.1f}s\n")

    run_caveat1_test(data['major'], data['minor'], data['total'])
