"""
10^9 correlation WITHOUT FFT. Uses np.dot on views = no temp arrays.
Should take ~1-2 hours instead of 9.
"""
import numpy as np
from math import log, gcd, pi
import time, hashlib, json
from scipy import stats

def von_mangoldt(N):
    import primesieve
    print(f"  Sieve to {N:.2e}...", end="", flush=True)
    t0 = time.time()
    lam = np.zeros(N+1, dtype=np.float64)
    for p in primesieve.primes(N):
        lp = log(p); pk = p
        while pk <= N: lam[pk] = lp; pk *= p
    print(f" {time.time()-t0:.1f}s")
    return lam

def get_major_indices(M, Q_max=20):
    s = set()
    for q in range(1, Q_max+1):
        for a in range(0, q+1):
            if gcd(a,q)!=1 and a>0: continue
            kc = round(a*M/q) % M
            for dk in range(-2,3): s.add((kc+dk)%M)
    return np.array(sorted(s), dtype=np.int64)

def run_1e9(N_center=10**9, window=2000, scan_step=20, Q_max=20):
    t_start = time.time()
    N_max = N_center + window + 2
    M = 2 * N_max + 1

    print(f"\n{'='*60}")
    print(f"  FAST NO-FFT: N = {N_center:.2e}")
    print(f"  RAM: ~8 GB (one Λ array + views only)")
    print(f"{'='*60}\n")

    lam = von_mangoldt(N_max)
    psi = float(np.sum(lam))
    print(f"  psi ratio: {psi/N_max:.4f}")
    print(f"  Array: {lam.nbytes/1e9:.1f} GB\n")

    lo = N_center - window + ((N_center-window)%2)
    hi = N_center + window
    scan = [n for n in range(max(6,lo), hi+1, scan_step) if n%2==0]

    # Step 2: Direct totals using np.dot on VIEWS (no copies!)
    print(f"  [1/3] Direct totals for {len(scan)} N values")
    print(f"         using np.dot (BLAS, zero temp arrays)...")
    t0 = time.time()

    totals = []
    for i, Ni in enumerate(scan):
        # lam[1:Ni] and lam[Ni-1:0:-1] are VIEWS, no memory allocated
        # np.dot computes inner product without materializing product array
        total = float(np.dot(lam[1:Ni], lam[Ni-1:0:-1]))
        totals.append(total)
        if (i+1) % 10 == 0 or i == 0:
            elapsed = time.time() - t0
            eta = elapsed/(i+1) * (len(scan)-i-1)
            print(f"    {i+1}/{len(scan)}"
                  f"  r_L({Ni})={total:.0f}"
                  f"  [{elapsed:.0f}s, ~{eta:.0f}s left]")

    print(f"  Done: {time.time()-t0:.1f}s\n")

    # Step 3: Major-arc S^2 at specific frequencies
    major_idx = get_major_indices(M, Q_max)
    print(f"  [2/3] S(α)² at {len(major_idx)} major frequencies...")
    t0 = time.time()

    S2_major = np.zeros(len(major_idx), dtype=np.complex128)
    chunk = 5 * 10**7  # process 50M at a time

    for j, k in enumerate(major_idx):
        phase = 2.0 * pi * k / M
        S_val = 0j
        for c_start in range(0, len(lam), chunk):
            c_end = min(c_start + chunk, len(lam))
            ns = np.arange(c_start, c_end, dtype=np.float64)
            S_val += np.sum(lam[c_start:c_end] * np.exp(1j * phase * ns))
        S2_major[j] = S_val ** 2
        if (j+1) % 100 == 0:
            print(f"    {j+1}/{len(major_idx)} ({time.time()-t0:.0f}s)")

    print(f"  Done: {time.time()-t0:.1f}s\n")

    # Step 4: Correlate
    print(f"  [3/3] Correlation...")
    maj_c, min_c = [], []
    for i, Ni in enumerate(scan):
        phases = -2.0 * pi * Ni * major_idx / M
        major_re = float(np.sum((S2_major * np.exp(1j * phases)).real))
        maj_c.append(major_re)
        min_c.append(totals[i] - major_re)

    ma, mi = np.array(maj_c), np.array(min_c)
    r, p_val = stats.pearsonr(ma, mi)
    r_sp, _ = stats.spearmanr(ma, mi)
    frac_pos = float(np.mean(np.array(min_c) > 0))

    total_time = time.time() - t_start
    tag = '✓ REINFORCING' if r > 0.3 else '~ MODERATE' if r > 0.1 \
          else '○ WEAK' if r > -0.1 else '✗ OPPOSING'

    print(f"""
  ┌──────────────────────────────────────────────┐
  │  N = {N_center:>14}                           │
  │  Pearson  r = {r:>+8.4f}       {tag:>14}  │
  │  Spearman ρ = {r_sp:>+8.4f}                       │
  │  p-value    = {p_val:>10.2e}                   │
  │  n_scan     = {len(scan):>6}                       │
  │  frac(minor > 0) = {frac_pos:.3f}                    │
  │  Time: {total_time:.0f}s ({total_time/60:.0f}min)                     │
  └──────────────────────────────────────────────┘""")

    cs = hashlib.sha256(json.dumps(
        {'N':N_center,'r':round(r,6)}).encode()).hexdigest()[:16]
    print(f"  Checksum: {cs}\n")
    return {'N':N_center,'r':float(r),'rho':float(r_sp),
            'p':float(p_val),'frac_pos':frac_pos,'time':total_time}

if __name__ == "__main__":
    run_1e9()
