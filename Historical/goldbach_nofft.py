"""
Goldbach correlation at 10^9 WITHOUT FFT.
Direct computation. ~16 GB RAM. No swap needed.
"""
import numpy as np
from math import log, gcd, isqrt, pi
import time, hashlib, json, gc as GC
from scipy import stats

def von_mangoldt(N):
    try:
        import primesieve
        print(f"  Sieve to {N:.2e}...", end="", flush=True)
        t0 = time.time()
        lam = np.zeros(N+1, dtype=np.float64)
        for p in primesieve.primes(N):
            lp = log(p)
            pk = p
            while pk <= N: lam[pk] = lp; pk *= p
        print(f" {time.time()-t0:.1f}s")
        return lam
    except ImportError:
        print("Need primesieve for 10^9"); return None

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
    print(f"  NO-FFT CORRELATION: N = {N_center:.2e}")
    print(f"  RAM needed: ~16 GB (no swap required)")
    print(f"{'='*60}\n")

    # Step 1: Sieve
    print("  [1/4] Sieve...")
    lam = von_mangoldt(N_max)
    if lam is None: return None
    psi = float(np.sum(lam))
    print(f"    psi = {psi:.0f} (ratio: {psi/N_max:.4f})")
    print(f"    Array: {lam.nbytes/1e9:.1f} GB")

    # Step 2: Direct total at scan points
    # r_total(N) = sum_{m=1}^{N-1} lam[m]*lam[N-m]
    lo = N_center - window + ((N_center-window)%2)
    hi = N_center + window
    scan = [n for n in range(max(6,lo), hi+1, scan_step) if n%2==0]

    print(f"\n  [2/4] Direct totals for {len(scan)} N values...")
    print(f"    (each is O(N) ~ 10^9 multiplies)")
    t0 = time.time()
    
    totals = []
    for i, Ni in enumerate(scan):
        # lam[1:Ni] * lam[Ni-1:0:-1] but we need careful indexing
        # r(Ni) = sum_{m=1}^{Ni-1} lam[m]*lam[Ni-m]
        m_arr = np.arange(1, min(Ni, len(lam)))
        complement = Ni - m_arr
        valid = complement < len(lam)
        total = float(np.sum(lam[m_arr[valid]] * lam[complement[valid]]))
        totals.append(total)
        if (i+1) % 20 == 0 or i == 0:
            elapsed = time.time() - t0
            eta = elapsed/(i+1) * (len(scan)-i-1)
            print(f"      {i+1}/{len(scan)}"
                  f"  r_L({Ni})={total:.0f}"
                  f"  [{elapsed:.0f}s elapsed, ~{eta:.0f}s remaining]")

    t_totals = time.time() - t0
    print(f"    Done: {t_totals:.1f}s")

    # Step 3: Major-arc S(alpha)^2 at specific frequencies
    major_idx = get_major_indices(M, Q_max)
    print(f"\n  [3/4] Major-arc S(a)^2 at {len(major_idx)} frequencies...")
    print(f"    (each frequency: dot product over 10^9 terms)")
    t0 = time.time()

    # For each k in major_idx, compute:
    #   S(k/M) = sum_{n<=N_max} lam[n] * exp(2*pi*i*n*k/M)
    # Then S^2, then for each Ni: Re[S^2 * exp(-2*pi*i*Ni*k/M)]
    
    n_arr = np.arange(len(lam), dtype=np.float64)  # 8 GB
    
    # Process major frequencies in batches to control memory
    batch_size = 64  # 64 frequencies at a time
    S2_major = np.zeros(len(major_idx), dtype=np.complex128)
    
    for b_start in range(0, len(major_idx), batch_size):
        b_end = min(b_start + batch_size, len(major_idx))
        batch_k = major_idx[b_start:b_end]
        
        for j, k in enumerate(batch_k):
            # S(k/M) = sum lam[n] * exp(2*pi*i*n*k/M)
            phase = 2.0 * pi * k / M
            # Process in chunks to avoid 8GB complex temp array
            chunk = 10**7
            S_val = 0j
            for c_start in range(0, len(lam), chunk):
                c_end = min(c_start + chunk, len(lam))
                ns = n_arr[c_start:c_end]
                S_val += np.sum(lam[c_start:c_end] *
                               np.exp(1j * phase * ns))
            S2_major[b_start + j] = S_val ** 2
        
        if b_start % 128 == 0:
            print(f"      {b_start}/{len(major_idx)} frequencies"
                  f"  ({time.time()-t0:.0f}s)")

    del n_arr; GC.collect()
    t_major = time.time() - t0
    print(f"    Done: {t_major:.1f}s")

    # Step 4: Correlate
    print(f"\n  [4/4] Computing correlations...")
    
    maj_c, min_c = [], []
    for i, Ni in enumerate(scan):
        # Major contribution: sum Re[S2 * exp(-2*pi*i*Ni*k/M)]
        phases = -2.0 * pi * Ni * major_idx / M
        major_re = float(np.sum((S2_major * np.exp(1j * phases)).real))
        minor_re = totals[i] - major_re
        maj_c.append(major_re)
        min_c.append(minor_re)

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
  │  Time: {total_time:.0f}s                              │
  └──────────────────────────────────────────────┘""")

    cs = hashlib.sha256(json.dumps(
        {'N':N_center,'r':round(r,6)}).encode()).hexdigest()[:16]
    print(f"  Checksum: {cs}\n")
    return {'N':N_center,'r':float(r),'rho':float(r_sp),
            'p':float(p_val),'frac_pos':frac_pos,'time':total_time}

if __name__ == "__main__":
    run_1e9()
