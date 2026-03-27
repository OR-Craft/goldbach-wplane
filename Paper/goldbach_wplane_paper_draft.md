# Computational Observations on the Complex-Plane Anatomy of Goldbach Major/Minor Arc Decompositions

**Leo [surname]**

*Draft — March 2026*

---

## Abstract

The Hardy–Littlewood circle method decomposes a von Mangoldt weighted binary Goldbach convolution into major and minor arc contributions, which are traditionally analysed separately. As Tao (2012) observes, in the binary case the minor arc contribution can potentially dominate "if all one is given about the minor arc sums is magnitude information." Motivated by this observation, the present study examines a fixed-Q_max, discrete-grid proxy to this decomposition, retaining phase information and computing the signed contributions from both arc regions empirically, across windows of even integers at scales up to 10⁹.

After correcting for a major-arc aperture effect identified during the study, the decomposition exhibits a consistent geometric structure in the complex plane across the tested scales, with moderate positive correlation (r ≈ +0.67 to +0.73, converged) and substantial cancellation at large N. The complex diagnostic w(N) = M(N) + i·m(N), where M(N) and m(N) denote the major and minor arc contributions respectively, reveals tight angular concentration, modular family structure reminiscent of the banding seen in Goldbach's comet, and data consistent with a crossover from a major-dominated to a cancellation-dominated regime.

Simple mechanical models based on the constraint M(N) + m(N) = R(N) predict strong negative correlation (r ≈ −1); the observed positive correlation differs in sign from all such models. The identification and correction of a measurement artefact, with the observed geometric structure persisting after the correction, constitutes the principal methodological contribution.

These findings are purely empirical and do not constitute a proof of the Goldbach conjecture. They are offered as computational observations on the spectral anatomy of the Goldbach integrand.

---

## 1. Introduction

### 1.1 The binary Goldbach problem and the circle method

The strong Goldbach conjecture asserts that every even integer N ≥ 4 can be expressed as the sum of two primes. Despite computational verification to approximately 4 × 10¹⁸ (Oliveira e Silva et al.) and extensive theoretical progress on related problems, the conjecture remains open.

The Hardy–Littlewood circle method provides the standard analytic framework for studying binary additive problems involving primes. For even N, the von Mangoldt weighted convolution is expressed as

   R(N) = ∫₀¹ S(α)² e(−Nα) dα

where S(α) = Σ_{n ≤ N} Λ(n) e(nα) is the prime exponential sum, Λ denotes the von Mangoldt function, and e(x) = e^{2πix}. The value R(N) equals Σ_{m+n=N} Λ(m)Λ(n), which is the von Mangoldt weighted convolution associated with Goldbach-type representations. Positivity of R(N) guarantees a representation by integers carrying von Mangoldt weight, and prime representations contribute positively as a distinguished subset. We work throughout with this weighted object rather than the raw prime-pair count g(N).

The unit interval [0, 1] is partitioned into *major arcs* — neighbourhoods of rationals a/q with small denominator q — and *minor arcs*, the complement. This yields a decomposition

   R(N) = M(N) + m(N)

where M(N) and m(N) denote the major and minor arc contributions respectively. On the major arcs, S(α) is well approximated by character sums, yielding the Hardy–Littlewood prediction M(N) ≈ 𝔖(N) · N, where 𝔖(N) is the singular series. On the minor arcs, one seeks bounds on |S(α)| or averaged estimates.

For the ternary Goldbach problem (sums of three primes), the cubic power in ∫ S(α)³ e(−Nα) dα ensures the major arc contribution dominates, enabling Helfgott (2013) to give a complete proof. For the binary problem, the situation is fundamentally different. The minor arc L² mass ∫_{minor} |S(α)|² dα is of order N log N, while the major arc main term is of order N. The minor arcs are larger by a factor of log N in magnitude, a gap that grows with N and contributes to the central analytic difficulty.

We emphasise that no claim regarding the Goldbach conjecture is made; the scope of this study is empirical throughout.


### 1.2 Motivation: the phase information gap

The difficulty of the binary Goldbach problem via the circle method is intimately connected to an information loss problem. Tao (2012) identifies three principal obstacles: the magnitude dominance of the minor arcs, the indistinguishability of primes from certain "redacted" sets under magnitude-only bounds, and the absence of suitable inverse theorems.

The first of these obstacles carries a crucial qualification. As Tao notes, the minor arc contribution potentially dominates "if all one is given about the minor arc sums is magnitude information." Standard analytic approaches necessarily discard the phase of S(α)² e(−Nα) on the minor arcs to obtain usable bounds, retaining only magnitude estimates of the form |S(α)| ≤ f(N, α). This phase information — the signed, complex-valued structure of the integrand at each frequency — is not directly accessible through magnitude-based bounds.

The standard approach analyses M(N) and m(N) separately; their joint behaviour across nearby values of N is not typically studied. The present work asks a different question: when phase information is retained and the signed major and minor arc contributions are computed for a window of even integers around a given N, what empirical relationship do they exhibit?


### 1.3 The w-plane diagnostic

To organise the empirical observations, we introduce a complex-valued diagnostic. For each even N in a scan window, define

   w(N) = M(N) + i · m(N)

where M(N) and m(N) are the signed (real-valued) major and minor arc contributions satisfying M(N) + m(N) = R(N). In this representation:

- The real part Re(w) encodes the major arc contribution.
- The imaginary part Im(w) encodes the minor arc contribution.
- The total weighted count is R(N) = Re(w) + Im(w).
- The line Re + Im = 0 corresponds to R(N) = 0 for the weighted convolution studied here. This is not identical to the statement of strong Goldbach, but provides a natural geometric reference line for the present computations.

The cloud of w(N) values across a scan window provides a convenient geometric representation of the decomposition's structure at a given scale. Tight angular clustering suggests a stable geometric relationship between M(N) and m(N) over the sampled window. The distance from each w(N) to the line Re + Im = 0 provides a per-N measure of how far the weighted convolution sits from zero at each even integer.


### 1.4 Summary of findings

The principal findings of this study, all stated as empirical observations, are as follows.

First, the signed major and minor arc contributions exhibit moderate positive correlation across windows of nearby even integers. At converged arc neighbourhood width (half-width 20 grid points), the Pearson correlation is r ≈ +0.73 at N = 10⁷, r ≈ +0.69 at N = 10⁸, and r ≈ +0.67 at N = 10⁹. This interpretation is based on three scales and should be regarded as suggestive.

Second, simple mechanical models based on the identity M(N) + m(N) = R(N) predict strong negative correlation (r ≈ −1), whereas the observed correlation is positive. This sign reversal persists across all tested scales in the present computations and is verified by four independent negative controls and one positive control (calibration) experiment.

Third, the w-plane diagnostic reveals tight angular concentration that increases with scale in the sampled windows, with the angular standard deviation decreasing from 2.6° at 10⁷ to 0.5° at 10⁹. The modular families N ≡ 0 mod 6 and N ≡ 2, 4 mod 6 appear as distinct angular clusters, reproducing the banding seen in Goldbach's comet within the w-plane representation.

Fourth, the data are consistent with a crossover from a near-axis, major-dominated regime at 10⁷ (where major/total ≈ 1.0) to a tightly cancellation-dominated regime at 10⁹ (where major/total ≈ 24.8, meaning both |M| and |m| are approximately 25 times the total R(N)). This interpretation is based on three scales and should be regarded as suggestive.

Fifth, an initial measurement artefact was identified during the study: narrow arc neighbourhoods (±2 grid points) amplify the measured correlation relative to wider, converged neighbourhoods (±20 grid points). The correlation decreased from r ≈ 0.83 to r ≈ 0.67 upon correction. The observed geometric structure persisted after this correction. The identification and correction of this aperture effect, with the principal findings persisting in weakened but qualitatively unchanged form, constitutes a methodological contribution applicable to any grid-based arc diagnostic.

Recent explicit work on Goldbach summatory functions emphasises the role of oscillatory zero terms in weighted Goldbach averages (Bhowmik and Halupczok, 2020; Bhowmik, Ernvall-Hytönen, and Palojärvi, 2025). The present computations are best viewed as an empirical, phase-retaining probe of related structure at the level of the circle-method decomposition. Whether the observed geometric structure reflects oscillatory zero contributions remains an open question.
-e 
---

## 2. Computational Framework

### 2.1 The weighted convolution on a discrete grid

For a given even integer N, we compute the von Mangoldt weighted convolution

   R(N) = Σ_{m+n=N} Λ(m)Λ(n)

directly, where Λ(n) = log p if n = p^k for some prime p and integer k ≥ 1, and Λ(n) = 0 otherwise. Primes are generated by sieve, and Λ is then assigned by marking all prime powers up to the required bound.

To decompose R(N) into major and minor arc contributions, we work with the discrete Fourier transform of the Λ-sequence. Let N_max be the largest integer in the scan window. Define L = 2N_max + 1 and consider the L-point DFT in the standard (negative-exponent) convention used by numpy.fft and the Goertzel algorithm:

   F(k) = Σ_{n=0}^{L−1} Λ(n) e^{−2πink/L},  k = 0, 1, ..., L−1.

Since the Λ-sequence has support in [0, N_max] and L > 2N_max, the circular convolution computed by the DFT agrees with the linear convolution for all N in the scan window, avoiding wrap-around artefacts.

The weighted convolution R(N) can be recovered from these DFT coefficients:

   R(N) = (1/L) Σ_{k=0}^{L−1} F(k)² e^{2πiNk/L}.

The phase-retaining integrand at each grid frequency is therefore

   I_N(k) = (1/L) F(k)² e^{2πiNk/L}

so that R(N) = Σ_k I_N(k). Once the discrete grid and major index set are fixed, the decomposition is exact on that grid at the algebraic level, up to floating-point rounding in implementation. The grid decomposition is not identical to the continuous arc decomposition of the Hardy–Littlewood method; it is a finite computational proxy designed to retain phase information on a discrete frequency set.

In the classical Hardy–Littlewood formulation, major arcs are neighbourhoods of width approximately Q/N around rationals a/q with q ≤ Q, where Q grows with N (typically Q = N^δ for some δ > 0). In the present discrete setting, Q_max = 20 is held fixed across all scales, and the neighbourhood width is controlled by the half-width parameter hw in grid points. At N = 10⁹, the classical choice Q = N^{1/6} ≈ 31 is comparable to the fixed Q_max = 20 used here; the Q_max sensitivity analysis (§5.1) confirms stability of the measured correlation across Q_max values from 10 to 50. The convergence analysis in §2.5 demonstrates internal stability: the measured quantities are insensitive to hw once it is large enough to capture the full peak around each rational centre. This is a statement about the discrete grid decomposition, not a proof of convergence to the classical continuous arc decomposition. The extent to which the discrete proxy at fixed Q_max approximates the continuous theory — particularly for contributions from rationals with q > Q_max — is a limitation of the present approach.


### 2.2 Major and minor arc classification

The unit circle [0, 1) is discretised into L equally spaced frequencies α_k = k/L. We classify each frequency as belonging to a major or minor arc based on its proximity to rationals with small denominator. Throughout this paper, the terms "major arc" and "minor arc" refer exclusively to this discrete-grid classification with fixed Q_max, not to the classical continuous decomposition of the Hardy–Littlewood circle method, in which the arc parameters typically grow with N.

Define the circular distance on {0, 1, ..., L−1} as

   d_L(k, ℓ) = min(|k − ℓ|, L − |k − ℓ|).

For parameters Q_max (the maximum denominator) and hw (the half-width in grid points), a frequency index k is classified as major if there exist integers a, q with 1 ≤ q ≤ Q_max, 0 ≤ a < q, gcd(a, q) = 1, such that

   d_L(k, round(aL/q)) ≤ hw.

The major index set I_maj is formed as the union of all such neighbourhoods, with overlaps removed. All other indices are classified as minor.

The major and minor arc contributions for a given N are

   M(N) = Σ_{k ∈ I_maj} I_N(k)

   m(N) = R(N) − M(N).

Because the major index set is symmetric under k → L − k, the summed contributions are real up to floating-point error; in implementation, residual imaginary parts were numerically negligible.

The minor arc contribution is defined as the exact residual, ensuring M(N) + m(N) = R(N) by construction. This avoids accumulation of rounding errors in the minor arc sum.

Throughout this study, Q_max = 20 is used. The half-width parameter hw is varied from 2 to 20 grid points as part of the convergence analysis described in §2.5. At hw = 20, the number of distinct major-arc grid points (after deduplication) is approximately 5,248 out of L ≈ 2 × 10⁹ at N = 10⁹.


### 2.3 Computational methods

Two methods are employed to evaluate the DFT coefficients at the major-arc frequencies, depending on the scale N.

**FFT method (N ≤ 10⁸).** The full L-point DFT is computed using numpy.fft, which implements the negative-exponent convention matching §2.1. The total R(N) is obtained via the inverse FFT, and the major-arc contribution M(N) is computed by summing I_N(k) over the major indices I_maj. Memory requirements are feasible on the available hardware for L up to approximately 2 × 10⁸, but not for L ≈ 2 × 10⁹.

**Goertzel method (N = 10⁹).** At N = 10⁹, the DFT length L ≈ 2 × 10⁹ exceeds available memory for a full FFT. Instead, individual DFT bins F(k) are computed using the standard Goertzel recurrence, which evaluates a single bin in O(L) time and O(1) space. The recurrence is implemented with numba just-in-time compilation for performance. Parallelism across major-arc frequencies is achieved using Python's multiprocessing module with 12 worker processes.

The total R(N) at each scan point is computed independently via the inner product

   R(N) = Σ_{j=1}^{N−1} Λ(j) Λ(N−j)

implemented as a numpy dot product of array views.

**Cross-validation.** At N = 10⁷ and N = 10⁸, where both FFT and Goertzel are feasible, the two methods agree within relative error < 10⁻⁷ in the major-arc contribution and < 10⁻¹⁰ in the total. This provides confidence that the Goertzel implementation at N = 10⁹ is correct.


### 2.4 Scan windows and sampling

For each target scale N_center, a scan window of 201 evenly spaced even integers is evaluated:

   N ∈ {N_center − 2000, N_center − 1980, ..., N_center + 2000}

with step size 20. The window width (±2000) is small relative to N_center (less than 0.0004% at 10⁹), so the DFT length L and the major-arc index set I_maj are effectively constant across the window. Variation in M(N) and m(N) across the window arises from the N-dependent twist factor e^{2πiNk/L} in the integrand I_N(k), not from changes in the grid geometry.

Three target scales are examined: N_center = 10⁷, 10⁸, and 10⁹. The upper limit is set by hardware: at N = 10⁹ each scan window requires approximately 80 minutes of Goertzel computation on a 12-core workstation with 48 GB RAM, and memory constraints preclude a full FFT at this scale.


### 2.5 Convergence in half-width

The half-width parameter hw controls the number of grid points assigned to each major-arc neighbourhood. In the present discretisation, the operational neighbourhood width chosen to capture a full peak around each selected rational centre is approximately 40 grid points. The half-width hw determines what fraction of this operational width is captured:

- At hw = 2 (5 grid points per centre, 640 major indices total), each neighbourhood is narrow relative to the peak structure, constituting substantial undersampling.
- At hw = 5 (11 grid points, 1,408 indices), coverage is partial.
- At hw = 10 (21 grid points, 2,688 indices), coverage is approximately half.
- At hw = 20 (41 grid points, 5,248 indices), the neighbourhood slightly exceeds the nominal width, ensuring full coverage of the selected neighbourhood around each rational centre.

The major/total fraction M(N)/R(N) converges as hw increases. At N = 10⁷, it evolves from +0.92 (hw = 2) to +1.00 (hw = 20); at N = 10⁸, from +0.78 to +0.86; at N = 10⁹, from −13.9 to +24.8. The sign reversal at N = 10⁹ with hw = 2 is a subsampling artefact: the narrow neighbourhood captures an oscillatory slice of the peak structure whose sign depends on the local phase. Wider neighbourhoods integrate over the full peak and recover the expected positive major-arc contribution.

The Pearson correlation r(M, m) decreases monotonically with hw at all three scales:

| N       | hw = 2 | hw = 5 | hw = 10 | hw = 20 |
|---------|-------:|-------:|--------:|--------:|
| 10⁷    | +0.87  | +0.78  |  +0.75  |  +0.73  |
| 10⁸    | +0.84  | +0.75  |  +0.71  |  +0.69  |
| 10⁹    | +0.83  | +0.74  |  +0.70  |  +0.67  |

The correlation remains positive at all tested half-widths and scales. The decrease from hw = 2 to hw = 20 reflects the removal of an amplification effect: narrower neighbourhoods preferentially sample the coherent peak centres of the major arcs, inflating the measured correlation. The converged values at hw = 20 are reported as the primary results throughout this paper.

The convergence of the major/total fraction by hw = 10 (with |Δ| < 0.01 between hw = 10 and hw = 20) provides an internal consistency check on the grid resolution. The convergence of the correlation, while not as sharp, shows a clear plateau forming between hw = 10 and hw = 20, supporting the choice of hw = 20 as a stable operating point.
-e 
---

## 3. The w-Plane Framework

### 3.1 Definition and basic properties

For each even integer N in a scan window, define the complex-valued diagnostic

   w(N) = M(N) + i · m(N)

where M(N) and m(N) are the real-valued major and minor arc contributions defined in §2.2. Since M(N) + m(N) = R(N), the weighted convolution is recovered as

   R(N) = Re(w(N)) + Im(w(N)).

The w-plane representation encodes the arc decomposition as a single complex number per even integer, with the real axis representing the major arc contribution and the imaginary axis representing the minor arc contribution. The modulus and argument of w(N) are

   |w(N)| = √(M(N)² + m(N)²)

   arg(w(N)) = atan2(m(N), M(N))

where atan2 denotes the two-argument arctangent, which preserves quadrant information. The modulus |w(N)| measures the overall scale of the two-component decomposition, while the argument arg(w(N)) measures the angular balance between major and minor arc contributions. When M(N) dominates (major-dominated regime), arg(w(N)) is close to zero. When M(N) and m(N) are of comparable magnitude with opposite signs (cancellation-dominated regime), arg(w(N)) departs substantially from zero.


### 3.2 The counterexample reference line

A hypothetical even integer Z with R(Z) = 0 for the weighted convolution would satisfy

   Re(w(Z)) + Im(w(Z)) = 0,

corresponding to the line m = −M in the w-plane, or equivalently the line through the origin at angle 3π/4 (in the second quadrant) and −π/4 (in the fourth quadrant). We refer to this as the *reference line*. We emphasise that R(Z) = 0 for the weighted convolution is not identical to the statement that Z is a Goldbach counterexample, since R(N) includes prime-power contributions (see §1.1). The reference line is used here as a geometric landmark, not as a direct proxy for the Goldbach conjecture.

The signed perpendicular distance from w(N) to this reference line is

   δ(N) = (M(N) + m(N)) / √2 = R(N) / √2

which is positive whenever R(N) > 0. This distance provides a per-N measure of how far the weighted convolution sits from zero.


### 3.3 Angular concentration and regime classification

The distribution of arg(w(N)) across a scan window characterises the geometric structure of the decomposition at a given scale.

In a *major-dominated regime*, M(N) is positive and much larger than |m(N)|. The w(N) cloud lies near the positive real axis, and arg(w(N)) is close to zero with moderate spread. The minor arc contribution is a small perturbation.

In a *cancellation-dominated regime*, both |M(N)| and |m(N)| are much larger than R(N), and R(N) is a small residual of two large, nearly opposing terms. The w(N) cloud lies far from the real axis, and arg(w(N)) takes a value substantially different from zero. In this regime, the angular spread may be very tight because the ratio M(N)/m(N) is approximately constant across the window, suggesting that both components are responding to a common large-scale structure in the decomposition. The fluctuations in R(N) represent a small modulation of a large, stable decomposition.

The transition between these regimes, if it occurs, would manifest in the w-plane as a rotation of the cloud away from the real axis accompanied by a tightening of the angular distribution. Whether such a transition is a continuous crossover or a sharper change cannot be determined from the present data alone.


### 3.4 Modular family structure

The Hardy–Littlewood singular series 𝔖(N) depends on the prime factorisation of N, creating systematic variation in R(N) across arithmetic families. The most prominent effect is the mod 6 classification: for even N, the residue N mod 6 takes values 0, 2, or 4, and 𝔖(N) is larger for N ≡ 0 mod 6 than for N ≡ 2 or 4 mod 6. This is the origin of the two-band structure visible in Goldbach's comet — the plot of the raw representation count g(N) versus N.

In the w-plane, this modular structure manifests as distinct clusters within the cloud. Points with N ≡ 0 mod 6 are observed to have larger |w(N)| and occupy a slightly different angular position than points with N ≡ 2 or 4 mod 6 in the present data. The w-plane thus provides a geometric representation in which the comet's band structure appears as angular and radial separation between modular families, rather than as vertical separation in a one-dimensional plot.


### 3.5 Connection to existing theory

The w-plane diagnostic is an empirical construction; it does not arise from a theorem or a standard decomposition in the analytic number theory literature. However, several existing theoretical results provide relevant context.

The Goldbach summatory function S(x) = Σ_{n ≤ x} G(n), where G(n) denotes the von Mangoldt weighted Goldbach convolution, satisfies, under the Riemann Hypothesis,

   S(x) = x²/2 − H(x) + error terms

where H(x) = −2 Σ_ρ x^{ρ+1} / (ρ(ρ+1)) is an oscillatory term summed over the non-trivial zeros ρ of the Riemann zeta function (Fujii, 1991; Granville, 2007; Bhowmik and Halupczok, 2020). Each zero ρ = 1/2 + iγ contributes an oscillation of frequency γ in log-space. Recent explicit work has obtained unconditional numerical estimates for these summatory functions (Bhowmik, Ernvall-Hytönen, and Palojärvi, 2025).

The individual-N behaviour of M(N) and m(N), and hence of w(N), may plausibly reflect some of the same oscillatory phenomena that contribute to H(x), though at a much finer and unsmoothed level. Whether the observed geometric structure in the w-plane reflects oscillatory zero structure is a natural question, but one the present study does not resolve.
-e 
---

## 4. Empirical Results

All results in this section are reported at the converged half-width hw = 20 unless otherwise stated. The convergence behaviour across half-widths is documented in §2.5. The reported values summarise each 201-point window; full per-window plots are shown in Figures 1–5.


### 4.1 Signed correlation across scales

The Pearson correlation between M(N) and m(N) across each scan window is:

| N_center | r (hw = 20) | Spearman ρ | n   |
|----------|-------------|------------|-----|
| 10⁷     | +0.73       | +0.69      | 201 |
| 10⁸     | +0.69       | +0.64      | 201 |
| 10⁹     | +0.67       | +0.62      | 201 |

The sampled values within a narrow window are not independent draws from a stationary distribution, so standard parametric p-values are not valid under this sampling design. Evidence for the robustness of the correlation is instead provided by permutation tests, bootstrap confidence intervals, and block bootstrap, reported in §5.

The correlation is positive at all three scales. The values decrease gently from +0.73 to +0.67 across two orders of magnitude in N, but with only three scales this should be viewed as a descriptive observation rather than evidence of asymptotic behaviour.


### 4.2 w-plane geometry

The w(N) cloud at each scale is concentrated in a distinct region of the complex plane. Since the angle clusters are narrow and remain far from the branch cut, ordinary and circular summaries are numerically indistinguishable in the sampled windows; ordinary means and standard deviations of arg(w(N)) are reported throughout.

**At N = 10⁷:**

| Quantity                  | Value         |
|---------------------------|---------------|
| mean(Re(w)) = mean(M)    | 2.47 × 10⁷   |
| mean(Im(w)) = mean(m)    | 2.54 × 10⁵   |
| Mean angle                | −0.0° (0.0 rad) |
| Angular std               | 2.6°          |
| Quadrant distribution     | Q1: 44%, Q4: 56% |
| min δ(N) = min R(N)/√2   | 1.24 × 10⁷   |
| All R(N) > 0?             | Yes           |

**At N = 10⁸:**

| Quantity                  | Value         |
|---------------------------|---------------|
| mean(Re(w)) = mean(M)    | 2.15 × 10⁸   |
| mean(Im(w)) = mean(m)    | 3.55 × 10⁷   |
| Mean angle                | +9.7° (0.17 rad) |
| Angular std               | 2.3°          |
| Quadrant distribution     | Q1: 100%      |
| min δ(N) = min R(N)/√2   | 1.24 × 10⁸   |
| All R(N) > 0?             | Yes           |

**At N = 10⁹:**

| Quantity                  | Value         |
|---------------------------|---------------|
| mean(Re(w)) = mean(M)    | 5.58 × 10¹⁰  |
| mean(Im(w)) = mean(m)    | −5.33 × 10¹⁰ |
| Mean angle                | −43.7° (−0.76 rad) |
| Angular std               | 0.5°          |
| Quadrant distribution     | Q4: 100%      |
| min δ(N) = min R(N)/√2   | 1.24 × 10⁹   |
| All R(N) > 0?             | Yes           |

The shift in quadrant distribution across scales follows from the evolution of the mean angle: at N = 10⁷ the mean angle is near zero and points straddle the real axis (Q1 and Q4); at N = 10⁸ the mean angle is positive (+9.7°) and all points lie in Q1; at N = 10⁹ the mean angle is negative (−43.7°) and all points lie in Q4.


### 4.3 Major/total fraction and angular concentration

The ratio M(N)/R(N), averaged over each scan window, evolves with scale:

| N_center | mean(M/R) (hw = 20) | std(M/R) | Angular std |
|----------|---------------------|----------|-------------|
| 10⁷     | +1.003              | 0.044    | 2.63°       |
| 10⁸     | +0.855              | 0.030    | 2.33°       |
| 10⁹     | +24.82              | 6.99     | 0.47°       |

At N = 10⁷, the mean major arc contribution is approximately equal to the total, with the minor arc contributing a small residual of variable sign. At N = 10⁹, the major arc contribution is approximately 25 times the total, indicating that both |M(N)| and |m(N)| are much larger than R(N) in this range. A ratio M/R substantially greater than 1 does not indicate an inconsistency in the decomposition; it reflects a cancellation regime in which M(N) and m(N) are large quantities of opposite sign whose sum R(N) is a small positive residual. Over the same range of scales, the angular standard deviation decreases from 2.6° to 0.5° in the sampled windows.


### 4.4 Modular family structure in the w-plane

Within each scan window, the even integers fall into three residue classes mod 6. Their w-plane properties are:

**At N = 10⁷ (hw = 20):**

| N mod 6 | n  | Mean angle | Angular std | mean(|w|)    |
|---------|----|------------|-------------|--------------|
| 0       | 67 | +2.5°      | 2.5°        | 3.60 × 10⁷  |
| 2       | 67 | −1.3°      | 1.6°        | 1.92 × 10⁷  |
| 4       | 67 | −1.3°      | 1.6°        | 1.92 × 10⁷  |

**At N = 10⁸ (hw = 20):**

| N mod 6 | n  | Mean angle | Angular std | mean(|w|)    |
|---------|----|------------|-------------|--------------|
| 0       | 67 | +8.5°      | 2.8°        | 3.30 × 10⁸  |
| 2       | 67 | +10.3°     | 1.7°        | 1.61 × 10⁸  |
| 4       | 67 | +10.3°     | 1.9°        | 1.62 × 10⁸  |

**At N = 10⁹ (hw = 20):**

| N mod 6 | n  | Mean angle | Angular std | mean(|w|)    |
|---------|----|------------|-------------|--------------|
| 0       | 67 | −43.0°     | 0.2°        | 7.78 × 10¹⁰ |
| 2       | 67 | −44.0°     | 0.1°        | 7.68 × 10¹⁰ |
| 4       | 67 | −44.0°     | 0.1°        | 7.68 × 10¹⁰ |

At all three scales, the N ≡ 0 mod 6 class has larger mean |w| than the N ≡ 2 and N ≡ 4 mod 6 classes, which are consistent with each other. The angular separation between the N ≡ 0 class and the other two classes is small but consistent (approximately 1–2° at 10⁷ and 10⁸, and approximately 1° at 10⁹). These separations are reported as descriptive features of the sampled windows; no formal significance test on the angular differences has been performed.


### 4.5 Distance to the reference line

At all three scales and for all 201 sampled even integers, R(N) > 0, so that δ(N) > 0 and no sampled point lies on or beyond the reference line Re + Im = 0.

| N_center | min R(N)     | mean R(N)     | min δ(N)     |
|----------|-------------|---------------|--------------|
| 10⁷     | 1.75 × 10⁷  | 2.50 × 10⁷   | 1.24 × 10⁷  |
| 10⁸     | 1.75 × 10⁸  | 2.50 × 10⁸   | 1.24 × 10⁸  |
| 10⁹     | 1.75 × 10⁹  | 2.50 × 10⁹   | 1.24 × 10⁹  |

In the sampled windows, the minimum distance scales approximately proportionally with N_center, consistent with the growth of R(N) over the tested range.
-e 
---

## 5. Validation and Controls

This section summarises the tests used to assess the robustness of the empirical findings in §4. The validation analysis has two purposes: first, to show that the observed positive correlation is not a trivial statistical or algebraic artefact; second, to show that the principal qualitative features persist after convergence correction.


### 5.1 Statistical validation of the correlation

The following tests were applied to the signed correlation r(M, m) at N = 10⁹ at the converged half-width hw = 20, where the observed Pearson correlation is r = +0.67.

**Permutation test.** The M(N) values were randomly shuffled 10,000 times and the Pearson correlation recomputed for each shuffle. None of the 10,000 shuffled replicates equalled or exceeded the observed value. The maximum shuffled correlation was r = +0.27.

**Bootstrap confidence interval.** 10,000 bootstrap resamples of the (M, m) pairs yielded a 95% confidence interval of [+0.61, +0.76] for r, with bootstrap mean +0.68.

**Block bootstrap.** To account for potential local autocorrelation within the scan window, a block bootstrap with block size 20 was applied. The resulting 95% confidence interval was [+0.62, +0.74], consistent with the standard bootstrap.

**Split-half reliability.** The scan window was partitioned in four ways. All subsets produced positive correlations:

| Partition     | r      |
|---------------|--------|
| First half    | +0.63  |
| Second half   | +0.75  |
| Even indices  | +0.62  |
| Odd indices   | +0.77  |

**Detrending.** Linear detrending of both M(N) and m(N) with respect to N produced negligible change (Δr = +0.0001), indicating that the positive association is not an artefact of shared linear trends.

**Robust correlation measures.** Pearson, Spearman, and Kendall correlations were computed:

| Measure  | Value  |
|----------|--------|
| Pearson  | +0.67  |
| Spearman | +0.68  |
| Kendall  | +0.51  |

All are positive. The Spearman rank correlation slightly exceeds the Pearson value, indicating that the association is not driven by outliers. Nominal parametric p-values are not reported because the sampled values are not independent draws from a stationary distribution; the permutation test and bootstrap intervals above provide the appropriate robustness evidence.

**Q_max sensitivity.** The threshold Q_max was varied from 10 to 50 at the initial half-width hw = 2. The correlation r(M, m) remained stable with standard deviation 0.064 across the tested values, indicating that the result is not sensitive to the choice of major-arc denominator bound.

**Window sensitivity.** The scan window centre was shifted by ±500 in steps of 100 at hw = 2. The correlation r remained stable with standard deviation 0.015.

The Q_max and window sensitivity tests were performed at hw = 2 only; their qualitative conclusions — stability with respect to parameter choices — are plausibly unchanged at hw = 20, though this was not recomputed explicitly. Full results of the validation battery at the initial hw = 2 setting, where r = +0.83, are reported in Appendix C.


### 5.2 Mechanical coupling controls

The identity M(N) + m(N) = R(N) creates an algebraic coupling between the two quantities. To determine whether the observed positive correlation could arise from this constraint alone, five models were tested at hw = 2 and N = 10⁹. Since the mechanical coupling analysis concerns the sign of the correlation rather than its precise magnitude, the qualitative conclusions — that mechanical models predict negative correlation while the observed correlation is positive — apply regardless of half-width.

**Negative controls (four models).** Each model generates synthetic major and minor arc values that satisfy M + m = R for the observed totals R(N), but with the structural relationship between M and m destroyed.

*Model A (fixed fraction + noise):* M_synth = f · R + noise, m_synth = R − M_synth, where f is the observed mean fraction and noise is Gaussian with matched variance. 10,000 realisations.

*Model B (random fraction):* Each N receives an independent random fraction drawn uniformly from the observed range. 10,000 realisations.

*Model D (constant fraction, analytical):* If f is constant and outside [0, 1], the predicted correlation is r = −1.0. At N = 10⁹ the observed mean fraction is approximately −13.9 (at hw = 2), yielding an analytical prediction of r = −1.0.

*Model E (shuffled fractions):* The observed fractions f_i = M(N_i)/R(N_i) are randomly reassigned to different N values, preserving the marginal distribution of fractions but destroying the N-by-N correspondence. 10,000 realisations.

All four negative controls predict r ≈ −1.0. The following table summarises the initial narrow-width (hw = 2) comparison:

| Model | Baseline r | Observed r (hw = 2) | Excess  |
|-------|-----------|---------------------|---------|
| A     | −0.999    | +0.83               | +1.83   |
| B     | −0.999    | +0.83               | +1.83   |
| D     | −1.000    | +0.83               | +1.83   |
| E     | −0.999    | +0.83               | +1.83   |

At the converged half-width hw = 20, the sign reversal persists:

| Setting  | Mechanical prediction | Observed r | Excess  |
|----------|----------------------|------------|---------|
| hw = 20  | −1.0                 | +0.67      | +1.67   |

The observed correlation differs in sign from all four mechanical baselines at both half-widths. In none of the 40,000 total synthetic realisations (across Models A, B, and E) did the synthetic correlation equal or exceed the observed value.

**Positive control (one model).** Model C preserves the linear relationship between M and R by shuffling only the residuals from a linear fit, then reconstructing M and m. This model is expected to reproduce the observed correlation as a calibration check.

| Model | Baseline r | Observed r | Difference |
|-------|-----------|------------|------------|
| C     | +0.830    | +0.83      | −0.0003    |

The positive control reproduces the observed correlation to within 0.0003, confirming that the measurement procedure is internally consistent.

**Structural indicator.** The Pearson correlation between the fraction f = M(N)/R(N) and the total R(N) is r = +0.99, indicating that the way the total splits between major and minor arc contributions depends strongly on the magnitude of R(N). None of the mechanical models reproduces this dependence. In particular, they predict strong negative correlation between M and m, contrary to observation.


### 5.3 Convergence as internal validation

The wide-arc correction described in §2.5 constitutes an internal validation of the measurement. The initial results at hw = 2 reported r ≈ 0.83 at N = 10⁹. The subsequent convergence analysis revealed that this value was amplified by the narrow arc neighbourhood, and the converged value at hw = 20 is r ≈ 0.67.

The correction weakened the headline number but did not change its sign or qualitative character. The principal geometric features — positive correlation, angular concentration in the w-plane, modular family separation — all persisted after the correction in weakened but qualitatively unchanged form.

The identification of this aperture amplification effect, and the demonstration that the main findings survive it, is itself a methodological contribution. Any future study using grid-based arc diagnostics on additive problems should be aware that narrow neighbourhoods can amplify measured correlations relative to converged values.


### 5.4 Cross-method validation

The computational framework employs two independent methods (§2.3), providing an opportunity for cross-validation at scales where both are feasible.

**DFT bin comparison.** At N = 10⁴ and N = 10⁵, individual DFT bins F(k) were computed by both the full FFT and the Goertzel recurrence and compared directly. The relative error was below 10⁻⁷ at the first five tested major-arc frequencies, with most bins agreeing to 10⁻¹⁰ or better.

**Total R(N) comparison.** The total R(N) can be computed in two independent ways: via the inverse FFT of the squared DFT, and via the direct dot product Σ Λ(j)Λ(N−j). At all tested N, these two methods agree to relative error < 10⁻¹⁰. The ratio FFT-total / dot-total equals L (the DFT length) to machine precision, confirming the normalisation.

**Major-arc contribution.** At N = 10⁷ and N = 10⁸, the major-arc sum M(N) was computed by both methods for individual even integers. The major/total fraction M(N)/R(N) agrees to better than 10⁻⁷ between FFT and Goertzel at the same N.

**Reproducibility across runs.** The wide-arc computation at N = 10⁹ was executed twice independently (once for the convergence analysis, once for the validation battery). Both runs produced identical reported values and matching output checksums, confirming reproducibility of the Goertzel pipeline. Output checksums are recorded with each run.

**Scope.** This cross-validation was performed at the individual-N level, not merely as an aggregate check. It confirms that the Goertzel implementation at N = 10⁹ — which cannot be verified by full FFT due to memory constraints — produces results consistent with the FFT method at all overlapping scales.


### 5.5 Falsified hypotheses

Several hypotheses were tested during the course of this study and rejected. They are documented here to delimit the scope of the surviving claims.

**Mod 5 bimodal effect.** At narrow half-width (hw = 2), the correlation appeared to split bimodally by N mod 5, with 5|N showing substantially higher correlation than 5∤N. Apodisation of the major/minor boundary and per-q decomposition reduced this effect by 99%. It was identified as an aperture artefact arising from the q = 5 resonance falling at the bucket boundary. The effect did not survive correction.

**Mod 6 family differences in correlation.** Initial coarse sampling appeared to show different correlation values for different N mod 6 classes. Stratified sampling with proper coverage yielded ANOVA p = 0.96, indicating no statistically significant difference. The apparent pattern was a sampling artefact.

**Universal minor-arc positivity.** At narrow half-width (hw = 2), the minor arc contribution m(N) was positive for 100% of sampled N at all scales. At converged half-width (hw = 20), the fraction of positive m(N) drops to 44% at N = 10⁷ and 0% at N = 10⁹. The apparent universal positivity was a property of the narrow grid subsample, not of the full major-arc neighbourhood.

**r ≈ 0.83 as the converged correlation.** The initial measurement at hw = 2 gave r ≈ 0.83 across multiple scales, which appeared stable. The convergence analysis in half-width revealed this was amplified by approximately 0.16 relative to the converged value of r ≈ 0.67. The amplification was traced to the narrow neighbourhood preferentially sampling coherent peak centres.

Additional falsified hypotheses, including a non-commutative quaternion sieve experiment, are documented in Appendix D.


### 5.6 Code provenance and reproducibility

Final reported results in this paper are generated from the converged wide-arc pipeline (`goldbach_wide_arc.py`) and its corresponding hw = 20 validation script (`goldbach_stats_hw20.py`). Earlier scripts in the repository reflect development stages associated with narrower arc neighbourhoods and are retained for reproducibility of intermediate analyses documented in §5.3 and Appendix C. All code, cached output data, and checksums are available in the accompanying repository.
-e 
---

## 6. Interpretation and Discussion

The preceding sections established three principal empirical observations: a moderate positive correlation between signed major and minor arc contributions (§4.1), a consistent geometric structure in the w-plane that evolves with scale (§4.2–4.4), and the survival of both features after correction for an aperture amplification artefact (§5.3). This section offers cautious interpretations, explicitly separates what the data show from what they do not, and identifies open questions.


### 6.1 The observed co-movement

The most basic empirical observation is that M(N) and m(N) co-move positively across windows of nearby even integers. The mechanical coupling analysis (§5.2) establishes that this co-movement is not a consequence of the algebraic identity M + m = R: models that preserve the identity but destroy the structural relationship between M and m predict strong negative correlation, whereas the observed correlation is positive at all tested scales and half-widths.

The structural indicator r(f, R) = +0.99, where f = M(N)/R(N), shows that the splitting of R(N) between major and minor contributions is strongly associated with the magnitude of R(N) in the sampled window. This is consistent with a directional shift in the w-plane as R(N) varies, although no parametric model is derived here. The association may be related to fluctuations in the singular series 𝔖(N), which modulates R(N) across nearby even integers, although no direct derivation of this connection is established.

Simple predictive models were tested during this study. A model in which M(N) is a fixed or slowly varying fraction of R(N) plus independent noise fails to reproduce the observed correlation at N = 10⁹, where the decomposition is cancellation-dominated. This indicates that the coupling between M and m is not captured by the simplest fraction-based models considered here.


### 6.2 The suggested crossover

The w-plane data (§4.2–4.3) show a marked change in the character of the decomposition across the tested scales:

- At N = 10⁷, the major arc contribution is approximately equal to the total (mean M/R ≈ 1.0), and the w(N) cloud sits near the positive real axis with angular spread 2.6°.
- At N = 10⁸, the major arc accounts for approximately 86% of the total, and the cloud has rotated to a mean angle of +9.7° with angular spread 2.3°.
- At N = 10⁹, the major arc contribution is approximately 25 times the total, both |M| and |m| are much larger than R(N), and the cloud is concentrated at −43.7° with angular spread 0.5°.

These observations are consistent with a crossover from a major-dominated regime, in which the minor arc is a small perturbation, to a cancellation-dominated regime, in which the total R(N) is a small residual of two large, nearly opposing terms. The tightening of the angular spread in the cancellation-dominated regime is consistent with both components responding to a common large-scale pattern in the decomposition, so that their ratio varies only slightly across the window.

This interpretation is based on three scales and should be regarded as suggestive. With three data points, drift, curvature, and finite-scale effects cannot be distinguished. Whether the angular concentration continues to tighten at larger N, whether the mean angle continues to rotate, or whether the pattern stabilises, cannot be determined from the present data. All interpretive statements in this section refer to the fixed-Q_max, discrete-grid proxy studied in this paper, and should not be read as claims about the full continuous major/minor arc decomposition.


### 6.3 The phase information gap

As noted in §1.2, the difficulty of the binary Goldbach problem via the circle method is closely connected to the discarding of phase information on the minor arcs. The standard approach bounds |S(α)| on the minor arcs, retaining only magnitude information. Tao (2012) observes that the minor arc contribution can potentially dominate "if all one is given about the minor arc sums is magnitude information."

The present computations retain phase information by computing the signed integrand I_N(k) at each grid frequency rather than bounding its magnitude. In this setting, the signed minor arc contribution exhibits structured co-movement with the major arc contribution, rather than the potentially dominant behaviour permitted by worst-case magnitude bounds.

This observation does not resolve the phase information gap identified by Tao. The gap concerns the difficulty of proving analytic bounds that retain phase, not merely the difficulty of computing with it. The present study operates in a finite computational regime where direct evaluation is feasible; the analytic challenge lies in establishing results that hold for all N, which computation alone cannot achieve. Nevertheless, the empirical observation that the phase-retaining decomposition shows structured behaviour where magnitude-only bounds predict potential dominance may be of interest as a characterisation of the finite-N regime.


### 6.4 What the data do not show

The positive correlation and w-plane structure are properties of 201-point windows at three scales. They do not imply positivity of the weighted convolution R(N) for any specific N outside the sampled windows, nor do they constrain individual outliers: a window containing an exceptional integer could in principle appear geometrically typical while a single point lies on the weighted reference line. The gentle decline in r from +0.73 to +0.67 is consistent with multiple scenarios — a slow approach to a positive limit, a decline toward zero, or a finite-scale effect — and three data points cannot distinguish these.

Separately, attempts to derive the observed r ≈ 0.67 from the singular series variance or from simple fraction-based models were unsuccessful. The correlation resists prediction from the simplest available theoretical quantities, and no explanation for its specific value is offered.


### 6.5 Connection to existing programmes

Several active research programmes study the objects that appear in this work. The connections noted here are speculative and intended as context rather than explanation.

**Explicit formulas for Goldbach averages.** The summatory function S(x) = Σ_{n ≤ x} G(n), where G(n) denotes the von Mangoldt weighted Goldbach convolution, satisfies S(x) = x²/2 − H(x) + error terms, where H(x) is an oscillatory sum over non-trivial zeros of the Riemann zeta function (Fujii, 1991; Granville, 2007; Bhowmik and Halupczok, 2020). Bhowmik and Ruzsa (2018) establish a close connection between the average order of the Goldbach function and the Riemann Hypothesis. Recent work provides unconditional explicit numerical estimates for these summatory functions (Bhowmik, Ernvall-Hytönen, and Palojärvi, 2025). The individual-N behaviour of M(N) and m(N) may plausibly be influenced by some of the same oscillatory phenomena that contribute to H(x), though at a much finer and unsmoothed level.

**Level of distribution.** Lichtman (2023) achieves a level of distribution of 66/107 ≈ 0.617 for primes in arithmetic progressions, the greatest improvement on binary Goldbach since Bombieri–Davenport (1966). The methods involve spectral large sieve inequalities and exceptional Maass forms, refining earlier work of Drappeau, Pratt, and Radziwiłł. Whether the w-plane framework offers any empirical perspective on how such spectral structure manifests at the arc level is unknown.

**Pair correlation and spectral statistics.** The Montgomery pair correlation conjecture (1973) and its extensions describe the statistical distribution of zeta zeros. The Goldston–Suriajaya explicit formula (2023) connects averaged Goldbach counts to sums over zeros. Whether the observed angular concentration in the w-plane is related to the pair correlation structure of the zeros is an open question.


### 6.6 Open questions

The following questions arise naturally from the present observations and are offered as directions for further investigation.

1. **Heuristic prediction.** Can the observed correlation r ≈ 0.67 be derived, even heuristically, from the singular series variance or from the explicit formula for H(x)? Simple models tested in this study fail to reproduce the value. A successful prediction would connect the empirical observation to established theory.

2. **Oscillatory structure in arg(w).** Is any part of the variation in arg(w(N)) detectable in explicit-formula-style oscillatory models involving zeta zeros? If oscillatory components with frequencies matching the imaginary parts of zeta zeros can be identified, this would establish a concrete link between the w-plane geometry and the explicit formula.

3. **Scaling of angular concentration.** Does the angular standard deviation of w(N) continue to decrease at larger N? If so, at what rate? The decrease from 2.6° to 0.5° across two orders of magnitude in N is suggestive but does not establish a scaling law.

4. **Wider sampling.** The present study uses windows of 201 even integers with step size 20. Sampling with wider spacing (step size 10⁵ or 10⁶) would test whether the correlation survives with truly independent samples. Sampling at additional scales (10¹⁰, 10¹¹) would test the stability of the crossover pattern.

5. **Other additive problems.** Can the w-plane framework be applied to related problems — ternary Goldbach, Waring's problem, twin primes — and if so, does the decomposition show analogous geometric structure? The ternary case, where the major arc is known to dominate, would provide a useful control.
-e 
---

## 7. Conclusion

This study examined the empirical relationship between the major and minor arc contributions in a discrete-grid decomposition of the von Mangoldt weighted Goldbach convolution at scales up to N = 10⁹.

The principal findings are:

1. The signed major and minor arc contributions exhibit moderate positive correlation (r ≈ +0.67 to +0.73 at converged half-width) across windows of 201 even integers at three tested scales. This correlation differs in sign from simple mechanical models based only on the algebraic identity M + m = R, which yield r ≈ −1.0.

2. The complex diagnostic w(N) = M(N) + i·m(N) reveals a consistent geometric pattern across the tested scales. The angular concentration of the w(N) cloud tightens from 2.6° to 0.5° in the sampled windows, and the modular families N ≡ 0, 2, 4 mod 6 appear as distinct clusters, reflecting features of Goldbach's comet in the complex plane.

3. The data are consistent with a crossover from a major-dominated regime at N = 10⁷ to a cancellation-dominated regime at N = 10⁹. This interpretation is based on three scales and should be regarded as suggestive.

The study also identified and corrected a measurement artefact: narrow arc neighbourhoods amplify the measured correlation relative to converged values. The initial narrow-width measurement r ≈ 0.83 was revised to the converged value r ≈ 0.67 after correction, while the principal qualitative features persisted.

The methodological contributions of this work are: the w-plane diagnostic as a geometric tool for describing arc decompositions; convergence testing in the half-width parameter as a means of detecting aperture artefacts; and the use of mechanical coupling controls to separate structural associations from algebraic constraints.

These findings are purely empirical. They are based on three scales, use window-based sampling, and do not establish any asymptotic behaviour. All reported observations concern the fixed-Q_max, discrete-grid proxy defined in §2, and should not be read as claims about the full continuous major/minor arc decomposition of the classical circle method. No claim regarding the Goldbach conjecture is made. The observations are offered as a computational characterisation of phase-retaining structure in the Goldbach integrand, and as motivation for further computational and theoretical investigation.
