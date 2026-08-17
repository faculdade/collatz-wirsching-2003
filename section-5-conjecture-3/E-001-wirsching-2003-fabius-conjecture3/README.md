# E-001; Certified test of Wirsching's (2003) Conjecture 3 via exact moments of φ

Related hypothesis: [`H-001`](../../notes/H-001.md)

## What was done

A numerical test with CERTIFIED (not heuristic) error of Wirsching's
(2003) Conjecture 3; the most concrete piece of a chain of 3
conjectures reducing to positive density of 3n+1 predecessors. The
central object, φ (the invariant density of the W₃ averaging operator,
which turns out to be the base-3 analogue of Rvachev's Fabius
function/"atomic function"), has EXACT RATIONAL moments via
self-similarity (X=d(2U+X)/3); allowing φ to be evaluated at
extreme-tail points (x~ℓ·3⁻ℓ) without iterating the W₃ operator and
without loss of precision, via a reduction using iterated
antiderivatives.

## Result

Run up to ℓ=500 (CLT window, u∈[−2,2]): the decisive ratio
L_ℓ=3^(1−ℓ)φ(3x_ℓ)/φ(x_ℓ) converges to the predicted value 2/3 with a
deficit of (0.580±0.001)/ℓ; a coefficient independently reproduced by
Berg-Krüppel's (1998) own φ₀ asymptotics. ln(φ/φ₀) converges to a
finite limit L=−0.619±0.001(statistical)±0.015(functional form), i.e.
c=e^L≈0.54 (honest interval 0.53-0.55). Strong uniformity across the
CLT window (dispersion <10⁻⁴ at ℓ=500, versus an individual variation
of ~100 units in ln φ). No visible log-periodic modulation.
**Conjecture 3 numerically SUPPORTED**; but it's the most concrete tip
of a chain of 3 conjectures; Conjectures 1 and 2 (above it) remain open
and are untouched by this test. See H-001 for the full analysis.

## Files

- `experiment_conjecture3.py`; main script: exact moments, reduction
  by antiderivatives, φ₀ via Berg-Krüppel asymptotics (symbolic
  constants α,β,γ,δ,ε), validation (φ(1/2)=3/2 exact, odd central
  moments=0), sweep up to ℓ=300.
- Extension to ℓ∈{350,400,450,500} (6+ decimals, out-of-sample test
  against the tail models fit at ℓ≤300) run inline, not persisted as a
  separate script; see H-001 for the numbers.

## Reproduce

```
python3 experiment_conjecture3.py
```

Cost: moments(310)~18s, moments(510)~5min (grows ~N^4-5). Don't raise
N_MAX without need; ℓ~1000+ would require tens of minutes to hours.

## Evidence (Rule 9a)

```
Command:      python3 experiment_conjecture3.py --max-ell 500
Commit:       e007f33
Date:         2026-08-17
Environment:  Linux, Python 3.12.3, mpmath (mp.dps=100), numpy 2.5.1
Exit:         0
Output:       experiments/E-001-wirsching-2003-fabius-conjecture3/conjecture3_shifted_ell500.log
Checked:      producer, against main.tex's Empirical Result 5.1 and
              Remark 5.2. Runtime 314.8 s. This is the run that closes the
              two evidence gaps critique round 1 raised, both found by
              both critics: the script now evaluates at the paper's
              x_l^+ = x_l + 3^-(l+1) rather than the bare x_l, and reaches
              l=500 rather than stopping at 300. Every figure the paper
              quotes now reproduces from committed output: deficit
              l*(2/3 - L_l) at u=0 gives 0.8026 at l=200 and 0.8021 at
              l=500; the spread of ln(phi/phi0) across the seven u at
              l=500 is 6.333e-5; the deficit ranges from -32.124 at u=-2
              to +28.290 at u=+2; the C/sqrt(l) fit over l=200..500 gives
              L_inf=-0.618860 and coeff=-0.7916, and the l>=350 sub-range
              moves L_inf by 2.99e-6.
              The run also reports the second model form, which the paper
              does not currently state: C/ln^2(l) gives L_inf=-0.599498
              with max residual 5.51e-5 against 2.14e-5, so both forms fit
              this range comparably. The resulting MODEL spread on L_inf
              is 0.0194, putting c in [0.5386, 0.5491]. That spread, not
              the 4e-4 sub-range stability, is the dominant systematic on
              c, exactly as this README's own earlier summary and
              notes/H-001.md both recorded.

Command:      python3 experiment_conjecture3.py
Commit:       a31aaec
Date:         2026-08-15
Environment:  Linux, Python 3.12.3, numpy 2.5.1
Exit:         0
Output:       output/experiment_conjecture3_20260815.txt
Checked:      fresh run reproduces the ℓ≤300 sweep; at ℓ=300, u=0,
              L_ℓ=0.664732 against the predicted 2/3≈0.66667, giving
              300·(2/3−L_ℓ)≈0.580, matching the deficit coefficient
              quoted above. The ℓ∈{350..500} extension quoted in
              Result is not reproduced by this script (see Files:
              it was run inline in the original session, not
              persisted); that part of the claim rests on the
              original, unpreserved run, not on this evidence block.
```
