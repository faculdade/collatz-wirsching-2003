# E-001; Certified test of Wirsching's (2003) Conjecture 3 via exact moments of φ

> **NOTE ON RECORDS.** Hypothesis identifiers such as H-013, the hypothesis tracker and the literature index refer to the framework repository, which is not public. Everything a reader needs to reproduce the paper's computations is in this repository.


Related hypothesis: ``H-001``

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

Updated 2026-08-17 after critique round 2 found this section contradicting
the manuscript on four points; the four corrections are named below rather
than silently applied, since this file is public and someone may have read
the old version.

Run to `ℓ=500` (CLT window, `u∈[−2,2]`), evaluated at Wirsching's shifted
point `x_ℓ⁺ = x_ℓ + 3^(−ℓ−1)`:

- The decisive ratio `L_ℓ = 3^(1−ℓ)·φ(3x_ℓ⁺)/φ(x_ℓ⁺)` approaches `2/3`
  at every `u` tested, monotonically in `ℓ`, with a central-sequence
  deficit `ℓ·(2/3 − L_ℓ)` of `0.8026` at `ℓ=200` and `0.8021` at
  `ℓ=500`.
- `ln(φ/φ₀)` has a limit consistent with a constant, near-uniformly in
  `u`: the spread across the seven tested `u` at `ℓ=500` is
  `6.33×10⁻⁵`.
- Extrapolating that limit is model-dependent, and this is the dominant
  uncertainty. `C/√ℓ` gives `L∞ = −0.618860`, so `c = 0.5386`, with
  maximum residual `2.1×10⁻⁵`; `C/ln²ℓ` gives `L∞ = −0.599498`, so
  `c = 0.5491`, with maximum residual `5.5×10⁻⁵`. The two fit this range
  comparably, so `c` is pinned only to `[0.539, 0.549]`.

**What the four corrections were.**

1. The deficit coefficient was quoted as `0.580/ℓ`. That is the value at
   the **bare** point `x_ℓ`, which the script used until 2026-08-17;
   Wirsching's (7.5) and the paper both evaluate at `x_ℓ⁺`, where the
   coefficient is `0.802`. Both are computable now (`--bare-point`
   recovers the old one) and they must not be mixed.
2. "No visible log-periodic modulation" is the wrong statement. `L_ℓ` is
   not merely free of visible modulation, it is **blind** to a
   `1`-periodic factor `Q(log₃x)` by construction: the same-phase pair
   satisfies `x_{3ℓ}/x_ℓ = 3^(1−2ℓ)`, an integer shift in `log₃`, so `Q`
   cancels identically at any amplitude.
3. "Conjectures 1 and 2 remain open" is out of date for Conjecture 1,
   which the paper proves. Conjecture 2 is open.
4. The script stopped at `ℓ=300` and the `ℓ=350..500` extension was run
   inline and not persisted. It now runs to `ℓ=500` and the output is
   committed as `conjecture3_shifted_ell500.log`.

Conjecture 3 remains numerically supported and unproved by this test. See
`H-001` in the framework repository for the full analysis.

## Files

- `experiment_conjecture3.py`; main script: exact moments, reduction
  by antiderivatives, φ₀ via Berg-Krüppel asymptotics (symbolic
  constants α,β,γ,δ,ε), validation (φ(1/2)=3/2 exact, odd central
  moments=0), sweep to ℓ=500 behind `--max-ell`, both extrapolation
  models, and the paper's shifted evaluation point by default.
- `conjecture3_shifted_ell500.log`; the committed output of the ℓ=500
  run, 314.8 s, from which every figure above is read.

## Reproduce

```
python3 experiment_conjecture3.py
```

Cost: moments(310)~18s, moments(510)~5min (grows ~N^4-5). Don't raise
N_MAX without need; ℓ~1000+ would require tens of minutes to hours.

## Evidence (Rule 9a)

```
Command:      python3 same_phase_drift.py
Commit:       8911937
Date:         2026-08-17
Environment:  Linux, Python 3.12.3, mpmath (mp.dps=100)
Exit:         0
Output:       section-6-conjecture-3/E-001-wirsching-2003-fabius-conjecture3/same_phase_drift_output.txt
Checked:      producer, against main.tex's the same-phase remark. All twelve figures
              reproduce: observed drifts -0.0128, +0.0387, +0.0320,
              +0.0271; predicted +0.1058, +0.0473, +0.0335, +0.0273;
              mismatches 112%, 18%, 4%, 1%. Written because critique round
              3 found all twelve quoted with no committed output behind
              any of them, and the (150,450) pair unrecoverable even from
              the older bare-point log, which stopped at l=300.

Command:      python3 experiment_conjecture3.py --max-ell 500
Commit:       e007f33
Date:         2026-08-17
Environment:  Linux, Python 3.12.3, mpmath (mp.dps=100), numpy 2.5.1
Exit:         0
Output:       section-6-conjecture-3/E-001-wirsching-2003-fabius-conjecture3/conjecture3_shifted_ell500.log
Checked:      producer, against main.tex's the numerical test of (*4) and (*5) and
              the same-phase remark. Runtime 314.8 s. This is the run that closes the
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
              The run also reports the second model form, which the paper states in the numerical test of (*4) and (*5): C/ln^2(l) gives L_inf=-0.599498
              with max residual 5.51e-5 against 2.14e-5, so both forms fit
              this range comparably. The resulting MODEL spread on L_inf
              is 0.0194, putting c in [0.5386, 0.5491]. That spread, not
              the 4e-4 shift from adding a 1/l term, is the dominant systematic on
              c, exactly as this README's own earlier summary and
              H-001 both recorded.

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
