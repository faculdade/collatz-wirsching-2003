# Section 5: the microcanonical decomposition and the ensemble regimes

Six experiments back this section, plus one included only as a dependency. Each folder has its own README with
the exact command, the environment, and the committed raw output.

| Folder | Backs |
|---|---|
| `E-004-wirsching-microcanonical-multiplicity/` | Empirical Result 4.8, the multiplicity and support diagnostics, including the first fully supported cost per level |
| `E-005-microcanonical-fourier/` | Empirical Result 4.11, the Fourier spectrum and the termwise absolute bound, which is negative and so vacuous at every level tested |
| `E-006-equivalence-ensembles-fixed-precision/` | Empirical Result 4.7, total-variation distances at fixed 3-adic precision, `l <= 13` |
| `E-007-wirsching-central-support/` | the targeted central-support predicate |
| `E-008-linear-block-nonequivalence/` | Theorem 4.6's Gaussian limit, checked against the closed form at `rho = 0.25, 0.5, 0.75` |
| `E-010-wirsching-conj2-central-zeros/` | Empirical Results 4.9 and 4.10, the quantile diagnosis, the central-cost zeros, and the fixed-set minima |
| `E-002-syracuse-collision-mass-k-ell/` | nothing in the paper directly. It is here because `E-006` imports its Syracuse level solver, and without it that script does not run |

## The one result to look at first

Theorem 4.1 gives the decomposition exactly, at every finite level. What
is open is the hypothesis `g_l(k,a) >= eta gbar_l(k)`, and
`E-004` shows why: the minimum multiplicity is zero at `k = l` at every
level computed, so no `eta > 0` works there yet.

## Cost warning

`E-010`'s heavy levels are not laptop work. `central_ratio.py` at
`l = 18` peaks near 119 GiB and needs `--scratch-dir` disk-backed mode;
`central_zeros.py` at `l = 19, 20` is comparable. Committed outputs for
both are in that folder, so nothing here requires re-running them.
Levels through `l = 16` are cheap.
