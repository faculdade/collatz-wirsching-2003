# Section 5: the microcanonical decomposition and the ensemble regimes

Six experiments back this section, plus one included only as a dependency. Each folder has its own README with
the exact command, the environment, and the committed raw output.

| Folder | Backs |
|---|---|
| `E-004-wirsching-microcanonical-multiplicity/` | the multiplicity and support diagnostics, including the first fully supported cost per level |
| `E-005-microcanonical-fourier/` | the Fourier spectrum result: the spectrum and the termwise absolute bound, which is negative and so vacuous at every level tested |
| `E-006-equivalence-ensembles-fixed-precision/` | the fixed-precision finite-level result: total-variation distances at fixed 3-adic precision, `l <= 13` |
| `E-007-wirsching-central-support/` | the targeted central-support predicate |
| `E-008-linear-block-nonequivalence/` | the linear-block nonequivalence theorem's Gaussian limit, checked against the closed form at `rho = 0.25, 0.5, 0.75` |
| `E-010-wirsching-conj2-central-zeros/` | the multiplicity and support diagnostics and the quantile diagnosis, the central-cost zeros, and the fixed-set minima |
| `E-002-syracuse-collision-mass-k-ell/` | nothing in the paper directly. It is here because `E-006` imports its Syracuse level solver, and without it that script does not run |

## The one result to look at first

the microcanonical decomposition theorem gives the decomposition exactly, at every finite level. What
is open is the hypothesis `g_l(k,a) >= eta gbar_l(k)`, and
`E-004` shows why: the minimum multiplicity is zero at `k = l` at every
level computed, so no `eta > 0` works there yet.

## Cost warning

`E-010`'s heavy levels are not laptop work. `central_ratio.py` at
`l = 18` peaks near 119 GiB and needs `--scratch-dir` disk-backed mode;
`central_zeros.py` at `l = 19, 20` is far cheaper: measured peaks of 9.09 and 27.11 GiB. Committed outputs for
both are in that folder, so nothing here requires re-running them.
Levels through `l = 16` are cheap.
