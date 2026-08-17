# Reproducibility material for "Wirsching's 2003 Conjectures on Positive Predecessor Density"

Code and data for R. A. Tavares, *Wirsching's 2003 Conjectures on
Positive Predecessor Density: A Proof, a Microcanonical Decomposition,
and a Certified Numerical Test*.

The paper resolves what can be resolved of the three conjectures left
open by G. J. Wirsching, *On the problem of positive predecessor density
in 3n+1 dynamics*, Discrete and Continuous Dynamical Systems 9(3) (2003),
771-787. Conjecture 1 is proved. Conjecture 2 is analysed and stays open.
Conjecture 3 is tested numerically with a zero-truncation-error
certificate.

## Layout

Folders follow the paper's sections. Each carries its own README saying
what it verifies, how to run it, and what output to expect.

| Folder | Paper section | Verifies |
|---|---|---|
| `section-3-conjecture-1/` | §3, Theorem 3.1 | the proof of Conjecture 1: the generating-function cancellation, in exact arithmetic, plus finite checks |
| `section-4-microcanonical/` | §4, Theorems 4.1 and 4.4 to 4.9, Empirical Results 4.7 to 4.11 | the microcanonical decomposition, the three ensemble regimes, and the quantile and support diagnostics |
| `section-5-conjecture-3/` | §5, Empirical Result 5.1 and Remark 5.2 | the exact evaluation of phi to depth 500 and the extrapolation of c |

## Running

```
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
```

CPython 3.12, `mpmath` and `numpy`. Nothing else. Every script takes
`--help`.

The fastest check on the paper's central claim needs no arguments and
runs in under a second:

```
python3 section-3-conjecture-1/cancellation_check.py
```

It verifies the cancellation `P_l(z) Q_l(z) = (1-z)^{-(l+1)}` coefficient
by coefficient, the support bound `sum_j (c_j - 1) = 3^l - l - 1`, and the
convolution identity `(p_l * gbar_l)(k) = C(k+l,l) / (2*3^(l-1))`, all in
exact integer and rational arithmetic, for `l = 1..8` and `k = 0..24`.

## What reproduces, and what does not

Every figure quoted in the paper comes from a run committed here, and
each experiment's README carries the exact command, the interpreter and
package versions, the exit code, and the raw output file.

Two runs are not reproducible on ordinary hardware, and the paper's
claims do not depend on re-running them:

- `E-010`'s `central_ratio.py` at `l = 18` peaks near 119 GiB and needed
  a disk-backed memmap mode. Its committed output is
  `ratio_ell18_output.txt`. Levels through `l = 16` run in minutes on a
  laptop.
- `E-010`'s `central_zeros.py` at `l = 19, 20` has the same shape of
  cost. Its committed output is `zeros_extended_output.txt`.

Both figures those runs produce are reported in the paper as computed at
those levels, with the peak memory stated as derived from the allocation
model rather than measured.

## Notation

`phi` is Wirsching's invariant density, the unique fixed point of
`W_3 f(x) = (3/2) int_{3x-2}^{3x} f`. `phi_0` is Berg and Kruppel's
explicit asymptotic for a particular solution of the truncated equation;
its closed form and all five constants are displayed in the paper, in
Section 2, so this repository does not restate them. `g_l(k,a)` are
Wirsching's generators and `gbar_l(k)` their residue average.

## License and citation

Cite the paper. Code here is provided so the paper's claims can be
checked and reused.
