# E-006: equivalence of ensembles at fixed 3-adic precision

This experiment projects the fixed-cost distribution
`g_ell(k,a)/sum_b g_ell(k,b)` modulo `3^r`, for fixed `r=1,2,3`, and
compares it in total variation with Tao's Syracuse law `mu_r` computed
by the independent E-002 recursion.

Run:

```sh
python3 fixed_precision_projection.py --max-ell 12 --max-precision 3
```

The tested costs are `k=ell` and `k=ell+floor(sqrt(ell))`. The theorem
in H-007 predicts convergence for every fixed precision and every
bounded central-limit offset. The finite computation checks the
direction and the implementation; it does not supply the proof.

The final two output rows with the default arguments are:

```text
ell=12 k=12 r=1:TV=0.020501489 r=2:TV=0.040237167 r=3:TV=0.06227528
ell=12 k=15 r=1:TV=0.043098385 r=2:TV=0.079022446 r=3:TV=0.12677471
```

The implementation also accepts growing diagnostic precision. For
example:

```sh
python3 fixed_precision_projection.py --max-ell 12 --max-precision 6
```

At `ell=12`, `k=ell`, it gives TV distances `0.06227528` at `r=3`
and `0.11369417` at `r=6`. These small levels do not test the theorem
for an asymptotic sequence `r=o(ell)`.

The full-precision diagnostic through level 13 is:

```sh
python3 fixed_precision_projection.py --max-ell 13 --max-precision 13
```

Its final rows give TV `0.25194861` for `k=ell` and `0.26078317` for
`k=ell+floor(sqrt(ell))` at `r=ell`. The script computes finite
distances only and makes no extrapolation from them.

At every reported level the script also checks the pointwise bound

```text
max_a p_ell,k^(r)(a)/mu_r(a) <= 1/P(K_ell=k),
```

using the exact bounded-composition count and the canonical folded-cost
normalization. It raises an exception if the domination from H-009
fails.

## Evidence (Rule 9a)

```
Command:      python3 fixed_precision_projection.py --max-ell 13 --max-precision 13
Commit:       73a0ada
Date:         2026-08-17
Environment:  Linux, Python 3.12.3, numpy 2.5.1
Exit:         0
Output:       experiments/E-006-equivalence-ensembles-fixed-precision/output/fixed_precision_projection_ell13_20260817.txt
Checked:      producer, against main.tex's Empirical Result 4.7. All four
              full-precision and precision-one values the paper quotes at
              ell=13 reproduce exactly: 0.25194861 -> 0.25195 (k=ell,
              r=13), 0.26078317 -> 0.26078 (k=ell+3, r=13),
              0.018847069 -> 0.01885 (k=ell, r=1), 0.040033188 ->
              0.04003 (k=ell+3, r=1). Added because critique round 1
              (2026-08-17, both critics independently) found these four
              values quoted in the paper with no committed output behind
              them: the block below certifies only the ell<=12, r<=3 run,
              so the ell=13 figures rested on README prose.

Command:      python3 fixed_precision_projection.py --max-ell 12 --max-precision 3
Commit:       a31aaec
Date:         2026-08-15
Environment:  Linux, Python 3.12.3, numpy 2.5.1
Exit:         0
Output:       output/fixed_precision_projection_20260815.txt
Checked:      fresh run's final two rows are
              "ell=12 k=12 r=1:TV=0.020501489 r=2:TV=0.040237167
              r=3:TV=0.06227528" and
              "ell=12 k=15 r=1:TV=0.043098385 r=2:TV=0.079022446
              r=3:TV=0.12677471", matching the values quoted above
              exactly (this run also backs the 0.04310/0.07902/0.12677
              figures cited in main.tex).
```
