# Section 4: the proof of Conjecture 1

Verifies Theorem 4.1, which proves Wirsching's Conjecture 1, the
implication `(*2) => (*1)`, and with it Corollary 4.2: condition `(*3)`
implies uniform positive predecessor density, unconditionally.

## `cancellation_check.py`

The proof's own step. Wirsching describes the averaged generators by an
urn count without a generating function; the proof turns that into two
generating functions over the coin set `c_0 = 1`, `c_j = 2*3^(j-1)`, and
uses their product collapsing:

```
sum_k q_l(k) z^k = prod_j (1 - z^{c_j}) / (1 - z)
sum_k p_l(k) z^k = prod_j (1 - z^{c_j})^{-1}
P_l(z) Q_l(z)    = (1 - z)^{-(l+1)}
```

Run it:

```
python3 cancellation_check.py                       # l=1..8, k=0..24
python3 cancellation_check.py --max-ell 10 --max-k 40
```

Expected output ends with `All checks passed`. Four things are checked,
in exact integer and rational arithmetic with no floating point:

1. the bounded-urn expansion against direct enumeration of occupancies,
2. the cancellation, coefficient by coefficient, against the binomial
   coefficients of `(1-z)^{-(l+1)}`,
3. the support bound `sum_j (c_j - 1) = 3^l - l - 1`, which is what
   forces the occupancy reading `0..c_j-1` of Wirsching's "capacity",
4. the convolution identity `(p_l * gbar_l)(k) = C(k+l,l)/(2*3^(l-1))`.

Runtime is under a second at the defaults and a few seconds at
`--max-ell 10 --max-k 40`.

## `partition_bound_check.py`

The proof's other arithmetic step. It needs the partition count for the
coins `1, 2, 6, 18, ...` to be `exp(O(log^2 m))`, and the paper proves
that in two clauses. This script computes the count exactly by dynamic
programming, reports the smallest constant for which the bound holds
(`0.372` up to `m = 4000`), and checks the paper's own elementary bound
against the exact count.

```
python3 partition_bound_check.py --max-m 4000
```

## `E-003-wirsching-conj1/`

The finite checks reported alongside the theorem. See its own README for
the command and committed output.
