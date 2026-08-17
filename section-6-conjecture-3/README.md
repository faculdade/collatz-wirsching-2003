# Section 6: the certified numerical test of Conjecture 3

`E-001-wirsching-2003-fabius-conjecture3/` holds everything.

## What is certified, and what is not

The evaluation of `phi` is exact. Moments are computed as exact rationals
from the self-similarity `X = (2U + X)/3`, combined with an
antiderivative reduction, so the truncation error is exactly zero at
every point the paper reports, not merely small. Working precision is
100 digits.

The extrapolation of `c` is not certified, and the paper says so. Two
model forms fit the tested range comparably:

```
C/sqrt(l)  : L_inf = -0.618860, c = 0.5386, max residual 2.14e-5
C/ln^2(l)  : L_inf = -0.599498, c = 0.5491, max residual 5.51e-5
```

so the measurement pins `c` only to `[0.539, 0.549]`. The script prints
both fits and their spread; that spread, not the sub-range stability, is
the dominant systematic.

## Running

```
cd E-001-wirsching-2003-fabius-conjecture3
python3 experiment_conjecture3.py --max-ell 500      # the paper's run, ~5 min
python3 experiment_conjecture3.py --max-ell 100      # quick, under a minute
```

The default evaluation point is the paper's `x_l^+ = x_l + 3^-(l+1)`,
following Wirsching's (7.5). `--bare-point` evaluates at `x_l` instead
and is kept only for comparison: it gives a deficit coefficient near
0.580 rather than 0.802, and the two should not be mixed.

Committed output for the paper's run is `conjecture3_shifted_ell500.log`.
