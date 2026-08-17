# E-005: Fourier spectrum of microcanonical multiplicities

For the fixed-cost distribution

    p_(ell,k)(a) = g_ell(k,a) / sum_b g_ell(k,b),

this experiment compares `p_(ell,k)` with uniform measure on the unit
residues. Additive Fourier inversion gives an exact sufficient lower
bound for the minimum probability in terms of the l1 norm of the
difference of their Fourier transforms.

Run:

```sh
python3 microcanonical_fourier.py --max-ell 12
```

The script also reports the normalized collision excess, the largest
individual coefficient, and the Fourier l1 mass grouped by conductor.
A negative lower bound only says that the triangle-inequality argument
is insufficient. It is not evidence for a negative probability or an
asymptotic obstruction.

## Evidence (Rule 9a)

```
Command:      python3 microcanonical_fourier.py --max-ell 12
Commit:       a31aaec
Date:         2026-08-15
Environment:  Linux, Python 3.12.3, numpy 2.5.1
Exit:         0
Output:       output/microcanonical_fourier_20260815.txt
Checked:      fresh run reproduces the ell=4..12 sweep at both tested
              costs; fourier_lower_bound is negative at every level
              shown (e.g. -194.53 at ell=12,k=12), consistent with the
              "not evidence for an obstruction" reading above.
```
