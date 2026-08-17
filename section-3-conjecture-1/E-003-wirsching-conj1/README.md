# E-003: finite checks for Wirsching Conjecture 1

The script independently constructs the coefficients

```text
P_ell(z) = product_j (1-z^c_j)^(-1),
Q_ell(z) = product_j (1-z^c_j)/(1-z),
```

and verifies coefficient by coefficient that
`P_ell(z)Q_ell(z)=(1-z)^(-ell-1)`. It also reports finite convolution
tail fractions and checks the elementary `exp(C log^2 m)` upper bound
for the infinite coin system through `m=1000`.

Run `python3 check_generating_identity.py`.

The computation is a diagnostic. The proof is the generating-function
cancellation and tail estimate recorded in H-002.

## Evidence (Rule 9a)

```
Command:      python3 check_generating_identity.py
Commit:       a31aaec
Date:         2026-08-15
Environment:  Linux, Python 3.12.3, numpy 2.5.1
Exit:         0
Output:       output/check_generating_identity_20260815.txt
Checked:      fresh run confirms "ok" (P_ell(z)Q_ell(z) identity holds
              coefficient by coefficient) at every level ell=2..12,
              and the finite subexponential check through m=1000
              gives C=0.372204.
```
