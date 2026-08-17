# E-007: targeted central support for Wirsching generators

This experiment evaluates the exact Boolean predicate
`g_ell(k,a)>0` without constructing the full table modulo `3^ell`.
The recursion is equation (2.1) of Wirsching (2003), truncated only by
the requested total cost, so its output is exact.

The motivating candidate obstruction was

```text
g_ell(ell, 2^(-1) mod 3^ell) = 0 for every ell.
```

It agrees with the exhaustive E-004 tables through `ell=16` and remains
true through `ell=21`. It fails at `ell=22`. The same residue family is
centrally reachable at every tested level from 22 through 60. Thus the
finite support holes do not provide this proposed counterexample to
condition `(?3)`.

Run:

```sh
python3 central_support.py
```

The assertions certify the transition at levels 21 and 22. This is a
diagnostic and neither proves eventual full central support nor proves
Wirsching's Conjecture 2.

The targeted predicate was independently compared with the complete
Boolean tables from E-004 for every unit residue, every cost through 12,
and every level through 8. All entries agreed.

## Evidence (Rule 9a)

```
Command:      python3 central_support.py
Commit:       a31aaec
Date:         2026-08-15
Environment:  Linux, Python 3.12.3, numpy 2.5.1
Exit:         0
Output:       output/central_support_20260815.txt
Checked:      fresh run gives inverse_two_center=0 for ell=1..21 and
              =1 starting at ell=22, matching the "fails at ell=22"
              transition claimed above exactly.
```
