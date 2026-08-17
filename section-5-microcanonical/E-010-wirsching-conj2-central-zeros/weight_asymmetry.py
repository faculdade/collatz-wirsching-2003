#!/usr/bin/env python3

> **NOTE ON PATHS.** References below to `notes/`, `literature/` or `HYPOTHESES.md` are to the framework repository, which is not public. Everything a reader needs to reproduce the paper's computations is in this repository.

"""Exact computation of Wirsching's weight sum S_ell(k,a).

Written 2026-08-16 to close a Rule 9a traceability gap: E-010's README
and notes/H-013.md both quote S_3(3,a) = 9/7 and 12/7 as the
counterexample showing the weights do not sum to 1 pointwise, and both
cited a producer script at a path that does not exist in this
repository. Documenting the values as hand-derived was not enough; the
critique that found this (H-015 round 4) was right that a quoted exact
value needs a producer, not a note saying someone once did it by hand.

The object. Dividing the g_ell recursion by gbar_ell(k) gives weights
whose sum over the admissible j's for a fixed residue a is

    S_ell(k,a) = sum over admissible j of
                 gbar_(ell-1)(k-j) / gbar_ell(k),

where gbar_ell(k) = (1/(2*3^(ell-1))) * sum over units b of g_ell(k,b)
is the Haar mean, and j is admissible for a when 2^(j+1)*a = 1 mod 3,
which fixes the parity of j from a mod 3 alone.

The point of the computation: S_ell(k,a) is NOT identically 1; it takes
exactly two values, one per residue class of a mod 3. That
non-constancy is what matters for H-013's closure, and it holds under
any normalization.

The average is 3/2, not 1, and this is a theorem, not a measurement.
main.tex's rem:no-monotone-certificate derives it: the bounded-urn
recursion gives sum_{j<c} gbar_(ell-1)(k-j) = 3*gbar_ell(k) for every
(ell,k), and every j is admissible for exactly one residue class of a
mod 3, so

    S_ell(k, a=1 mod 3) + S_ell(k, a=2 mod 3) = 3   identically,

hence the Haar average over the units is 3/2 at every level. A convex
combination averaging to 1 would need that sum to be 2. This script
checks the identity rather than asserting an average.

An earlier version of this script printed "the MEAN is 1". That came
from averaging over ALL residues mod 3^ell, where the inadmissible
class contributes an artificial zero and the units are exactly 2/3 of
the total, so the all-residue average is always 2/3 of the unit
average. Obtaining 1 that way is a consequence of the measure chosen,
not an identity. Found by H-015 round 5; see notes/H-022.md, closed.

On the factor of 3 in notes/H-013.md's displayed weight. It is
correct under the reading that bar_g there means the raw
bounded-composition count B, since B_ell = 2*3^(ell-1)*gbar_ell makes
3*B_(ell-1)(k-j)/B_ell(k) identically equal to
gbar_(ell-1)(k-j)/gbar_ell(k), which is the weight this script uses.
The notation is overloaded, not wrong. See notes/H-022.md, opened and
closed 2026-08-16.

Everything is done in exact rationals; no floating point anywhere.
"""
from __future__ import annotations

import argparse
from fractions import Fraction


def bounded_composition_counts(ell: int, k_max: int) -> list[int]:
    """B_ell(k): number of (j_1,...,j_ell) with 0 <= j_i < 2*3^(i-1)
    summing to k. This is 2*3^(ell-1) * gbar_ell(k), Wirsching's own
    identity, and is the same quantity central_ratio.py checks its rows
    against."""
    counts = [0] * (k_max + 1)
    counts[0] = 1
    for index in range(1, ell + 1):
        cap = 2 * 3 ** (index - 1)
        prefix = [0] * (k_max + 2)
        for cost in range(k_max + 1):
            prefix[cost + 1] = prefix[cost] + counts[cost]
        counts = [prefix[cost + 1] - prefix[max(0, cost - cap + 1)]
                  for cost in range(k_max + 1)]
    return counts


def gbar(ell: int, k: int, k_max: int) -> Fraction:
    """B_ell(k) / (2*3^(ell-1)), extended to ell = 0 by the same formula.

    At ell = 0 that denominator is 2*3^(-1) = 2/3, so gbar_0(0) = 3/2,
    NOT 1. An earlier version returned 1 here, which broke the identity
    S(a=1) + S(a=2) = 3 at ell = 1 (it gave 2). Harmless for the levels
    this project quotes, ell >= 3, but wrong, and it is exactly the kind
    of edge case a hardcoded conclusion hides; found by H-015 round 5.
    """
    if ell == 0:
        return Fraction(3, 2) if k == 0 else Fraction(0)
    return Fraction(bounded_composition_counts(ell, k_max)[k],
                    2 * 3 ** (ell - 1))


def admissible_js(ell: int, a: int) -> list[int]:
    """j in [0, 2*3^(ell-1)) with 2^(j+1) * a = 1 (mod 3).

    Since 2 = -1 mod 3, this is (-1)^(j+1) * a = 1 mod 3, so a = 1 mod 3
    forces j odd and a = 2 mod 3 forces j even. Nothing else about a
    matters, which is why S_ell(k,a) can take only two values per level.
    """
    cap = 2 * 3 ** (ell - 1)
    return [j for j in range(cap) if (pow(2, j + 1, 3) * a) % 3 == 1]


def weight_sum(ell: int, k: int, a: int) -> Fraction:
    k_max = k
    denominator = gbar(ell, k, k_max)
    if denominator == 0:
        return Fraction(0)
    total = Fraction(0)
    for j in admissible_js(ell, a):
        if 0 <= k - j <= k_max:
            total += gbar(ell - 1, k - j, k_max)
    return total / denominator


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ell", type=int, default=3)
    parser.add_argument("--k", type=int, default=3)
    args = parser.parse_args()

    ell, k = args.ell, args.k
    modulus = 3**ell
    print(f"S_{ell}({k},a) for every residue a mod {modulus}, exact rationals")
    print(f"gbar_{ell}({k}) = {gbar(ell, k, k)}")
    print()

    by_class: dict[int, set[Fraction]] = {0: set(), 1: set(), 2: set()}
    for a in range(modulus):
        value = weight_sum(ell, k, a)
        by_class[a % 3].add(value)

    for residue_class in (0, 1, 2):
        values = by_class[residue_class]
        assert len(values) == 1, (
            f"a = {residue_class} mod 3 gave more than one value: {values}"
        )
        print(f"  a = {residue_class} mod 3:  S = {values.pop()}")

    print()
    one = weight_sum(ell, k, 1)
    two = weight_sum(ell, k, 2)
    total = one + two
    units = [a for a in range(modulus) if a % 3]
    unit_mean = sum(weight_sum(ell, k, a) for a in units) / len(units)
    print(f"S(a=1 mod 3) + S(a=2 mod 3) = {total}")
    assert total == 3, (
        f"the identity S(a=1)+S(a=2)=3 failed at ell={ell}, k={k}: "
        f"got {total}. See main.tex rem:no-monotone-certificate."
    )
    print(f"mean over the {len(units)} units    = {unit_mean}")
    assert unit_mean == Fraction(3, 2), (
        f"unit mean should be 3/2, got {unit_mean}"
    )
    print()
    print("Both check out. The identity is proved in main.tex's")
    print("rem:no-monotone-certificate; this script verifies it exactly")
    print("at the given level rather than asserting it.")


if __name__ == "__main__":
    main()
