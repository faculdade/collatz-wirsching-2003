#!/usr/bin/env python3
"""Exact targeted support recursion for Wirsching generators.

The predicate hit(ell, k, a) is true exactly when g_ell(k,a)>0.  It
uses equation (2.1) of Wirsching (2003), but follows only one requested
residue instead of materializing all residues modulo 3**ell.
"""

from functools import lru_cache


@lru_cache(maxsize=None)
def hit(level, cost, residue):
    if level == 0:
        return cost == 0
    if cost < 0:
        return False

    modulus = 3**level
    old_modulus = modulus // 3
    cap = 2 * old_modulus
    residue %= modulus
    for increment in range(min(cost, cap - 1) + 1):
        numerator = pow(2, increment + 1, modulus) * residue - 1
        if numerator % 3:
            continue
        predecessor = (numerator // 3) % old_modulus
        if hit(level - 1, cost - increment, predecessor):
            return True
    return False


def least_cost(level, residue, lower, upper):
    for cost in range(lower, upper + 1):
        if hit(level, cost, residue):
            return cost
    return None


def main():
    transition = None
    for level in range(1, 61):
        modulus = 3**level
        inverse_two = (modulus + 1) // 2
        central = hit(level, level, inverse_two)
        minimum = (
            least_cost(level, inverse_two, 0, level + 7)
            if level <= 22
            else None
        )
        print(
            f"ell={level:2d} inverse_two_center={int(central)} "
            f"least_cost={minimum if minimum is not None else 'not-computed'}",
            flush=True,
        )
        if central and transition is None:
            transition = level

    assert transition == 22
    assert hit(21, 21, (3**21 + 1) // 2) is False
    assert hit(22, 22, (3**22 + 1) // 2) is True


if __name__ == "__main__":
    main()
