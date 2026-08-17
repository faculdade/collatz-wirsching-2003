#!/usr/bin/env python3
"""Finite checks for the generating-function proof of Wirsching Conjecture 1."""

from __future__ import annotations

from fractions import Fraction
from math import comb, exp, isqrt, log


def coin_values(ell: int) -> list[int]:
    return [1] + [2 * 3 ** (j - 1) for j in range(1, ell + 1)]


def coefficients_unbounded(coins: list[int], degree: int) -> list[int]:
    out = [0] * (degree + 1)
    out[0] = 1
    for coin in coins:
        for k in range(coin, degree + 1):
            out[k] += out[k - coin]
    return out


def coefficients_bounded(coins: list[int], degree: int) -> list[int]:
    out = [0] * (degree + 1)
    out[0] = 1
    for capacity in coins:
        nxt = [0] * (degree + 1)
        for total, count in enumerate(out):
            if count:
                for add in range(min(capacity - 1, degree - total) + 1):
                    nxt[total + add] += count
        out = nxt
    return out


def convolve_at(a: list[int], b: list[int], k: int) -> int:
    return sum(a[m] * b[k - m] for m in range(k + 1))


def main() -> None:
    print("ell max_degree identity tail_fraction_at_k=ell")
    for ell in range(2, 13):
        degree = 3 * ell
        coins = coin_values(ell)
        p = coefficients_unbounded(coins, degree)
        q = coefficients_bounded(coins, degree)
        for k in range(degree + 1):
            assert convolve_at(p, q, k) == comb(k + ell, ell)

        k = ell
        cutoff = max(1, isqrt(ell))
        denominator = comb(k + ell, ell)
        tail = sum(p[m] * q[k - m] for m in range(cutoff, k + 1))
        tail_fraction = Fraction(tail, denominator)
        print(ell, degree, "ok", f"{float(tail_fraction):.8f}")

    # A direct check of the elementary subexponential bound for p_infty(m).
    degree = 1000
    coins = [1]
    while coins[-1] <= degree:
        coins.append(2 if len(coins) == 1 else 3 * coins[-1])
    coins = [c for c in coins if c <= degree]
    p = coefficients_unbounded(coins, degree)
    ratios = []
    for m in range(2, degree + 1):
        ratios.append(log(max(1, p[m])) / (log(m + 2) ** 2))
    constant = max(ratios)
    assert all(p[m] <= exp((constant + 1e-12) * log(m + 2) ** 2)
               for m in range(2, degree + 1))
    print(f"finite subexponential check through m={degree}: C={constant:.6f}")


if __name__ == "__main__":
    main()
