#!/usr/bin/env python3
"""Finite check of linear-block microcanonical nonequivalence."""

import argparse
from math import ceil, lgamma, log, sqrt

import numpy as np
from scipy.integrate import quad
from scipy.stats import norm


def negative_binomial_sum(number, maximum):
    """P(sum of number Geom_0(1/2) variables equals s), 0<=s<=maximum."""
    if number == 0:
        result = np.zeros(maximum + 1)
        result[0] = 1.0
        return result
    s = np.arange(maximum + 1, dtype=np.float64)
    logs = (
        np.array([lgamma(number + value) for value in s])
        - lgamma(number)
        - np.array([lgamma(value + 1) for value in s])
        - (number + s) * log(2.0)
    )
    return np.exp(logs)


def folded_sum_pmf(start, number, maximum):
    """Exact coefficients through maximum for consecutive folded costs."""
    explicit = []
    index = start
    while index < start + number and 2 * 3**index <= maximum:
        explicit.append(2 * 3**index)
        index += 1

    result = np.array([1.0])
    for cap in explicit:
        atom = 2.0 ** (-np.arange(1, cap + 1, dtype=np.float64))
        atom /= 1.0 - 2.0 ** (-cap)
        result = np.convolve(result, atom)[: maximum + 1]

    remaining = number - len(explicit)
    bulk = negative_binomial_sum(remaining, maximum)
    for bulk_index in range(index, start + number):
        cap = 2 * 3**bulk_index
        if cap > 1074:
            break
        bulk /= 1.0 - 2.0 ** (-cap)
    return np.convolve(result, bulk)[: maximum + 1]


def normal_limit(rho, offset):
    canonical_scale = sqrt(2.0 * rho)
    conditioned_scale = sqrt(2.0 * rho * (1.0 - rho))
    integral, _ = quad(
        lambda value: abs(
            norm.pdf(value, loc=rho * offset, scale=conditioned_scale)
            - norm.pdf(value, loc=0.0, scale=canonical_scale)
        ),
        -np.inf,
        np.inf,
        epsabs=1e-12,
        points=None,
    )
    return 0.5 * integral


def finite_tv(level, rho, offset):
    block = int(round(rho * level))
    remainder = level - block
    cost = int(round(level + offset * sqrt(level)))
    maximum = max(cost, int(ceil(block + 15 * sqrt(2 * block) + 100)))
    block_pmf = folded_sum_pmf(remainder, block, maximum)
    remainder_pmf = folded_sum_pmf(0, remainder, cost)

    weights = np.zeros(maximum + 1)
    for value in range(min(cost, maximum) + 1):
        weights[value] = block_pmf[value] * remainder_pmf[cost - value]
    conditioning_mass = weights.sum()
    conditioned = weights / conditioning_mass
    tail = max(0.0, 1.0 - block_pmf.sum())
    tv = 0.5 * (np.abs(conditioned - block_pmf).sum() + tail)
    return block, cost, tv, conditioning_mass, conditioned.sum(), block_pmf.sum()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--levels", type=int, nargs="*", default=(100, 200, 500, 1000))
    parser.add_argument("--rhos", type=float, nargs="*", default=(0.25, 0.5, 0.75))
    parser.add_argument("--offset", type=float, default=0.0)
    args = parser.parse_args()

    for rho in args.rhos:
        if not 0.0 < rho < 1.0:
            parser.error("every rho must lie strictly between zero and one")
        limit = normal_limit(rho, args.offset)
        print(f"rho={rho:.6g} gaussian_limit={limit:.12f}")
        for level in args.levels:
            block, cost, tv, mass, qsum, psum = finite_tv(
                level, rho, args.offset
            )
            print(
                f"  ell={level:4d} r={block:4d} k={cost:4d} "
                f"TV={tv:.12f} error={tv-limit:+.3e} "
                f"P(K=k)={mass:.6e} sums=({qsum:.12f},{psum:.12f})"
            )
            assert abs(qsum - 1.0) < 1e-12
            assert abs(psum - 1.0) < 1e-12


if __name__ == "__main__":
    main()
