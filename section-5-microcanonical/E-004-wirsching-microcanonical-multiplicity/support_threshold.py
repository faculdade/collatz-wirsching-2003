#!/usr/bin/env python3
"""Boolean support thresholds for fixed-cost Wirsching generators."""

import argparse
from math import isqrt

import numpy as np


def extend_support(previous, ell, k_max):
    old_modulus = previous.shape[1]
    modulus = 3**ell
    current = np.zeros((k_max + 1, modulus), dtype=np.bool_)
    residues = np.arange(old_modulus, dtype=np.int64)
    cap = 2 * 3 ** (ell - 1)
    inv2 = pow(2, -1, modulus)
    inverse_power = 1
    for increment in range(min(k_max, cap - 1) + 1):
        inverse_power = (inverse_power * inv2) % modulus
        targets = ((3 * residues + 1) * inverse_power) % modulus
        for old_cost in range(k_max - increment + 1):
            current[old_cost + increment, targets] |= previous[old_cost]
    return current


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-ell", type=int, default=16)
    parser.add_argument("--k-max", type=int, default=30)
    args = parser.parse_args()
    support = np.zeros((args.k_max + 1, 1), dtype=np.bool_)
    support[0, 0] = True

    for ell in range(1, args.max_ell + 1):
        support = extend_support(support, ell, args.k_max)
        units = np.arange(3**ell) % 3 != 0
        rows = support[:, units]
        cumulative = np.logical_or.accumulate(rows, axis=0)
        cumulative_levels = np.flatnonzero(np.all(cumulative, axis=1))
        exact_levels = np.flatnonzero(np.all(rows, axis=1))
        center = ell
        upper = ell + isqrt(ell)
        print(
            f"ell={ell:2d} cumulative_first="
            f"{int(cumulative_levels[0]) if len(cumulative_levels) else 'none'} "
            f"exact_first={int(exact_levels[0]) if len(exact_levels) else 'none'} "
            f"center_fraction={rows[center].mean():.8f} "
            f"upper_fraction={rows[upper].mean():.8f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
