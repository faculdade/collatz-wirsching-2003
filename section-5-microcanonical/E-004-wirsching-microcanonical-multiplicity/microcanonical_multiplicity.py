#!/usr/bin/env python3
"""Exact finite-level multiplicities g_ell(k,a) in Wirsching (2003).

The recursion is equation (2.1):

  g_(ell+1)(k,a) = sum_{0 <= j < 2*3^ell}
                   g_ell(k-j, (2^(j+1)*a-1)/3).

We run it forward.  If b is a residue modulo 3^(ell-1), then the unique
residue a modulo 3^ell contributing at increment j is

  a = (3*b+1) * 2^(-(j+1)) mod 3^ell.

The script reports the minimum positive-residue multiplicity divided by
the Haar average, together with the fraction of residues hit.  It is a
finite diagnostic and does not establish Wirsching's Conjecture 2.
"""

import argparse
from math import isqrt

import numpy as np


def extend(previous, ell, k_max):
    previous_modulus = previous.shape[1]
    modulus = 3**ell
    current = np.zeros((k_max + 1, modulus), dtype=np.int64)
    residues = np.arange(previous_modulus, dtype=np.int64)
    cap = 2 * 3 ** (ell - 1)
    inv2 = pow(2, -1, modulus)
    inverse_power = 1

    for increment in range(min(k_max, cap - 1) + 1):
        inverse_power = (inverse_power * inv2) % modulus
        targets = ((3 * residues + 1) * inverse_power) % modulus
        for old_cost in range(k_max - increment + 1):
            values = previous[old_cost]
            if np.any(values):
                current[old_cost + increment, targets] += values
    return current


def bounded_composition_counts(ell, k_max):
    counts = [0] * (k_max + 1)
    counts[0] = 1
    for stage in range(ell):
        cap = 2 * 3**stage
        updated = [0] * (k_max + 1)
        for total, count in enumerate(counts):
            for increment in range(min(cap - 1, k_max - total) + 1):
                updated[total + increment] += count
        counts = updated
    return counts


def diagnostics(table, ell, offsets):
    units = 2 * 3 ** (ell - 1)
    unit_mask = np.arange(3**ell) % 3 != 0
    for offset in offsets:
        cost = ell + offset * isqrt(ell)
        if not 0 <= cost < table.shape[0]:
            continue
        row = table[cost, unit_mask]
        total = int(row.sum())
        hit = int(np.count_nonzero(row))
        minimum = int(row.min())
        maximum = int(row.max())
        min_ratio = minimum * units / total if total else float("nan")
        max_ratio = maximum * units / total if total else float("nan")
        collision = units * float(np.sum(row.astype(np.float64) ** 2)) / total**2 if total else float("nan")
        print(
            f"ell={ell:2d} offset={offset:+d} k={cost:2d} total={total} "
            f"hit={hit}/{units} hit_fraction={hit/units:.8f} "
            f"min_over_mean={min_ratio:.8g} max_over_mean={max_ratio:.8g} "
            f"normalized_collision={collision:.8g}"
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-ell", type=int, default=12)
    parser.add_argument("--offsets", type=int, nargs="*", default=(-1, 0, 1))
    args = parser.parse_args()
    if args.max_ell < 1:
        parser.error("--max-ell must be positive")

    k_max = args.max_ell + max(abs(value) for value in args.offsets) * isqrt(args.max_ell) + 2
    table = np.ones((k_max + 1, 1), dtype=np.int64)
    table[1:, 0] = 0

    for ell in range(1, args.max_ell + 1):
        table = extend(table, ell, k_max)
        expected = bounded_composition_counts(ell, k_max)
        totals = table.sum(axis=1)
        if not np.array_equal(totals, np.array(expected, dtype=np.int64)):
            raise RuntimeError(f"bounded-composition check failed at ell={ell}")
        diagnostics(table, ell, args.offsets)


if __name__ == "__main__":
    main()
