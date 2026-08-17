#!/usr/bin/env python3
"""Bin-by-bin check of H-005 against the independent Syracuse recursion."""

import numpy as np

from microcanonical_multiplicity import extend


def solve_syracuse_level(previous, ell, terms=100):
    modulus = 3**ell
    orbit_size = 2 * 3 ** (ell - 1)
    forcing_law = np.zeros(modulus)
    forcing_law[1 + 3 * np.arange(len(previous))] = previous

    orbit = np.empty(orbit_size, dtype=np.int64)
    orbit[0] = 1
    for index in range(1, orbit_size):
        orbit[index] = (2 * orbit[index - 1]) % modulus
    forcing = forcing_law[orbit]

    values = np.zeros(orbit_size)
    for shift in range(1, terms + 1):
        values += 0.5**shift * np.roll(forcing, -shift)
    law = np.zeros(modulus)
    law[orbit] = values
    return law


def microcanonical_law(ell):
    max_cost = 3**ell - ell - 1
    table = np.ones((max_cost + 1, 1), dtype=np.int64)
    table[1:, 0] = 0
    for stage in range(1, ell + 1):
        table = extend(table, stage, max_cost)

    constant = 1.0
    for stage in range(ell):
        cap = 2 * 3**stage
        constant *= 0.5 / (1.0 - 2.0 ** (-cap))
    weights = 2.0 ** (-np.arange(max_cost + 1))
    return constant * np.sum(table * weights[:, None], axis=0)


def main():
    syracuse = np.zeros(3)
    syracuse[1] = 1.0 / 3.0
    syracuse[2] = 2.0 / 3.0
    for ell in range(1, 5):
        if ell > 1:
            syracuse = solve_syracuse_level(syracuse, ell)
        canonical = microcanonical_law(ell)
        discrepancy = float(np.max(np.abs(canonical - syracuse)))
        print(f"ell={ell} mass={canonical.sum():.17g} max_difference={discrepancy:.3e}")
        assert discrepancy < 1e-14


if __name__ == "__main__":
    main()
