#!/usr/bin/env python3
"""Fixed-precision projections of Wirsching's microcanonical law."""

import argparse
import importlib.util
from math import isqrt, prod
from pathlib import Path

import numpy as np


EXPERIMENTS = Path(__file__).resolve().parents[1]


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, EXPERIMENTS / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def project(row, ell, precision):
    modulus = 3**precision
    projected = np.zeros(modulus)
    residues = np.arange(3**ell) % modulus
    np.add.at(projected, residues, row.astype(np.float64))
    return projected / projected.sum()


def canonical_cost_probability(total, ell, cost):
    normalization = prod(
        1.0 / (2.0 * (1.0 - 2.0 ** (-(2 * 3**stage))))
        for stage in range(ell)
    )
    return total * 2.0 ** (-cost) * normalization


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-ell", type=int, default=12)
    parser.add_argument("--max-precision", type=int, default=3)
    args = parser.parse_args()
    micro = load(
        "micro",
        "E-004-wirsching-microcanonical-multiplicity/microcanonical_multiplicity.py",
    )
    syracuse = load(
        "syracuse",
        "E-002-syracuse-collision-mass-k-ell/experiment_k_ell.py",
    )

    references = {}
    law = np.zeros(3)
    law[1] = 1.0 / 3.0
    law[2] = 2.0 / 3.0
    references[1] = law.copy()
    for precision in range(2, args.max_precision + 1):
        law = syracuse.solve_level(law, precision)
        references[precision] = law.copy()

    k_max = args.max_ell + isqrt(args.max_ell) + 2
    table = np.ones((k_max + 1, 1), dtype=np.int64)
    table[1:, 0] = 0
    for ell in range(1, args.max_ell + 1):
        table = micro.extend(table, ell, k_max)
        if ell < args.max_precision:
            continue
        for cost in (ell, ell + isqrt(ell)):
            values = []
            row = table[cost]
            cost_probability = canonical_cost_probability(
                int(row.sum()), ell, cost
            )
            for precision in range(1, args.max_precision + 1):
                marginal = project(row, ell, precision)
                tv = 0.5 * float(np.abs(marginal - references[precision]).sum())
                positive = references[precision] > 0
                likelihood_ratio = float(
                    np.max(marginal[positive] / references[precision][positive])
                )
                if likelihood_ratio > (1.0 + 1e-10) / cost_probability:
                    raise RuntimeError(
                        "microcanonical-to-canonical domination failed "
                        f"at ell={ell}, k={cost}, r={precision}"
                    )
                values.append(f"r={precision}:TV={tv:.8g}")
            print(f"ell={ell:2d} k={cost:2d} " + " ".join(values))


if __name__ == "__main__":
    main()
