#!/usr/bin/env python3
"""Fourier diagnostics for fixed-cost Wirsching multiplicities.

For p_(ell,k)(a)=g_ell(k,a)/sum_b g_ell(k,b), compare p with uniform
measure on the units modulo 3^ell.  Additive Fourier inversion gives

  |p(a)-u(a)| <= 3^(-ell) sum_xi |p_hat(xi)-u_hat(xi)|.

The script computes this exact sufficient bound and groups its Fourier
mass by conductor.  Finite values do not imply an asymptotic estimate.
"""

import argparse
import importlib.util
from math import isqrt
from pathlib import Path

import numpy as np


def load_microcanonical_module():
    path = Path(__file__).resolve().parents[1] / (
        "E-004-wirsching-microcanonical-multiplicity/microcanonical_multiplicity.py"
    )
    spec = importlib.util.spec_from_file_location("microcanonical_multiplicity", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def conductor(xi, ell):
    if xi == 0:
        return 0
    valuation = 0
    while xi % 3 == 0:
        xi //= 3
        valuation += 1
    return ell - valuation


def analyze(row, ell, cost):
    modulus = 3**ell
    unit_mask = np.arange(modulus) % 3 != 0
    total = float(row.sum())
    probability = row.astype(np.float64) / total
    uniform = unit_mask.astype(np.float64) / unit_mask.sum()
    difference_transform = np.fft.fft(probability - uniform)
    absolute = np.abs(difference_transform)
    spectral_l1 = float(absolute.sum())
    spectral_linf = float(absolute[1:].max())
    relative_lower_bound = 1.0 - (2.0 / 3.0) * spectral_l1
    actual_min_ratio = float(probability[unit_mask].min() * unit_mask.sum())
    collision_excess = float(unit_mask.sum() * np.sum(probability**2) - 1.0)

    by_conductor = {}
    for xi in range(1, modulus):
        level = conductor(xi, ell)
        by_conductor[level] = by_conductor.get(level, 0.0) + float(absolute[xi])
    largest = sorted(by_conductor.items(), key=lambda item: item[1], reverse=True)[:4]
    groups = ",".join(f"r{level}:{mass:.5g}" for level, mass in largest)
    primitive = [xi for xi in range(1, modulus) if xi % 3 != 0]
    primitive_linf = float(absolute[primitive].max())
    print(
        f"ell={ell:2d} k={cost:2d} min_ratio={actual_min_ratio:.8g} "
        f"collision_excess={collision_excess:.8g} spectral_l1={spectral_l1:.8g} "
        f"spectral_linf={spectral_linf:.8g} primitive_linf={primitive_linf:.8g} "
        f"fourier_lower_bound={relative_lower_bound:.8g} largest_groups={groups}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-ell", type=int, default=12)
    args = parser.parse_args()
    micro = load_microcanonical_module()
    k_max = args.max_ell + isqrt(args.max_ell) + 2
    table = np.ones((k_max + 1, 1), dtype=np.int64)
    table[1:, 0] = 0
    for ell in range(1, args.max_ell + 1):
        table = micro.extend(table, ell, k_max)
        if ell >= 4:
            analyze(table[ell], ell, ell)
            analyze(table[ell + isqrt(ell)], ell, ell + isqrt(ell))


if __name__ == "__main__":
    main()
