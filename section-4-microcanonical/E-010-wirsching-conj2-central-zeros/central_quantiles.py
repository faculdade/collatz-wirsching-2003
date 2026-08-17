#!/usr/bin/env python3
"""Quantile version of the quantitative half of Wirsching's condition (?3).

H-013 asks whether the decay of

    R_ell(k) = min_a g_ell(k,a) / gbar_ell(k)

(computed by `central_ratio.py`) is extreme-value statistics of a growing
index set (2*3^(ell-1) units at level ell) or a genuine deterioration of
the left tail of the distribution of g_ell(k,a)/gbar_ell(k) over units a.
A fixed quantile of that distribution is not subject to the extreme-value
confound: the number of units contributing to, say, the 10^-3 quantile
grows with ell (there are more units below that quantile in absolute
count), but the quantile itself asks "what value does a fixed *fraction*
of the population fall below," which does not mechanically shrink just
because the population grows, the way a minimum does.

If a fixed low quantile stays flat in ell while the minimum keeps
decaying, the minimum's decay is extreme-value statistics and the tail
itself is not deteriorating. If the quantile also decays at comparable
rate to the minimum, the tail is genuinely thinning and H-013's question
(is inf_a liminf_ell R_ell(ell+d,a) positive) leans toward "no."

Uses the same exact forward recursion as `central_ratio.py` (recurrence
audited there against Wirsching's row-sum identity); this script adds
no new arithmetic, only a different statistic over the same table.
"""

from __future__ import annotations

import argparse
import time
from math import isqrt

import numpy as np


def bounded_composition_counts(ell: int, k_max: int) -> np.ndarray:
    """Number of (j_1,...,j_ell) with 0 <= j_i < 2*3^(i-1) and sum = k."""
    counts = np.zeros(k_max + 1, dtype=object)
    counts[0] = 1
    for index in range(1, ell + 1):
        cap = 2 * 3 ** (index - 1)
        prefix = np.zeros(k_max + 2, dtype=object)
        running = 0
        for cost in range(k_max + 1):
            running += counts[cost]
            prefix[cost + 1] = running
        updated = np.zeros(k_max + 1, dtype=object)
        for cost in range(k_max + 1):
            low = max(0, cost - cap + 1)
            updated[cost] = prefix[cost + 1] - prefix[low]
        counts = updated
    return counts


def extend(previous: np.ndarray, ell: int, k_max: int) -> np.ndarray:
    modulus = 3**ell
    old_modulus = modulus // 3
    current = np.zeros((k_max + 1, modulus), dtype=np.int64)
    residues = np.arange(old_modulus, dtype=np.int64)
    cap = 2 * old_modulus
    inv2 = pow(2, -1, modulus)
    inverse_power = 1
    for increment in range(min(k_max, cap - 1) + 1):
        inverse_power = (inverse_power * inv2) % modulus
        targets = ((3 * residues + 1) % modulus) * inverse_power % modulus
        for old_cost in range(k_max - increment + 1):
            current[old_cost + increment, targets] += previous[old_cost]
    return current


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-ell", type=int, default=16)
    parser.add_argument("--offsets", type=int, nargs="+", default=[5, 12],
                        help="reported costs are ell + offset; both "
                             "defaults have complete support from ell=4 "
                             "(d=+12) or ell>=10 (d=+5) per central_zeros.py")
    parser.add_argument("--quantiles", type=float, nargs="+",
                        default=[1e-4, 1e-3, 1e-2, 1e-1, 0.5],
                        help="fraction of the unit population; the "
                             "reported value is the largest v such that "
                             "at most this fraction of units has "
                             "g_ell(k,a)/gbar_ell(k) <= v")
    parser.add_argument("--thresholds", type=float, nargs="*",
                        default=[0.05, 0.1, 0.2, 0.3],
                        help="report P[R <= v] for these fixed ratio "
                             "values v, separating bulk broadening from "
                             "tail pinching to zero")
    parser.add_argument("--bucket-quantile", type=float, default=1e-3,
                        help="bottom fraction of the unit population "
                             "used for the small-representative "
                             "composition check")
    parser.add_argument("--small-bound", type=int, default=3**8,
                        help="an integer a below this bound is 'small'; "
                             "reports what fraction of the bottom "
                             "bucket has an integer representative "
                             "below it, connecting the population "
                             "statistic to H-013's fixed-a question")
    args = parser.parse_args()

    offsets = sorted(args.offsets)
    quantiles = sorted(args.quantiles)
    k_max = args.max_ell + offsets[-1]
    table = np.zeros((k_max + 1, 1), dtype=np.int64)
    table[0, 0] = 1

    print("Quantiles of g_ell(k,a)/gbar_ell(k) over unit residues a mod "
          "3^ell, k = ell + d.")
    print("Row sum checked against Wirsching's bounded-composition "
          "identity at every (ell, d). One forward pass in ell; all "
          "offsets are read from the same table per level.")
    header = ("ell   d " + "   n_units "
              + " ".join(f"  q={q:<9.1e}" for q in quantiles)
              + "      min")
    print(header)
    print("-" * len(header))

    for ell in range(1, args.max_ell + 1):
        started = time.time()
        table = extend(table, ell, k_max)
        modulus = 3**ell
        unit = np.tile(np.array([False, True, True]), modulus // 3)
        n_units = int(unit.sum())
        unit_values = np.nonzero(unit)[0]

        for offset in offsets:
            cost = ell + offset
            if cost > k_max:
                continue
            reference = bounded_composition_counts(ell, cost)
            row = table[cost][unit]
            total = int(row.sum())
            if total != int(reference[cost]):
                raise AssertionError(
                    f"row total mismatch ell={ell} k={cost}: "
                    f"{total} != {int(reference[cost])}"
                )
            if total == 0:
                # cost not yet reachable at this level (below the minimum
                # possible cost of ell); nothing to report.
                continue
            mean = total / n_units
            order = np.argsort(row, kind="stable")
            sorted_row = row[order]
            fields = []
            for q in quantiles:
                idx = min(int(q * n_units), n_units - 1)
                fields.append(f"{sorted_row[idx] / mean:11.5f}")
            smallest = int(sorted_row[0])
            argmin_a = int(unit_values[order[0]])
            print(f"{ell:3d} {offset:+3d} {n_units:10d} "
                  + " ".join(fields) + f"  {smallest / mean:7.5f}"
                  + f"  argmin_a={argmin_a}",
                  flush=True)

            thr_fields = []
            for v in args.thresholds:
                frac = float(np.count_nonzero(row <= v * mean)) / n_units
                thr_fields.append(f"P[R<={v:.2f}]={frac:8.6f}")
            print("      thresholds: " + "  ".join(thr_fields), flush=True)

            bucket_n = max(1, int(args.bucket_quantile * n_units))
            bucket_idx = order[:bucket_n]
            bucket_a = unit_values[bucket_idx]
            small_frac = float(np.count_nonzero(
                bucket_a < args.small_bound)) / bucket_n
            print(f"      bottom {args.bucket_quantile:.1e} bucket "
                  f"({bucket_n} residues): fraction with a < "
                  f"{args.small_bound} is {small_frac:.6f}", flush=True)

            del sorted_row, order
        print(f"    [ell={ell} total {time.time() - started:6.1f}s]",
              flush=True)


if __name__ == "__main__":
    main()
