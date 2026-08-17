#!/usr/bin/env python3
"""Quantitative part of Wirsching's condition (?3).

Condition (?3) has two independent halves.  The first is support: every
unit residue must satisfy g_ell(k,a) > 0.  The second is quantitative:
the ratio

    R_ell(k) = min_a g_ell(k,a) / gbar_ell(k),
    gbar_ell(k) = (sum_a g_ell(k,a)) / (2*3^(ell-1)),

must stay above a fixed mu > 0.  `central_zeros.py` handles the first
half.  This script handles the second, using exact integer counts from
the same forward recursion.

Costs are reported at k = ell (the centre of Wirsching's window, where
the support is still incomplete on every level reached here) and at the
first cost whose support is complete, which through ell = 18 is ell + 5.
The point of the second column is that support and positivity are
different obstructions: R_ell can decay to zero even after the support
has filled in, and no analytic statement about the limiting operator
W_3 can see either half.

The row sum is checked at every level against the independent count of
bounded compositions of k with capacities 2, 6, 18, ..., 2*3^(ell-1),
which is Wirsching's identity 2*3^(ell-1)*gbar_ell(k).
"""

from __future__ import annotations

import argparse
import os
import time
from math import isqrt
from pathlib import Path

import numpy as np


def available_ram_bytes() -> int | None:
    """Physical RAM currently available, or None if not determinable
    (e.g. non-Linux). Used only for a preflight warning, not a hard
    cross-platform guarantee."""
    try:
        return os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    except (ValueError, OSError, AttributeError):
        return None


def new_array(shape: tuple[int, int], scratch_dir: Path | None,
              name: str) -> np.ndarray:
    """Allocate the (k_max+1, modulus) table either in RAM (scratch_dir
    is None, the original behaviour, fast but bounded by physical
    memory) or as a disk-backed memmap under scratch_dir (slower, but
    peak physical RAM use is governed by the OS page cache, not by the
    array's full size, so a run whose array does not fit in RAM can
    still complete rather than crash with MemoryError)."""
    if scratch_dir is None:
        array = np.zeros(shape, dtype=np.int64)
    else:
        path = scratch_dir / f"{name}.dat"
        array = np.memmap(path, dtype=np.int64, mode="w+", shape=shape)
        array[:] = 0
    return array


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


def extend(previous: np.ndarray, ell: int, k_max: int,
           scratch_dir: Path | None) -> np.ndarray:
    modulus = 3**ell
    old_modulus = modulus // 3
    current = new_array((k_max + 1, modulus), scratch_dir, f"table_ell{ell}")
    residues = np.arange(old_modulus, dtype=np.int64)
    cap = 2 * old_modulus
    inv2 = pow(2, -1, modulus)
    inverse_power = 1
    for increment in range(min(k_max, cap - 1) + 1):
        inverse_power = (inverse_power * inv2) % modulus
        targets = ((3 * residues + 1) % modulus) * inverse_power % modulus
        for old_cost in range(k_max - increment + 1):
            current[old_cost + increment, targets] += previous[old_cost]
    if scratch_dir is not None:
        current.flush()
    return current


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-ell", type=int, default=14)
    parser.add_argument("--offsets", type=int, nargs="+",
                        default=[0, 4, 5, 6, 7, 8, 10, 12],
                        help="reported costs are ell + offset")
    parser.add_argument("--sqrt-multiples", type=int, nargs="+",
                        default=[1, 2, 3],
                        help="extra columns at ell + u*isqrt(ell)")
    parser.add_argument("--fixed-residues", type=int, nargs="+",
                        default=[1, 2, 4, 5, 7, 8, 10, 11],
                        help="small integers a tracked separately; "
                             "Theorem 1 only consumes integer a")
    parser.add_argument("--fixed-offset", type=int, default=0,
                        help="cost offset used for the fixed-residue column")
    parser.add_argument("--random-fixed", type=int, default=0,
                        help="add this many random integers below "
                             "--random-bound to the tracked set")
    parser.add_argument("--random-bound", type=int, default=3**8)
    parser.add_argument("--seed", type=int, default=20260809)
    parser.add_argument("--scratch-dir", type=str, default=None,
                         help="if set, back the (k_max+1) x 3^ell tables "
                              "with disk-backed memmap files under this "
                              "directory instead of allocating them in "
                              "RAM. Slower (I/O bound instead of memory "
                              "bound) but lets a run whose array does not "
                              "fit in physical RAM complete instead of "
                              "raising MemoryError. Only the current and "
                              "previous level's files are kept on disk at "
                              "once; each is deleted as soon as the next "
                              "level's file is built and flushed.")
    parser.add_argument("--ram-safety-margin", type=float, default=0.9,
                         help="in RAM mode (no --scratch-dir), abort before "
                              "allocating a table whose size exceeds this "
                              "fraction of currently available RAM, rather "
                              "than risk locking up the machine (Rule 9b)")
    args = parser.parse_args()

    scratch_dir = Path(args.scratch_dir) if args.scratch_dir else None
    # Remember whether we created it, so the free-space guard's abort path
    # only removes a directory this run made, never one the caller
    # already had (found in critique round 3).
    scratch_dir_created = False
    if scratch_dir is not None:
        scratch_dir_created = not scratch_dir.exists()
        scratch_dir.mkdir(parents=True, exist_ok=True)

    tracked = list(args.fixed_residues)
    if args.random_fixed:
        rng = np.random.default_rng(args.seed)
        pool = np.arange(1, args.random_bound)
        pool = pool[pool % 3 != 0]
        extra = rng.choice(pool, size=args.random_fixed, replace=False)
        tracked = sorted(set(tracked) | {int(v) for v in extra})

    offsets = sorted(args.offsets)
    multiples = args.sqrt_multiples
    k_max = args.max_ell + max(
        offsets[-1], max(multiples) * isqrt(args.max_ell)
    )

    # extend() builds the new (final) level's array while main()'s loop
    # only deletes the PREVIOUS level's array after extend() returns, so
    # both coexist for the full duration of the last transition: the
    # estimate below is the final level's array plus the previous (1/3
    # the size) level's array, not the final level's array alone.
    # This counts the two big tables only; extend() also holds a few
    # O(8 * 3^(ell-1)) working arrays (residues, targets, the gather),
    # roughly 2/(k_max+1) more, which the RAM margin below absorbs. The
    # disk figure is exact, since only the memmapped tables hit disk.
    final_bytes = (k_max + 1) * 3**args.max_ell * 8
    peak_bytes = final_bytes + final_bytes // 3
    if scratch_dir is None:
        ram = available_ram_bytes()
        if ram is not None and peak_bytes > args.ram_safety_margin * ram:
            raise SystemExit(
                f"Refusing to run in RAM mode: the last level transition "
                f"peaks at about {peak_bytes / 2**30:.1f} GiB (final level "
                f"{final_bytes / 2**30:.1f} GiB plus the still-live "
                f"previous level), more than {args.ram_safety_margin:.0%} "
                f"of the {ram / 2**30:.1f} GiB currently available. Pass "
                f"--scratch-dir to run disk-backed instead, or lower "
                f"--max-ell."
            )
        print(f"[preflight] last transition peaks at ~{peak_bytes / 2**30:.1f} "
              f"GiB; {ram / 2**30:.1f} GiB available now (RAM mode)."
              if ram is not None else
              f"[preflight] last transition peaks at ~{peak_bytes / 2**30:.1f} "
              f"GiB; available RAM could not be determined (RAM mode).")
    else:
        free_disk = None
        try:
            usage = os.statvfs(scratch_dir)
            free_disk = usage.f_bavail * usage.f_frsize
        except OSError:
            pass
        if free_disk is not None and peak_bytes > free_disk:
            # scratch_dir had to exist for statvfs; don't leave an empty
            # directory behind on the abort path, but only remove one this
            # run actually created.
            if scratch_dir_created:
                try:
                    scratch_dir.rmdir()
                except OSError:
                    pass
            raise SystemExit(
                f"Refusing to run disk-backed: the last level transition "
                f"peaks at about {peak_bytes / 2**30:.1f} GiB on disk, "
                f"more than the {free_disk / 2**30:.1f} GiB free at "
                f"{scratch_dir}."
            )
        print(f"[preflight] last transition peaks at ~{peak_bytes / 2**30:.1f} "
              f"GiB on disk at {scratch_dir} (disk-backed mode)"
              + (f", {free_disk / 2**30:.1f} GiB free." if free_disk is not None
                 else ", free disk space could not be determined."))

    table = new_array((k_max + 1, 1), scratch_dir, "table_ell0")
    table[0, 0] = 1
    if scratch_dir is not None:
        table.flush()

    print("min_a g_ell(ell+d, a) / gbar_ell(ell+d); "
          "a dash marks an incomplete support")
    header = ("ell " + " ".join(f"   d={d:<+3d}" for d in offsets)
              + "  |" + " ".join(f"  d={u}rt " for u in multiples)
              + f"  || min over {len(tracked)} fixed integers"
              + f" at d={args.fixed_offset:+d}")
    print(header)
    print("-" * len(header))

    for ell in range(1, args.max_ell + 1):
        started = time.time()
        previous_path = (Path(table.filename)
                          if scratch_dir is not None else None)
        table = extend(table, ell, k_max, scratch_dir)
        if previous_path is not None:
            previous_path.unlink(missing_ok=True)
        modulus = 3**ell
        unit = np.tile(np.array([False, True, True]), modulus // 3)
        n_units = int(unit.sum())

        reference = bounded_composition_counts(ell, k_max)
        fields = []
        columns = list(offsets) + [None] + [u * isqrt(ell) for u in multiples]
        for offset in columns:
            if offset is None:
                fields.append(" |")
                continue
            cost = ell + offset
            if cost > k_max:
                fields.append("      .")
                continue
            row = table[cost][unit]
            total = int(row.sum())
            if total != int(reference[cost]):
                raise AssertionError(
                    f"row total mismatch ell={ell} k={cost}: "
                    f"{total} != {int(reference[cost])}"
                )
            mean = total / n_units
            smallest = int(row.min())
            if smallest == 0:
                fields.append("      -")
            else:
                fields.append(f"{smallest / mean:7.4f}")

        cost = ell + args.fixed_offset
        if 0 <= cost <= k_max:
            mean = int(reference[cost]) / n_units
            picks = [a for a in tracked if a % 3 and a < modulus]
            values = [int(table[cost, a]) for a in picks]
            worst = min(values) / mean if picks and mean else float("nan")
            report = f"  || {worst:8.4f} ({len(picks):3d})"
        else:
            report = "  ||        ."

        print(f"{ell:3d} " + " ".join(fields) + report
              + f"   [{time.time() - started:6.1f}s]", flush=True)


if __name__ == "__main__":
    main()
