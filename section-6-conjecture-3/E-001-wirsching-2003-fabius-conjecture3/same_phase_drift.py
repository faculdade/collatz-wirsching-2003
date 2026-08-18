#!/usr/bin/env python3
"""The same-phase drift of ln(phi/phi_0), for the same-phase remark of paper 05.

the same-phase remark compares ln(phi/phi_0) at the pairs (l, 3l) evaluated at the
BARE points x_l = l/3^l, not at the shifted x_l^+ the rest of the paper
uses. The pair shares phase exactly: x_{3l}/x_l = 3^(1-2l), an integer
shift in log_3, so a 1-periodic Q(log_3 x) cancels identically from the
difference at any amplitude. The drift therefore measures the finite-size
correction to phi_0 and nothing else.

Reported for each pair:
  observed  = ln(phi/phi_0)(x_{3l}) - ln(phi/phi_0)(x_l)
  predicted = b (1/sqrt(3l) - 1/sqrt(l)) = b (3^{-1/2} - 1)/sqrt(l)
              under the smooth C/sqrt(l) fit, with b read from the
              committed shifted-point run
  mismatch  = |observed - predicted| / |predicted|

Written 2026-08-17: critique round 3 found these four observed values,
four predicted values and four mismatch percentages quoted in the paper
with no committed output behind any of them. Three of the observed values
were recoverable from an older bare-point log; the (150,450) pair was not,
because that run stopped at l=300.

Run:  python3 same_phase_drift.py [--b -0.7916]
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from fractions import Fraction as Fr
from math import sqrt


def load_main():
    spec = importlib.util.spec_from_file_location(
        "e", "experiment_conjecture3.py")
    mod = importlib.util.module_from_spec(spec)
    saved, sys.argv = sys.argv, ["x"]
    spec.loader.exec_module(mod)
    sys.argv = saved
    return mod


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--b", type=float, default=-0.7916,
                    help="slope of the C/sqrt(l) fit, from the committed "
                         "shifted-point run")
    args = ap.parse_args()

    m = load_main()
    pairs = [(10, 30), (50, 150), (100, 300), (150, 450)]
    N = max(b for _, b in pairs) + 10
    print(f"computing exact moments to M_{N} ...", flush=True)
    M = m.moments(N)

    def lnr(ell: int) -> float:
        """ln(phi/phi_0) at the BARE point x_l = l/3^l."""
        x = Fr(ell, 3 ** ell)
        logphi, rel, _, _, _ = m.log_phi_reduced(M, x)
        assert rel == 0, f"nonzero truncation error at l={ell}: {rel}"
        return float(logphi - m.log_phi0(x))

    print(f"\nsame-phase drift at the bare points, b = {args.b}\n")
    print(f"{'(l,3l)':>12} {'observed':>10} {'predicted':>10} {'mismatch':>10}")
    for a, b in pairs:
        obs = lnr(b) - lnr(a)
        pred = args.b * (3 ** -0.5 - 1) / sqrt(a)
        mis = abs(obs - pred) / abs(pred)
        print(f"{'(%d,%d)' % (a, b):>12} {obs:+10.4f} {pred:+10.4f} "
              f"{100 * mis:9.0f}%")
    print("\nA 1-periodic Q cancels from every row above by construction, so")
    print("these measure the finite-size correction to phi_0 alone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
