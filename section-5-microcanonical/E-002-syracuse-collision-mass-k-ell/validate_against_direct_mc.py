#!/usr/bin/env python3
"""
Validacao pedida pelo advisor: a recursao de solve_level() em
k_ell_check.py so foi conferida no caso base ell=1 (que nao exercita
a maquinaria). Aqui comparamos p_dp (da recursao) contra uma
amostragem Monte Carlo DIRETA de Syrac_ell (definicao primaria,
F_ell(a) = sum_j 3^(j-1) * 2^{-a_[1,j]} mod 3^ell, a_i~Geom(2) iid),
para ell=3 e ell=4, bin a bin (nao so K_ell).
"""
import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from experiment_k_ell import solve_level


def direct_mc_syrac(ell, n_samples, rng):
    mod = 3 ** ell
    counts = np.zeros(mod, dtype=np.int64)
    inv2 = pow(2, -1, mod)  # 2^{-1} mod 3^ell
    batch = 200_000
    done = 0
    while done < n_samples:
        b = min(batch, n_samples - done)
        # Geom(2): P(a=k)=2^-k, k>=1  => a = 1 + floor(log2(1/U)) equivalente a
        # amostrar via geometric distribution numpy (p=0.5, valores >=0) +1
        a = rng.geometric(p=0.5, size=(b, ell)).astype(np.int64)  # numpy geometric: k>=1, P(k)=p(1-p)^{k-1} com p=0.5 -> igual a Geom(2) do paper
        cumsum_a = np.cumsum(a, axis=1)  # a_[1,j] para j=1..ell
        F = np.zeros(b, dtype=object)  # usar Python int (mod grande, mas ell pequeno aqui, cabe em int64 tambem)
        F = np.zeros(b, dtype=np.int64)
        pow3 = 1
        for j in range(ell):
            # 2^{-a_[1,j]} mod mod = inv2^{a_[1,j]} mod mod
            exps = cumsum_a[:, j]
            # calcular inv2^exps mod mod vetorizado via power modular rapido (loop pequeno, ell<=4)
            term = pow(inv2, 1, mod)
            vals = np.array([pow(int(inv2), int(e), mod) for e in exps], dtype=np.int64)
            F = (F + pow3 * vals) % mod
            pow3 = (pow3 * 3) % mod
        idx, c = np.unique(F, return_counts=True)
        counts[idx] += c
        done += b
    return counts / counts.sum()


def main():
    rng = np.random.default_rng(123)
    p1 = np.zeros(3)
    p1[1] = 1 / 3
    p1[2] = 2 / 3

    p = p1
    for ell in [2, 3, 4]:
        p = solve_level(p, ell)
        if ell in (3, 4):
            n = 3_000_000 if ell == 3 else 3_000_000
            p_mc = direct_mc_syrac(ell, n, rng)
            maxdiff = np.max(np.abs(p - p_mc))
            print(f"ell={ell}: max|p_dp - p_mc| = {maxdiff:.6f}  "
                  f"(ruido MC esperado ~ 1/sqrt(N)~{1/np.sqrt(n):.6f})")
            # imprime alguns valores lado a lado
            support = np.where(p > 1e-9)[0][:6]
            for y in support:
                print(f"    y={y:4d}  p_dp={p[y]:.6f}  p_mc={p_mc[y]:.6f}")


if __name__ == "__main__":
    main()
