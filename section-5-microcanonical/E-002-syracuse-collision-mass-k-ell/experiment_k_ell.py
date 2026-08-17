#!/usr/bin/env python3
"""
Checagem finita (pre-escrita, pedida pelo advisor antes de formalizar o
Lema 2' de H-126): como se comporta K_ell nos niveis ell <= 17?

Nenhum intervalo finito decide se K_inf e' finito. A saida deste script
serve apenas como evidencia numerica; ver H-140 para a correcao de status.

K_ell := 3^ell * sum_y p_ell(y)^2 = 3^ell * P(Syrac_ell = Syrac'_ell)
para copias independentes.

Recursao usada (derivada nesta sessao a partir da identidade de Tao
F_n(a) = 2^{-a_1}*(1+3*F_{n-1}(a_2,...,a_n)), a_1~Geom(2) memoryless):

  mu_n(y) = (1/2)*nu(2y mod 3^n) + (1/2)*mu_n(2y mod 3^n)

onde nu = lei de 1+3*F_{n-1} (suportada nos residuos = 1 mod 3).
F_n e' sempre invertivel mod 3 (nunca multiplo de 3, Remark 1.15 de
Tao), e 2 e' raiz primitiva mod 3^n para todo n (fato classico), entao
o suporte inteiro de mu_n vive numa unica orbita ciclica de tamanho
P=2*3^(n-1) sob multiplicacao por 2. Ao longo da orbita (indexada por
k, y_k = 2^k mod 3^n), a equacao acima e' uma recursao linear circular
que se resolve por serie geometrica truncada (kernel (1/2)^s decai tao
rapido que s>100 e' irrelevante em ponto flutuante).
"""
import numpy as np
import time


def build_nu(p_prev, ell):
    mod_prev = len(p_prev)
    mod = 3 ** ell
    nu = np.zeros(mod)
    idx = 1 + 3 * np.arange(mod_prev)
    nu[idx] = p_prev
    return nu


def solve_level(p_prev, ell, s_max=100):
    mod = 3 ** ell
    P = 2 * (mod // 3)
    nu = build_nu(p_prev, ell)

    y = np.empty(P, dtype=np.int64)
    y[0] = 1
    for k in range(1, P):
        y[k] = (y[k - 1] * 2) % mod
    n_k = nu[y]

    x = np.zeros(P)
    w = 1.0
    for s in range(s_max):
        w = 0.5 ** (s + 1)
        if w == 0.0:
            break
        x += w * np.roll(n_k, -(s + 1))

    p = np.zeros(mod)
    p[y] = x
    return p


def main():
    p = np.zeros(3)
    p[1] = 1 / 3
    p[2] = 2 / 3
    assert abs(p.sum() - 1.0) < 1e-12

    K1 = 3 * np.sum(p ** 2)
    print(f"ell=1  sum(p)={p.sum():.10f}  K_ell={K1:.6f}  (esperado exato 5/3={5/3:.6f})")

    results = [(1, K1)]
    for ell in range(2, 18):
        t0 = time.time()
        p = solve_level(p, ell)
        K = (3 ** ell) * np.sum(p ** 2)
        dt = time.time() - t0
        print(f"ell={ell:2d}  sum(p)={p.sum():.8f}  K_ell={K:10.5f}  "
              f"K_ell/ell={K/ell:8.5f}  tempo={dt:6.2f}s", flush=True)
        results.append((ell, K))

    print("\n=== Diagnostico ===")
    print("Se K_ell estabiliza/satura -> Hipotese L2 (K_inf<inf) suportada -> Lema 2' de H-126 tem conteudo.")
    print("Se K_ell cresce no intervalo testado -> evidencia finita contra L2; nao prova divergencia.")
    ells = np.array([r[0] for r in results[-6:]])
    Ks = np.array([r[1] for r in results[-6:]])
    # ajuste log-log tosco para ver se K cresce como ell^p
    if len(ells) >= 4 and np.all(Ks > 0):
        coeffs = np.polyfit(np.log(ells), np.log(Ks), 1)
        print(f"Ajuste K_ell ~ ell^p nos ultimos pontos: p={coeffs[0]:.3f}")


if __name__ == "__main__":
    main()


# Migration note: references H-126, H-140, out of this repo's scope, left unrenumbered.
