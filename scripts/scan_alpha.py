# scripts/scan_alpha.py
"""
Скрипт для вычисления alpha в 5D архивной модели
по заданному наблюдаемому отношению Omega_DM / Omega_b.
"""

from archival_dm.theory import alpha_from_ratio, omega_dm_over_ob


def main():
    # Входные параметры модели
    N = 1e122
    eps = 1e-70
    gstar_sm = 106.75
    M5_over_Mpl = 1.0
    beta = 1.0
    S_entropy = 0.31
    ratio_target = 5.4  # наблюдаемое Omega_DM / Omega_b

    alpha = alpha_from_ratio(
        ratio_target=ratio_target,
        N=N,
        eps=eps,
        gstar_sm=gstar_sm,
        M5_over_Mpl=M5_over_Mpl,
        beta=beta,
        S_entropy=S_entropy,
    )

    print("=== Archival 5D model: solving for alpha ===")
    print(f"N           = {N:.2e}")
    print(f"eps         = {eps:.2e}")
    print(f"g*_SM       = {gstar_sm}")
    print(f"(M5/Mpl)^β  = {M5_over_Mpl**beta:.2e}")
    print(f"S_entropy   = {S_entropy}")
    print(f"target Ω_DM/Ω_b = {ratio_target}")

    print(f"\n=> alpha ≈ {alpha:.3f}")

    # Проверка: подставим alpha обратно в формулу
    ratio_check = omega_dm_over_ob(
        N=N,
        eps=eps,
        alpha=alpha,
        gstar_sm=gstar_sm,
        M5_over_Mpl=M5_over_Mpl,
        beta=beta,
        S_entropy=S_entropy,
    )

    print(f"Check: Ω_DM/Ω_b(alpha) = {ratio_check:.3f}")


if __name__ == "__main__":
    main()