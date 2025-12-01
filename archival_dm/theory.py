# archival_dm/theory.py
"""
Analytical utilities for the 5D archival dark matter model.

Содержит:
- V_eff(eps) для степенного warp-профиля Phi(z) = k * |z|^p
- параметрическую формулу для Omega_DM/Omega_b
- обратную задачу: alpha из заданного Omega_DM/Omega_b
- оценку отношения Gamma_grav / H (заморозка при чисто гравитационном обмене)
"""

import math
import numpy as np

# Reduced Planck mass (в ГэВ), если понадобится в оценках
MPL = 2.435e18

# Стандартное g_* для SM при высоких T
GSTAR_SM_DEFAULT = 106.75


def veff_powerlaw(eps, k, p):
    """
    V_eff(eps) для warp-профиля Phi(z) = k * |z|^p.

    Определение:
      V_eff(eps) = ∫_{-∞}^{+∞} dz exp[-4 * Phi(z) / eps]
                  = ∫ dz exp[-4 k |z|^p / eps]

    При p > 0 получаем аналитически:
      V_eff(eps) = (2/p) * (eps / (4k))^(1/p) * Gamma(1/p)

    eps : float
        Параметр подавления (ε).
    k   : float
        Масштаб варп-фактора.
    p   : float
        Показатель степени в Phi(z) = k |z|^p.

    Возвращает:
      float : значение V_eff(eps) (в единицах 1/k при k>0).
    """
    eps = float(eps)
    k = float(k)
    p = float(p)

    if k <= 0.0:
        raise ValueError("k must be positive for veff_powerlaw.")
    if p <= 0.0:
        raise ValueError("p must be positive for veff_powerlaw.")

    base = eps / (4.0 * k)
    exponent = 1.0 / p
    pref = 2.0 / p
    # math.gamma для вещественного аргумента > 0
    return pref * (base ** exponent) * math.gamma(exponent)


def omega_dm_over_ob(
    N,
    eps,
    alpha,
    gstar_sm=GSTAR_SM_DEFAULT,
    M5_over_Mpl=1.0,
    beta=1.0,
    S_entropy=0.31,
):
    """
    Параметрическая формула для Omega_DM/Omega_b:

      Omega_DM / Omega_b ~ (N * eps^alpha / g*_SM)
                           * (M5/Mpl)^beta
                           * S_entropy

    Здесь это АНЗАЦ, а не строгий вывод из полной 5D-теории.

    Параметры:
      N          : голографическое число степеней свободы (~1e122)
      eps        : подавление четвертого измерения (~1e-70)
      alpha      : показатель, связанный с геометрией (V_eff ~ eps^alpha)
      gstar_sm   : эффективное число степеней свободы SM при T_dec
      M5_over_Mpl: отношение 5D и 4D планковских масс (M5/Mpl)
      beta       : показатель для scale-фактора (обычно 0–2)
      S_entropy  : энтропийный фактор ~(g_*^SM,dec / g_*^SM,0)^(-1/3)

    Возвращает:
      float : оценка отношения Omega_DM / Omega_b.
    """
    N = float(N)
    eps = float(eps)
    alpha = float(alpha)
    gstar_sm = float(gstar_sm)
    M5_over_Mpl = float(M5_over_Mpl)
    beta = float(beta)
    S_entropy = float(S_entropy)

    factor_geom = N * (eps ** alpha) / gstar_sm
    factor_scale = M5_over_Mpl ** beta
    ratio = factor_geom * factor_scale * S_entropy
    return ratio


def alpha_from_ratio(
    ratio_target,
    N,
    eps,
    gstar_sm=GSTAR_SM_DEFAULT,
    M5_over_Mpl=1.0,
    beta=1.0,
    S_entropy=0.31,
):
    """
    Обратная задача: по заданному Omega_DM/Omega_b (= ratio_target)
    найти alpha при фиксированных N, eps, g*_SM, (M5/Mpl)^beta, S_entropy.

    Решаем:
      ratio_target = (N * eps^alpha / g*_SM) * (M5/Mpl)^beta * S_entropy
      => eps^alpha = ratio_target * g*_SM / (N * (M5/Mpl)^beta * S_entropy)
      => alpha     = log10( RHS ) / log10(eps)
    """
    ratio_target = float(ratio_target)
    N = float(N)
    eps = float(eps)
    gstar_sm = float(gstar_sm)
    M5_over_Mpl = float(M5_over_Mpl)
    beta = float(beta)
    S_entropy = float(S_entropy)

    rhs = (
        ratio_target
        * gstar_sm
        / (N * (M5_over_Mpl ** beta) * S_entropy)
    )
    if rhs <= 0.0:
        raise ValueError(
            f"RHS <= 0 in alpha_from_ratio (rhs={rhs}), "
            "невозможно взять логарифм."
        )

    log10_rhs = math.log10(rhs)
    log10_eps = math.log10(eps)  # eps ~ 1e-70 => log10_eps ~ -70

    alpha = log10_rhs / log10_eps
    return alpha


def gamma_over_H_grav(T, eps, Mpl=MPL):
    """
    Оценка отношения Gamma_grav / H для чисто гравитационного
    теплообмена между архивным и SM-сектором.

    Приблизительно:
      Gamma_grav ~ (T^5 / Mpl^4) * eps^2
      H          ~ T^2 / Mpl      (радиационно-доминированная эпоха)

    => Gamma_grav / H ~ (T^3 / Mpl^3) * eps^2

    Параметры:
      T   : температура (в ГэВ), скаляр или numpy-массив
      eps : подавление ε
      Mpl : планковская масса (reduced Planck mass), по умолчанию 2.435e18 ГэВ

    Возвращает:
      numpy.ndarray или float : Gamma_grav / H как функция T.
    """
    T = np.array(T, dtype=float)
    ratio = (T ** 3) / (Mpl ** 3) * (eps ** 2)
    return ratio