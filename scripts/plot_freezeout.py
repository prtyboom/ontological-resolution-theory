# scripts/plot_freezeout.py
"""
График отношения Gamma_grav / H для чисто гравитационного теплообмена
в 5D архивной модели при eps = 1e-70.
"""

import numpy as np
import matplotlib.pyplot as plt

from archival_dm.theory import gamma_over_H_grav


def main():
    eps = 1e-70

    # Диапазон температур: от 1 ГэВ до M_P ~ 1e19 ГэВ
    T_min = 1.0       # GeV
    T_max = 2.0e19    # GeV
    T_vals = np.logspace(np.log10(T_min), np.log10(T_max), 200)

    ratio = gamma_over_H_grav(T_vals, eps)

    plt.figure(figsize=(6, 4))
    plt.loglog(T_vals, ratio, label=r'$\Gamma_{\rm grav}/H$')
    plt.axhline(1.0, color='r', ls='--', label=r'$\Gamma/H = 1$')

    plt.xlabel(r'$T\ \mathrm{[GeV]}$')
    plt.ylabel(r'$\Gamma_{\rm grav}/H$')
    plt.title(r'Gravitational heat exchange, $\varepsilon=10^{-70}$')
    plt.legend()
    plt.grid(True, which='both', ls=':')

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()