"""
OPERATION SPECTRAL GAP v2
=========================
Поиск (pi - 3) в относительном спектральном дефекте
ORT Project | June 2025
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt


def dirac_spectrum_continuous(N_modes):
    """
    Спектр Дирака на S^1, длина L=2*pi
    lambda_k = k + 1/2, где k от -N до N
    """
    k = np.arange(-N_modes, N_modes + 1)
    return k + 0.5


def dirac_matrix_graph(N):
    """
    Оператор Дирака на цикле C_N
    Антипериодические граничные условия
    """
    a = 2 * np.pi / N
    
    D = np.zeros((N, N), dtype=complex)
    
    for j in range(N):
        j_right = (j + 1) % N
        j_left = (j - 1) % N
        
        phase_right = -1 if j == N - 1 else 1
        phase_left = -1 if j == 0 else 1
        
        D[j, j_right] = phase_right / (2 * a)
        D[j, j_left] = -phase_left / (2 * a)
    
    return D


def heat_trace(eigenvalues, t):
    """
    Heat kernel trace: Tr(exp(-t*D^2))
    """
    return np.sum(np.exp(-t * eigenvalues**2))


def compute_all(N, t_values):
    """
    Вычисляет спектры и heat traces
    """
    # Дискретный спектр
    D = dirac_matrix_graph(N)
    eigs_disc = np.real(linalg.eigvals(D))
    eigs_disc = np.sort(eigs_disc)
    
    # Непрерывный спектр (столько же мод)
    eigs_cont = dirac_spectrum_continuous(N // 2)
    
    # Heat traces для разных t
    H_disc = [heat_trace(eigs_disc, t) for t in t_values]
    H_cont = [heat_trace(eigs_cont, t) for t in t_values]
    
    return eigs_disc, eigs_cont, np.array(H_disc), np.array(H_cont)


def main():
    PI_MINUS_3 = np.pi - 3
    
    print("=" * 60)
    print("OPERATION SPECTRAL GAP v2")
    print("Target: (pi - 3) = {:.8f}".format(PI_MINUS_3))
    print("=" * 60)
    
    # Параметры
    N_list = [6, 8, 10, 12, 14, 16, 20, 24, 30, 36, 48, 60]
    t_values = np.linspace(0.01, 2.0, 100)
    
    results = []
    
    print("\n--- Heat Trace Analysis ---\n")
    print("{:>5} {:>15} {:>15} {:>15}".format(
        "N", "H_cont(t=1)", "H_disc(t=1)", "Relative Δ"))
    print("-" * 55)
    
    for N in N_list:
        eigs_d, eigs_c, H_d, H_c = compute_all(N, t_values)
        
        # Значения при t=1
        idx_t1 = np.argmin(np.abs(t_values - 1.0))
        H_cont_1 = H_c[idx_t1]
        H_disc_1 = H_d[idx_t1]
        
        # Относительный дефект
        if H_cont_1 != 0:
            rel_defect = (H_cont_1 - H_disc_1) / H_cont_1
        else:
            rel_defect = 0
        
        results.append({
            'N': N,
            'H_cont': H_cont_1,
            'H_disc': H_disc_1,
            'rel_defect': rel_defect,
            'H_cont_full': H_c,
            'H_disc_full': H_d
        })
        
        print("{:5d} {:15.6f} {:15.6f} {:15.8f}".format(
            N, H_cont_1, H_disc_1, rel_defect))
    
    # Фокус на N=12
    print("\n" + "=" * 60)
    print("FOCUS: N = 12")
    print("=" * 60)
    
    r12 = [r for r in results if r['N'] == 12][0]
    
    print(f"\nRelative defect at N=12: {r12['rel_defect']:.8f}")
    print(f"(pi - 3) =                {PI_MINUS_3:.8f}")
    print(f"Ratio:                    {r12['rel_defect'] / PI_MINUS_3:.6f}")
    
    # Поиск t где дефект = pi-3
    print("\n--- Поиск t где |Δ| ≈ (π-3) ---")
    
    eigs_d, eigs_c, H_d, H_c = compute_all(12, t_values)
    rel_defects = (H_c - H_d) / np.where(H_c != 0, H_c, 1)
    
    # Где abs(rel_defect) ближе всего к pi-3?
    distances = np.abs(np.abs(rel_defects) - PI_MINUS_3)
    best_idx = np.argmin(distances)
    best_t = t_values[best_idx]
    best_defect = rel_defects[best_idx]
    
    print(f"Best match at t = {best_t:.4f}")
    print(f"Relative defect = {best_defect:.8f}")
    print(f"|(π-3)| = {PI_MINUS_3:.8f}")
    print(f"Difference: {abs(abs(best_defect) - PI_MINUS_3):.2e}")
    
    # Графики
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Относительный дефект vs N при t=1
    ax1 = axes[0, 0]
    Ns = [r['N'] for r in results]
    rel_defs = [r['rel_defect'] for r in results]
    ax1.plot(Ns, rel_defs, 'bo-', lw=2, ms=8)
    ax1.axhline(PI_MINUS_3, color='r', ls='--', lw=2, label=f'π-3 = {PI_MINUS_3:.4f}')
    ax1.axhline(-PI_MINUS_3, color='r', ls='--', lw=2)
    ax1.axhline(0, color='gray', ls='-', lw=1)
    ax1.axvline(12, color='gold', ls=':', lw=2, label='N=12')
    ax1.set_xlabel('N', fontsize=12)
    ax1.set_ylabel('Relative Defect (H_c - H_d)/H_c', fontsize=12)
    ax1.set_title('Относительный дефект vs N (t=1)', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Heat trace для N=12 vs t
    ax2 = axes[0, 1]
    ax2.plot(t_values, H_c, 'b-', lw=2, label='Continuous S^1')
    ax2.plot(t_values, H_d, 'r--', lw=2, label='Discrete C_12')
    ax2.set_xlabel('t', fontsize=12)
    ax2.set_ylabel('Heat Trace', fontsize=12)
    ax2.set_title('Heat Kernel Trace: N=12', fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Относительный дефект vs t для N=12
    ax3 = axes[1, 0]
    ax3.plot(t_values, rel_defects, 'g-', lw=2)
    ax3.axhline(PI_MINUS_3, color='r', ls='--', lw=2, label=f'π-3')
    ax3.axhline(-PI_MINUS_3, color='r', ls='--', lw=2)
    ax3.axvline(best_t, color='gold', ls=':', lw=2, label=f't*={best_t:.3f}')
    ax3.set_xlabel('t', fontsize=12)
    ax3.set_ylabel('Relative Defect', fontsize=12)
    ax3.set_title('Rel. Defect vs t (N=12) — поиск π-3', fontsize=14)
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Спектры
    ax4 = axes[1, 1]
    eigs_d_12, eigs_c_12, _, _ = compute_all(12, [1.0])
    ax4.stem(range(len(eigs_c_12)), np.sort(eigs_c_12), 'b', markerfmt='bo', 
             label='Continuous', basefmt=' ')
    ax4.stem(range(len(eigs_d_12)), np.sort(eigs_d_12), 'r', markerfmt='rx',
             label='Discrete', basefmt=' ')
    ax4.set_xlabel('Mode index', fontsize=12)
    ax4.set_ylabel('Eigenvalue λ', fontsize=12)
    ax4.set_title('Спектры Дирака: N=12', fontsize=14)
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    out_path = r'C:\Users\admin\Documents\ontological-resolution-theory\figures\spectral_defect_v2.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"\nGraph saved: {out_path}")
    
    plt.show()
    
    print("\n" + "=" * 60)
    print("OPERATION COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()