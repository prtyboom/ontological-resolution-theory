"""
SPECTRAL TORUS TEST
===================
Проверка: появляется ли (pi-3) на торе T^2?
Сравниваем непрерывный тор с дискретной решёткой N x N

ORT Project | June 2025
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import os

# Гарантируем существование папки
os.makedirs(r'C:\Users\admin\Documents\ontological-resolution-theory\figures', exist_ok=True)

PI_MINUS_3 = np.pi - 3


def dirac_continuous_torus(N_modes):
    """
    Спектр |D| на непрерывном торе T^2 с антипериодическими условиями.
    λ_{k,m} = sqrt((k+1/2)^2 + (m+1/2)^2), k,m ∈ [-N_modes, N_modes]
    Возвращаем |λ| (для действия e^{-λ²/Λ²} знак не важен).
    """
    spectrum = []
    for k in range(-N_modes, N_modes + 1):
        for m in range(-N_modes, N_modes + 1):
            lam = np.sqrt((k + 0.5)**2 + (m + 0.5)**2)
            spectrum.append(lam)
    return np.sort(spectrum)


def shift_1d_antiperiodic(N):
    """1D shift matrix with antiperiodic BC: ψ_N = -ψ_0"""
    S = np.roll(np.eye(N), -1, axis=1)
    S[-1, 0] = -1
    return S


def dirac_discrete_torus(N):
    """
    Correct Dirac operator on N×N torus with antiperiodic BC.
    Uses tensor product of 1D shifts.
    """
    a = 2 * np.pi / N  # lattice spacing
    N2 = N * N

    # 1D shifts with twist
    Sx_1d = shift_1d_antiperiodic(N)
    Sy_1d = shift_1d_antiperiodic(N)

    # 2D shifts
    Sx = np.kron(Sx_1d, np.eye(N))
    Sy = np.kron(np.eye(N), Sy_1d)

    # Discrete derivatives: ∇_x = (Sx - Sx^T) / (2a), same for y
    Dx = (Sx - Sx.T) / (2 * a)
    Dy = (Sy - Sy.T) / (2 * a)

    # Dirac operator: D = [[0, Dx - i Dy], [Dx + i Dy, 0]]
    zero = np.zeros((N2, N2), dtype=complex)
    D_plus = Dx + 1j * Dy
    D_minus = Dx - 1j * Dy

    D_full = np.block([[zero, D_minus],
                       [D_plus, zero]])
    return D_full


def spectral_action(spectrum, Lambda):
    """S = Tr(exp(-D^2 / Lambda^2)) = sum exp(-λ² / Λ²)"""
    return np.sum(np.exp(-spectrum**2 / Lambda**2))


def compute_defect_torus(N, Lambda):
    """
    Computes defect = S_cont - S_disc for T^2 vs C_N × C_N
    Uses consistent number of modes: (2N+1)^2 for continuous ≈ 4N²
    """
    # Discrete spectrum
    D_disc = dirac_discrete_torus(N)
    D2 = D_disc @ D_disc.conj().T  # D is Hermitian → D2 = D^2
    eigvals_D2 = linalg.eigvalsh(D2)
    spectrum_disc = np.sqrt(np.maximum(eigvals_D2, 0.0))

    # Continuous: take enough modes to cover similar range
    N_modes = N  # (2N+1)^2 ≈ 4N², same order as disc
    spectrum_cont = dirac_continuous_torus(N_modes)

    # To compare fairly, use same number of modes (or cap continuous)
    min_len = min(len(spectrum_cont), len(spectrum_disc))
    spec_c = spectrum_cont[:min_len]
    spec_d = spectrum_disc[:min_len]

    S_cont = spectral_action(spec_c, Lambda)
    S_disc = spectral_action(spec_d, Lambda)

    return S_cont - S_disc, S_cont, S_disc, spec_d, spec_c


# === MAIN EXECUTION ===
print("=" * 70)
print("SPECTRAL TORUS TEST: Поиск (pi-3) на T^2")
print(f"pi - 3 = {PI_MINUS_3}")
print("=" * 70)

N_list = [3, 4, 5, 6, 7, 8, 10, 12]
Lambda = 10.0

print(f"\nLambda = {Lambda}")
print("\n  N     N^2    |Defect|   |Def|/(pi-3)   |Def|/(4*N^2)")
print("-" * 65)

results = []
for N in N_list:
    try:
        defect, S_cont, S_disc, _, _ = compute_defect_torus(N, Lambda)
        abs_def = abs(defect)
        ratio = abs_def / PI_MINUS_3
        norm4 = abs_def / (4 * N * N)
        results.append((N, abs_def, ratio, norm4))
        print(f"{N:3d}   {N*N:4d}   {abs_def:9.4f}   {ratio:11.4f}   {norm4:12.6f}")
    except Exception as e:
        print(f"{N:3d}   ERROR: {e}")

# FOCUS: best divisor search
print("\n" + "=" * 70)
print("FOCUS: Поиск оптимального делителя")
print("=" * 70)

for N in [4, 6, 8, 12]:
    try:
        defect, _, _, _, _ = compute_defect_torus(N, Lambda)
        d = abs(defect)
        print(f"\nN = {N} (N^2 = {N*N}):")

        best_div, best_match = 1, 0
        for div in range(1, 500):
            candidate = d / div
            match = 100 * (1 - abs(candidate - PI_MINUS_3) / PI_MINUS_3)
            if match > best_match:
                best_match = match
                best_div = div

        print(f"  |Defect| = {d:.6f}")
        print(f"  Best divisor: {best_div}")
        print(f"  |Defect|/{best_div} = {d/best_div:.6f}")
        print(f"  Match to (pi-3): {best_match:.2f}%")

        if best_match > 95:
            # Factorization
            n = best_div
            factors = []
            temp = n
            p = 2
            while p * p <= temp:
                while temp % p == 0:
                    factors.append(p)
                    temp //= p
                p += 1
            if temp > 1:
                factors.append(temp)
            print(f"  Factorization: {best_div} = {' × '.join(map(str, factors))}")

    except Exception as e:
        print(f"  N={N} failed: {e}")

# Lambda scan for N=4
print("\n" + "=" * 70)
print("LAMBDA SCAN для N=4")
print("=" * 70)

N_test = 4
Lambda_list = [2, 3, 4, 5, 6, 8, 10, 12, 15, 20]

print("\n Lambda   |Defect|   Best Div   Match%")
print("-" * 45)

for L in Lambda_list:
    try:
        d, _, _, _, _ = compute_defect_torus(N_test, L)
        d = abs(d)
        best_div, best_match = 1, 0
        for div in range(1, 200):
            match = 100 * (1 - abs(d/div - PI_MINUS_3) / PI_MINUS_3)
            if match > best_match:
                best_match = match
                best_div = div
        print(f"  {L:4.0f}    {d:8.4f}      {best_div:4d}    {best_match:5.1f}%")
    except:
        print(f"  {L:4.0f}    ERROR")

# Plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Defect ratio vs N²
if results:
    Ns = [r[0] for r in results]
    N2s = [r[0]**2 for r in results]
    ratios = [r[2] for r in results]
    axes[0].plot(N2s, ratios, 'bo-', lw=2, ms=8)
    axes[0].set_xlabel('N² (число узлов)', fontsize=12)
    axes[0].set_ylabel('|Defect| / (π − 3)', fontsize=12)
    axes[0].set_title('Тор T²: дефект vs размер решётки', fontsize=14)
    axes[0].grid(True, alpha=0.3)

# Comparison with S¹ (approximate values from your run)
circle_N = [3, 4, 6, 8, 12]
circle_ratio = [1.1, 2.6, 8.5, 18.8, 51.7]
axes[1].plot(circle_N, circle_ratio, 'bs--', lw=2, ms=8, label='Circle S¹')
if results:
    axes[1].plot(Ns, ratios, 'ro-', lw=2, ms=8, label='Torus T²')
axes[1].set_xlabel('N', fontsize=12)
axes[1].set_ylabel('|Defect| / (π − 3)', fontsize=12)
axes[1].set_title('Сравнение: окружность vs тор', fontsize=14)
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
output_path = r'C:\Users\admin\Documents\ontological-resolution-theory\figures\spectral_torus.png'
plt.savefig(output_path, dpi=150)
print(f"\nSaved: {output_path}")
plt.show()

print("\n" + "=" * 70)
print("TORUS TEST COMPLETE")
print("=" * 70)