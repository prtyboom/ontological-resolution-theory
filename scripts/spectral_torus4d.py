"""
SPECTRAL TORUS 4D TEST
=======================
Проверка: появляется ли (pi - 3) на T^4?
ORT Project | December 2025
"""

import numpy as np
from scipy import linalg
import os

PI_MINUS_3 = np.pi - 3

def shift_1d_antiperiodic(N):
    """1D shift with antiperiodic BC: ψ_N = -ψ_0"""
    S = np.roll(np.eye(N), -1, axis=1)
    S[-1, 0] = -1
    return S

def dirac_continuous_torus4d(N_modes):
    """Continuous |D| spectrum on T^4 with antiperiodic BC"""
    spectrum = []
    for k in range(-N_modes, N_modes + 1):
        for m in range(-N_modes, N_modes + 1):
            for n in range(-N_modes, N_modes + 1):
                for p in range(-N_modes, N_modes + 1):
                    lam = np.sqrt((k + 0.5)**2 + (m + 0.5)**2 + (n + 0.5)**2 + (p + 0.5)**2)
                    spectrum.append(lam)
    return np.sort(spectrum)

def dirac_discrete_torus4d(N):
    """Dirac operator on N^4 lattice with antiperiodic BC"""
    a = 2 * np.pi / N
    N4 = N**4

    S1 = shift_1d_antiperiodic(N)
    I = np.eye(N)

    # 4D shifts via kron
    Sx = np.kron(np.kron(np.kron(S1, I), I), I)
    Sy = np.kron(np.kron(np.kron(I, S1), I), I)
    Sz = np.kron(np.kron(np.kron(I, I), S1), I)
    St = np.kron(np.kron(np.kron(I, I), I), S1)

    # Derivatives
    Dx = (Sx - Sx.T) / (2 * a)
    Dy = (Sy - Sy.T) / (2 * a)
    Dz = (Sz - Sz.T) / (2 * a)
    Dt = (St - St.T) / (2 * a)

    # Chiral coupling (simplified 2-component)
    zero = np.zeros((N4, N4), dtype=complex)
    D_plus = Dx + 1j * (Dy + Dz + Dt)
    D_minus = Dx - 1j * (Dy + Dz + Dt)

    return np.block([[zero, D_minus],
                     [D_plus, zero]])

def spectral_action(spectrum, Lambda):
    return np.sum(np.exp(-spectrum**2 / Lambda**2))

def compute_defect_torus4d(N, Lambda):
    # Discrete spectrum
    D = dirac_discrete_torus4d(N)
    D2 = D @ D.conj().T
    eigs = np.sqrt(np.maximum(np.real(linalg.eigvalsh(D2)), 0.0))

    # Continuous spectrum
    N_modes = N  # (2N+1)^4 ~ 16N^4, matching order
    spec_cont = dirac_continuous_torus4d(N_modes)

    min_len = min(len(eigs), len(spec_cont))
    S_disc = spectral_action(eigs[:min_len], Lambda)
    S_cont = spectral_action(spec_cont[:min_len], Lambda)

    return S_cont - S_disc, S_cont, S_disc

# === RUN TEST ===
os.makedirs("figures", exist_ok=True)
Lambda = 10.0
N_list = [4, 5]

print("=" * 60)
print("SPECTRAL DEFECT ON T^4")
print(f"pi - 3 = {PI_MINUS_3:.8f}")
print("=" * 60)

for N in N_list:
    try:
        defect, _, _ = compute_defect_torus4d(N, Lambda)
        abs_def = abs(defect)
        print(f"\nN = {N} (N^4 = {N**4})")
        print(f"  |Defect| = {abs_def:.6f}")

        best_div, best_match = 1, 0
        # Search up to d = 1000 (expect larger defect in 4D)
        for d in range(1, 1001):
            ratio = abs_def / d
            match = 100 * (1 - abs(ratio - PI_MINUS_3) / PI_MINUS_3)
            if match > best_match:
                best_match = match
                best_div = d

        print(f"  Best divisor: {best_div}")
        print(f"  |Defect|/{best_div} = {abs_def/best_div:.8f}")
        print(f"  Match to (pi-3): {best_match:.3f}%")
    except Exception as e:
        print(f"N={N}: ERROR — {e}")

print("\n=== DONE ===")