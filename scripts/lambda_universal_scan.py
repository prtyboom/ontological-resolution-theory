"""
LAMBDA UNIVERSAL SCAN
=====================
Поиск оптимальной Lambda для каждой топологии
Проверка гипотезы: Lambda_opt ≈ N или Lambda_opt = const?

ORT Project | December 2025
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

PI_MINUS_3 = np.pi - 3

# === ОБЩИЕ ФУНКЦИИ ===

def shift_1d_antiperiodic(N):
    """1D сдвиг с антипериодичностью"""
    S = np.roll(np.eye(N), -1, axis=1)
    S[-1, 0] = -1
    return S

def spectral_action(eigenvalues, Lambda):
    """Heat kernel trace"""
    return np.sum(np.exp(-eigenvalues**2 / Lambda**2))

def find_best_divisor(defect, max_div=1000):
    """Ищет делитель дающий лучший match к pi-3"""
    best_div = 1
    best_match = 0
    for d in range(1, max_div + 1):
        ratio = abs(defect) / d
        match = 100 * (1 - abs(ratio - PI_MINUS_3) / PI_MINUS_3)
        if match > best_match:
            best_match = match
            best_div = d
    return best_div, best_match

# === S^1 (ОКРУЖНОСТЬ) ===

def spectrum_S1_continuous(N_modes):
    k = np.arange(-N_modes, N_modes + 1)
    return k + 0.5

def spectrum_S1_discrete(N):
    a = 2 * np.pi / N
    S = shift_1d_antiperiodic(N)
    nabla = (S - S.T) / (2 * a)
    zero = np.zeros((N, N))
    D = np.block([[zero, nabla], [nabla.T, zero]])
    eigs = np.sqrt(np.maximum(linalg.eigvalsh(D @ D), 0))
    return eigs

def defect_S1(N, Lambda):
    spec_d = spectrum_S1_discrete(N)
    spec_c = np.abs(spectrum_S1_continuous(N))[:len(spec_d)]
    return spectral_action(spec_c, Lambda) - spectral_action(spec_d, Lambda)

# === T^2 (2-ТОР) ===

def spectrum_T2_continuous(N_modes):
    spectrum = []
    for k in range(-N_modes, N_modes + 1):
        for m in range(-N_modes, N_modes + 1):
            lam = np.sqrt((k + 0.5)**2 + (m + 0.5)**2)
            spectrum.append(lam)
    return np.array(sorted(spectrum))

def spectrum_T2_discrete(N):
    a = 2 * np.pi / N
    S1 = shift_1d_antiperiodic(N)
    I_N = np.eye(N)
    
    Sx = np.kron(S1, I_N)
    Sy = np.kron(I_N, S1)
    
    Dx = (Sx - Sx.T) / (2 * a)
    Dy = (Sy - Sy.T) / (2 * a)
    
    N2 = N * N
    zero = np.zeros((N2, N2), dtype=complex)
    D_plus = Dx + 1j * Dy
    D_minus = Dx - 1j * Dy
    
    D = np.block([[zero, D_minus], [D_plus, zero]])
    eigs = np.sqrt(np.maximum(np.real(linalg.eigvalsh(D @ D.conj().T)), 0))
    return eigs

def defect_T2(N, Lambda):
    spec_d = spectrum_T2_discrete(N)
    spec_c = spectrum_T2_continuous(N)[:len(spec_d)]
    return spectral_action(spec_c, Lambda) - spectral_action(spec_d, Lambda)

# === T^3 (3-ТОР) ===

def spectrum_T3_continuous(N_modes):
    spectrum = []
    for k in range(-N_modes, N_modes + 1):
        for m in range(-N_modes, N_modes + 1):
            for n in range(-N_modes, N_modes + 1):
                lam = np.sqrt((k + 0.5)**2 + (m + 0.5)**2 + (n + 0.5)**2)
                spectrum.append(lam)
    return np.array(sorted(spectrum))

def spectrum_T3_discrete(N):
    a = 2 * np.pi / N
    S1 = shift_1d_antiperiodic(N)
    I_N = np.eye(N)
    
    Sx = np.kron(np.kron(S1, I_N), I_N)
    Sy = np.kron(np.kron(I_N, S1), I_N)
    Sz = np.kron(np.kron(I_N, I_N), S1)
    
    Dx = (Sx - Sx.T) / (2 * a)
    Dy = (Sy - Sy.T) / (2 * a)
    Dz = (Sz - Sz.T) / (2 * a)
    
    N3 = N ** 3
    zero = np.zeros((N3, N3), dtype=complex)
    
    # Simplified Dirac: D = [[0, Dx+iDy+Dz], [h.c., 0]]
    D_upper = Dx + 1j * Dy + Dz
    D = np.block([[zero, D_upper], [D_upper.conj().T, zero]])
    
    eigs = np.sqrt(np.maximum(np.real(linalg.eigvalsh(D @ D.conj().T)), 0))
    return eigs

def defect_T3(N, Lambda):
    spec_d = spectrum_T3_discrete(N)
    spec_c = spectrum_T3_continuous(N)[:len(spec_d)]
    return spectral_action(spec_c, Lambda) - spectral_action(spec_d, Lambda)

# === T^4 (4-ТОР) ===

def spectrum_T4_continuous(N_modes):
    spectrum = []
    for k in range(-N_modes, N_modes + 1):
        for m in range(-N_modes, N_modes + 1):
            for n in range(-N_modes, N_modes + 1):
                for p in range(-N_modes, N_modes + 1):
                    lam = np.sqrt((k + 0.5)**2 + (m + 0.5)**2 + (n + 0.5)**2 + (p + 0.5)**2)
                    spectrum.append(lam)
    return np.array(sorted(spectrum))

def spectrum_T4_discrete(N):
    a = 2 * np.pi / N
    S1 = shift_1d_antiperiodic(N)
    I_N = np.eye(N)
    
    Sx = np.kron(np.kron(np.kron(S1, I_N), I_N), I_N)
    Sy = np.kron(np.kron(np.kron(I_N, S1), I_N), I_N)
    Sz = np.kron(np.kron(np.kron(I_N, I_N), S1), I_N)
    Sw = np.kron(np.kron(np.kron(I_N, I_N), I_N), S1)
    
    Dx = (Sx - Sx.T) / (2 * a)
    Dy = (Sy - Sy.T) / (2 * a)
    Dz = (Sz - Sz.T) / (2 * a)
    Dw = (Sw - Sw.T) / (2 * a)
    
    N4 = N ** 4
    zero = np.zeros((N4, N4), dtype=complex)
    
    D_upper = Dx + 1j * Dy + Dz + 1j * Dw
    D = np.block([[zero, D_upper], [D_upper.conj().T, zero]])
    
    eigs = np.sqrt(np.maximum(np.real(linalg.eigvalsh(D @ D.conj().T)), 0))
    return eigs

def defect_T4(N, Lambda):
    spec_d = spectrum_T4_discrete(N)
    spec_c = spectrum_T4_continuous(N)[:len(spec_d)]
    return spectral_action(spec_c, Lambda) - spectral_action(spec_d, Lambda)

# === MAIN ===

print("=" * 70)
print("LAMBDA UNIVERSAL SCAN")
print(f"pi - 3 = {PI_MINUS_3:.8f}")
print("=" * 70)

# Конфигурации для тестирования
configs = [
    ("S1", 12, defect_S1),
    ("T2", 8, defect_T2),
    ("T3", 5, defect_T3),
    ("T4", 4, defect_T4),
]

Lambda_values = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 20]

results = {}

for name, N, defect_fn in configs:
    print(f"\n{'='*70}")
    print(f"{name}, N = {N}")
    print(f"{'='*70}")
    print(f"\n Lambda    |Defect|    Best Div    |Def|/Div      Match%")
    print("-" * 65)
    
    best_lambda = None
    best_match_overall = 0
    scan_results = []
    
    for L in Lambda_values:
        try:
            d = defect_fn(N, L)
            div, match = find_best_divisor(d)
            ratio = abs(d) / div
            scan_results.append((L, abs(d), div, ratio, match))
            
            marker = ""
            if match > best_match_overall:
                best_match_overall = match
                best_lambda = L
                marker = " <-- BEST"
            
            print(f"  {L:4.0f}    {abs(d):10.4f}      {div:4d}      {ratio:.6f}    {match:5.2f}%{marker}")
        except Exception as e:
            print(f"  {L:4.0f}    ERROR: {e}")
    
    results[name] = {
        'N': N,
        'best_lambda': best_lambda,
        'best_match': best_match_overall,
        'scan': scan_results
    }
    
    print(f"\n>>> {name}: Best Lambda = {best_lambda}, Match = {best_match_overall:.2f}%")

# === SUMMARY ===

print("\n" + "=" * 70)
print("SUMMARY: Optimal Lambda for each topology")
print("=" * 70)
print(f"\n{'Topology':<8} {'N':<4} {'Lambda_opt':<12} {'Match%':<10} {'Lambda/N':<10}")
print("-" * 50)

for name, data in results.items():
    N = data['N']
    L_opt = data['best_lambda']
    match = data['best_match']
    ratio = L_opt / N if L_opt else 0
    print(f"{name:<8} {N:<4} {L_opt:<12} {match:<10.2f} {ratio:<10.2f}")

# === PLOT ===

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for idx, (name, data) in enumerate(results.items()):
    ax = axes[idx // 2, idx % 2]
    scan = data['scan']
    if scan:
        Ls = [s[0] for s in scan]
        matches = [s[4] for s in scan]
        ax.plot(Ls, matches, 'bo-', lw=2, ms=8)
        ax.axhline(99, color='r', ls='--', lw=1, label='99%')
        ax.axvline(data['best_lambda'], color='g', ls=':', lw=2, 
                   label=f"Best Λ={data['best_lambda']}")
        ax.axvline(data['N'], color='orange', ls=':', lw=2, alpha=0.7,
                   label=f"N={data['N']}")
        ax.set_xlabel('Lambda', fontsize=12)
        ax.set_ylabel('Match %', fontsize=12)
        ax.set_title(f"{name}, N={data['N']}", fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_ylim([min(matches) - 5, 102])

plt.tight_layout()
plt.savefig(r'C:\Users\admin\Documents\ontological-resolution-theory\figures\lambda_universal_scan.png', dpi=150)
print("\nSaved: lambda_universal_scan.png")
plt.show()

print("\n" + "=" * 70)
print("SCAN COMPLETE")
print("=" * 70)