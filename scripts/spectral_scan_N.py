"""
SPECTRAL SCAN: Поиск аномалий при разных N
Проверяем простые числа и около-зодиакальные
ORT Project | June 2025
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

PI_MINUS_3 = np.pi - 3

def dirac_discrete(N):
    a = 2 * np.pi / N
    shift = np.roll(np.eye(N), -1, axis=1)
    shift[-1, 0] = -1
    nabla = (shift - shift.T) / (2 * a)
    zero = np.zeros((N, N))
    D = np.block([[zero, nabla], [nabla.T, zero]])
    return D

def dirac_continuous(N_modes):
    k = np.arange(-N_modes, N_modes + 1)
    return np.sort(k + 0.5)

def spectral_action(spectrum, Lambda):
    return np.sum(np.exp(-spectrum**2 / Lambda**2))

def compute(N, Lambda):
    D = dirac_discrete(N)
    eigs = np.sqrt(np.maximum(linalg.eigvalsh(D @ D), 0))
    spec_c = np.abs(dirac_continuous(N)[:2*N])
    S_c = spectral_action(spec_c, Lambda)
    S_d = spectral_action(eigs, Lambda)
    return S_c - S_d

print("=" * 70)
print("SPECTRAL SCAN: Поиск аномалий")
print("Target: (pi - 3) =", PI_MINUS_3)
print("=" * 70)

Lambda = 10.0

# Расширенный список N
N_list = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 24, 30, 36]

print("\n  N    |Defect|   |Def|/(pi-3)   |Def|/N    Ratio/N   Note")
print("-" * 70)

results = []
for N in N_list:
    d = compute(N, Lambda)
    ratio = abs(d) / PI_MINUS_3
    per_N = abs(d) / N
    ratio_N = ratio / N
    
    # Помечаем особые числа
    note = ""
    if N == 12:
        note = "<-- ZODIAC"
    elif N in [7, 11, 13, 17]:
        note = "<-- PRIME"
    elif N in [4, 6, 8]:
        note = "<-- 2^n or 2*3"
    
    # Проверка близости к целому
    nearest_int = round(ratio)
    if abs(ratio - nearest_int) < 0.1:
        note += f" ~{nearest_int}"
    
    results.append((N, abs(d), ratio, per_N, ratio_N))
    print(f"{N:3d}  {abs(d):9.4f}  {ratio:11.4f}  {per_N:9.4f}  {ratio_N:9.4f}   {note}")

print("\n" + "=" * 70)
print("АНАЛИЗ ПАТТЕРНОВ")
print("=" * 70)

# Ищем линейную зависимость ratio от N
Ns = np.array([r[0] for r in results])
ratios = np.array([r[2] for r in results])

# Линейный fit: ratio = a*N + b
coeffs = np.polyfit(Ns, ratios, 1)
print(f"\nЛинейный fit: |Defect|/(pi-3) = {coeffs[0]:.4f}*N + {coeffs[1]:.4f}")

# Квадратичный fit
coeffs2 = np.polyfit(Ns, ratios, 2)
print(f"Квадратичный: |Defect|/(pi-3) = {coeffs2[0]:.4f}*N^2 + {coeffs2[1]:.4f}*N + {coeffs2[2]:.4f}")

# Специальный анализ N=12
print("\n" + "=" * 70)
print("ФОКУС: Сравнение N=11, 12, 13")
print("=" * 70)

for N in [11, 12, 13]:
    d = compute(N, Lambda)
    ratio = abs(d) / PI_MINUS_3
    divisor = ratio / N
    print(f"\nN={N}:")
    print(f"  |Defect| = {abs(d):.6f}")
    print(f"  |Defect|/(pi-3) = {ratio:.4f}")
    print(f"  |Defect|/(N*(pi-3)) = {divisor:.4f}")
    print(f"  |Defect|/{N*4} = {abs(d)/(N*4):.6f}")
    if abs(abs(d)/(N*4) - PI_MINUS_3) / PI_MINUS_3 < 0.1:
        print(f"  --> MATCH: |Defect| ~ 4*N*(pi-3)")

# График
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 1. |Defect|/(pi-3) vs N
ax1 = axes[0, 0]
ax1.plot(Ns, ratios, 'bo-', lw=2, ms=8)
ax1.plot(Ns, np.polyval(coeffs, Ns), 'r--', lw=2, label=f'Linear fit')
ax1.axvline(12, color='gold', ls=':', lw=2, label='N=12')
for N in [7, 11, 13, 17]:
    if N in Ns:
        ax1.axvline(N, color='green', ls=':', lw=1, alpha=0.5)
ax1.set_xlabel('N', fontsize=12)
ax1.set_ylabel('|Defect| / (pi-3)', fontsize=12)
ax1.set_title('Отношение к (pi-3) vs N', fontsize=14)
ax1.legend()
ax1.grid(True, alpha=0.3)

# 2. Отклонение от линейного fit
ax2 = axes[0, 1]
residuals = ratios - np.polyval(coeffs, Ns)
ax2.bar(Ns, residuals, color='purple', alpha=0.7)
ax2.axhline(0, color='black', lw=1)
ax2.axvline(12, color='gold', ls=':', lw=2)
ax2.set_xlabel('N', fontsize=12)
ax2.set_ylabel('Residual', fontsize=12)
ax2.set_title('Отклонение от линейного закона', fontsize=14)
ax2.grid(True, alpha=0.3)

# 3. |Defect|/(4*N) vs N — проверка гипотезы
ax3 = axes[1, 0]
normalized = [abs(r[1])/(4*r[0]) for r in results]
ax3.plot(Ns, normalized, 'go-', lw=2, ms=8)
ax3.axhline(PI_MINUS_3, color='r', ls='--', lw=2, label=f'pi-3 = {PI_MINUS_3:.4f}')
ax3.axvline(12, color='gold', ls=':', lw=2)
ax3.set_xlabel('N', fontsize=12)
ax3.set_ylabel('|Defect| / (4*N)', fontsize=12)
ax3.set_title('Гипотеза: |Defect| = 4*N*(pi-3)?', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3)

# 4. Ratio/N vs N — ищем константу
ax4 = axes[1, 1]
ratio_per_N = [r[2]/r[0] for r in results]
ax4.plot(Ns, ratio_per_N, 'mo-', lw=2, ms=8)
ax4.axhline(4, color='r', ls='--', lw=2, label='4')
ax4.axhline(np.mean(ratio_per_N), color='blue', ls=':', lw=2, 
            label=f'Mean = {np.mean(ratio_per_N):.2f}')
ax4.axvline(12, color='gold', ls=':', lw=2)
ax4.set_xlabel('N', fontsize=12)
ax4.set_ylabel('|Defect| / (N*(pi-3))', fontsize=12)
ax4.set_title('Нормировка на N: ищем константу', fontsize=14)
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(r'C:\Users\admin\Documents\ontological-resolution-theory\figures\spectral_scan_N.png', 
            dpi=150)
print("\nSaved: spectral_scan_N.png")
plt.show()

print("\n" + "=" * 70)
print("SCAN COMPLETE")
print("=" * 70)