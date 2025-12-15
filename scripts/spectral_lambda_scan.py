"""
SPECTRAL LAMBDA SCAN: Устойчивость при разных UV-cutoff
Проверяем: сохраняется ли 52*(pi-3) при разных Lambda?
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
print("LAMBDA SCAN: Устойчивость отношения |Defect|/52 к (pi-3)")
print("=" * 70)

N = 12  # фиксируем зодиакальную конфигурацию

Lambda_list = [1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 100]

print(f"\nN = {N} (фиксировано)")
print(f"Target: pi - 3 = {PI_MINUS_3:.8f}")
print()
print("  Lambda    |Defect|    |Def|/52    pi-3       Match%")
print("-" * 60)

results = []
for L in Lambda_list:
    d = abs(compute(N, L))
    ratio_52 = d / 52
    match_pct = 100 * (1 - abs(ratio_52 - PI_MINUS_3) / PI_MINUS_3)
    results.append((L, d, ratio_52, match_pct))
    print(f"{L:7.1f}  {d:10.6f}  {ratio_52:10.8f}  {PI_MINUS_3:.8f}  {match_pct:6.2f}%")

print("\n" + "=" * 70)
print("АНАЛИЗ")
print("=" * 70)

# Лучший Lambda
best = max(results, key=lambda x: x[3])
print(f"\nЛучшее совпадение при Lambda = {best[0]}")
print(f"  |Defect|/52 = {best[2]:.8f}")
print(f"  Match = {best[3]:.2f}%")

# Среднее
avg_match = np.mean([r[3] for r in results])
print(f"\nСреднее совпадение по всем Lambda: {avg_match:.2f}%")

# Стабильность
std_ratio = np.std([r[2] for r in results])
mean_ratio = np.mean([r[2] for r in results])
print(f"Среднее |Defect|/52: {mean_ratio:.8f}")
print(f"Стандартное отклонение: {std_ratio:.8f}")
print(f"Относительная стабильность: {100*std_ratio/mean_ratio:.2f}%")

# Проверка других делителей
print("\n" + "=" * 70)
print("ПРОВЕРКА ДРУГИХ ДЕЛИТЕЛЕЙ при Lambda=10")
print("=" * 70)

d_10 = abs(compute(12, 10))
for divisor in [48, 50, 51, 52, 53, 54, 56]:
    ratio = d_10 / divisor
    match = 100 * (1 - abs(ratio - PI_MINUS_3) / PI_MINUS_3)
    marker = " <-- BEST" if divisor == 52 else ""
    print(f"  |Defect|/{divisor} = {ratio:.8f}  Match: {match:.2f}%{marker}")

# Графики
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

Lambdas = [r[0] for r in results]
defects = [r[1] for r in results]
ratios = [r[2] for r in results]
matches = [r[3] for r in results]

# 1. |Defect| vs Lambda
ax1 = axes[0, 0]
ax1.semilogx(Lambdas, defects, 'bo-', lw=2, ms=8)
ax1.set_xlabel('Lambda (UV cutoff)', fontsize=12)
ax1.set_ylabel('|Defect|', fontsize=12)
ax1.set_title('Спектральный дефект vs Lambda (N=12)', fontsize=14)
ax1.grid(True, alpha=0.3)

# 2. |Defect|/52 vs Lambda
ax2 = axes[0, 1]
ax2.semilogx(Lambdas, ratios, 'go-', lw=2, ms=8)
ax2.axhline(PI_MINUS_3, color='r', ls='--', lw=2, label=f'pi-3 = {PI_MINUS_3:.4f}')
ax2.set_xlabel('Lambda', fontsize=12)
ax2.set_ylabel('|Defect| / 52', fontsize=12)
ax2.set_title('Нормированный дефект vs Lambda', fontsize=14)
ax2.legend()
ax2.grid(True, alpha=0.3)

# 3. Match % vs Lambda
ax3 = axes[1, 0]
ax3.semilogx(Lambdas, matches, 'mo-', lw=2, ms=8)
ax3.axhline(100, color='r', ls='--', lw=2, label='100%')
ax3.axhline(99, color='orange', ls=':', lw=2, label='99%')
ax3.set_xlabel('Lambda', fontsize=12)
ax3.set_ylabel('Match %', fontsize=12)
ax3.set_title('Точность совпадения с (pi-3)', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3)
ax3.set_ylim([min(matches)-5, 105])

# 4. Сравнение делителей
ax4 = axes[1, 1]
divisors = [48, 49, 50, 51, 52, 53, 54, 55, 56]
matches_div = []
for div in divisors:
    ratio = d_10 / div
    match = 100 * (1 - abs(ratio - PI_MINUS_3) / PI_MINUS_3)
    matches_div.append(match)

ax4.bar(divisors, matches_div, color='cyan', alpha=0.7)
ax4.axvline(52, color='red', ls='-', lw=3, label='52')
ax4.set_xlabel('Divisor', fontsize=12)
ax4.set_ylabel('Match %', fontsize=12)
ax4.set_title('Какой делитель даёт лучшее совпадение?', fontsize=14)
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(r'C:\Users\admin\Documents\ontological-resolution-theory\figures\spectral_lambda_scan.png', 
            dpi=150)
print("\nSaved: spectral_lambda_scan.png")
plt.show()

print("\n" + "=" * 70)
print("LAMBDA SCAN COMPLETE")
print("=" * 70)