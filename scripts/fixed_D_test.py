"""
Fixed D Test v1.0
Фиксируем D = N², смотрим при каком α ошибка |defect/D - (π-3)| минимальна.
Сравниваем с baseline (случайные дефекты).
"""

import numpy as np
import matplotlib.pyplot as plt
from oscillator_sim import Membrane, PI_MINUS_3, PI

# ============================================================
# ТЕСТ С ФИКСИРОВАННЫМ D
# ============================================================

def test_fixed_D(N, D, steps=200, num_runs=10):
    """
    Для данного N и фиксированного D сканируем α.
    Возвращаем ошибку |defect/D - (π-3)| для каждого α.
    """
    
    alpha_values = [1/x for x in [300, 250, 200, 180, 160, 150, 145, 140, 
                                   137, 135, 130, 125, 120, 110, 100, 80, 60]]
    
    results = []
    
    for alpha in alpha_values:
        defects = []
        
        for run in range(num_runs):
            np.random.seed(run * 1000 + N)
            membrane = Membrane(N, alpha)
            
            for _ in range(steps):
                membrane.step()
            
            defect = membrane.compute_spectral_defect()
            defects.append(defect)
        
        avg_defect = np.mean(defects)
        std_defect = np.std(defects)
        
        ratio = avg_defect / D
        error = abs(ratio - PI_MINUS_3)
        rel_error = error / PI_MINUS_3
        
        results.append({
            'alpha': alpha,
            'inv_alpha': 1/alpha,
            'defect': avg_defect,
            'std': std_defect,
            'ratio': ratio,
            'error': error,
            'rel_error': rel_error
        })
    
    return results


def baseline_fixed_D(D, num_samples=2000, defect_range=(0.5, 8.0)):
    """
    Baseline: случайные дефекты с тем же фиксированным D.
    """
    
    errors = []
    
    for _ in range(num_samples):
        defect = np.random.uniform(*defect_range)
        ratio = defect / D
        error = abs(ratio - PI_MINUS_3)
        errors.append(error)
    
    return {
        'mean': np.mean(errors),
        'std': np.std(errors),
        'median': np.median(errors),
        'p5': np.percentile(errors, 5),
        'p95': np.percentile(errors, 95)
    }


def run_test(N_values, steps=200, num_runs=10):
    """
    Запуск теста для всех N.
    """
    
    print("=" * 80)
    print("FIXED D TEST: D = N²")
    print(f"π - 3 = {PI_MINUS_3:.8f}")
    print("=" * 80)
    
    all_results = {}
    
    for N in N_values:
        D = N * N
        
        print(f"\n{'='*60}")
        print(f"N = {N}, D = N² = {D}")
        print("=" * 60)
        
        # Тест динамики
        results = test_fixed_D(N, D, steps, num_runs)
        
        # Baseline
        # Масштабируем диапазон дефектов под типичные значения
        typical_defects = [r['defect'] for r in results]
        defect_min = min(typical_defects) * 0.5
        defect_max = max(typical_defects) * 1.5
        baseline = baseline_fixed_D(D, defect_range=(defect_min, defect_max))
        
        # Находим лучший α
        best = min(results, key=lambda x: x['error'])
        
        print(f"\nBaseline (random defects, same D={D}):")
        print(f"  Error mean ± std: {baseline['mean']:.6f} ± {baseline['std']:.6f}")
        print(f"  Error median: {baseline['median']:.6f}")
        print(f"  Error 5-95%: [{baseline['p5']:.6f}, {baseline['p95']:.6f}]")
        
        print(f"\nDynamics scan:")
        print(f"  {'1/α':>8}  {'defect':>10}  {'def/D':>10}  {'error':>10}  {'vs baseline':>12}")
        print(f"  {'-'*55}")
        
        for r in results:
            # Сравнение с baseline
            if r['error'] < baseline['p5']:
                vs_base = "BETTER"
            elif r['error'] < baseline['median']:
                vs_base = "good"
            elif r['error'] > baseline['p95']:
                vs_base = "WORSE"
            else:
                vs_base = "~"
            
            marker = " <-- BEST" if r == best else ""
            print(f"  {r['inv_alpha']:>8.1f}  {r['defect']:>10.4f}  {r['ratio']:>10.6f}  "
                  f"{r['error']:>10.6f}  {vs_base:>12}{marker}")
        
        print(f"\nBest α = 1/{best['inv_alpha']:.0f}")
        print(f"  defect = {best['defect']:.4f}")
        print(f"  defect/D = {best['ratio']:.6f}")
        print(f"  error = {best['error']:.6f}")
        print(f"  (π-3) = {PI_MINUS_3:.6f}")
        
        # Сигнал?
        if best['error'] < baseline['p5']:
            print(f"  >>> SIGNAL: error < baseline 5th percentile!")
        elif best['error'] < baseline['median']:
            print(f"  >>> Weak signal: error < baseline median")
        else:
            print(f"  >>> No signal: error within baseline noise")
        
        all_results[N] = {
            'D': D,
            'results': results,
            'baseline': baseline,
            'best': best
        }
    
    return all_results


def plot_results(all_results, filename):
    """
    График: error vs α для каждого N, с baseline.
    """
    
    n_plots = len(all_results)
    fig, axes = plt.subplots(2, (n_plots + 1) // 2, figsize=(14, 10))
    axes = axes.flatten()
    
    for idx, (N, data) in enumerate(all_results.items()):
        ax = axes[idx]
        
        inv_alphas = [r['inv_alpha'] for r in data['results']]
        errors = [r['error'] for r in data['results']]
        
        ax.semilogy(inv_alphas, errors, 'bo-', markersize=6, label='Dynamics')
        
        # Baseline
        ax.axhline(y=data['baseline']['median'], color='r', linestyle='-', 
                   label=f"Baseline median", alpha=0.7)
        ax.axhline(y=data['baseline']['p5'], color='r', linestyle='--', 
                   label=f"Baseline 5%", alpha=0.5)
        ax.axhline(y=data['baseline']['p95'], color='r', linestyle='--', 
                   alpha=0.5)
        
        # 1/137
        ax.axvline(x=137, color='g', linestyle=':', alpha=0.7, label='α=1/137')
        
        ax.set_xlabel('1/α')
        ax.set_ylabel('|defect/D - (π-3)|')
        ax.set_title(f'N = {N}, D = {data["D"]}')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    
    # Убираем лишние subplot'ы
    for idx in range(len(all_results), len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    print(f"\nSaved: {filename}")


def summary(all_results):
    """
    Итоговая сводка.
    """
    
    print("\n" + "=" * 80)
    print("SUMMARY: Fixed D = N² Test")
    print("=" * 80)
    
    print(f"\n{'N':>4}  {'D=N²':>6}  {'1/α_opt':>8}  {'error':>10}  {'baseline_med':>12}  {'signal?':>10}")
    print("-" * 70)
    
    signals = []
    
    for N, data in all_results.items():
        best = data['best']
        base_med = data['baseline']['median']
        base_p5 = data['baseline']['p5']
        
        if best['error'] < base_p5:
            signal = "STRONG"
            signals.append((N, 'strong'))
        elif best['error'] < base_med:
            signal = "weak"
            signals.append((N, 'weak'))
        else:
            signal = "-"
        
        print(f"{N:>4}  {data['D']:>6}  {best['inv_alpha']:>8.0f}  {best['error']:>10.6f}  "
              f"{base_med:>12.6f}  {signal:>10}")
    
    print("\n" + "-" * 70)
    
    if signals:
        strong = [s for s in signals if s[1] == 'strong']
        weak = [s for s in signals if s[1] == 'weak']
        
        if strong:
            print(f"Strong signal at N = {[s[0] for s in strong]}")
        if weak:
            print(f"Weak signal at N = {[s[0] for s in weak]}")
    else:
        print("No signal detected: dynamics indistinguishable from random baseline")
    
    # Проверяем, есть ли α = 1/137 среди лучших
    alpha_137_count = sum(1 for N, data in all_results.items() 
                          if abs(data['best']['inv_alpha'] - 137) < 5)
    print(f"\nα_opt ≈ 1/137 in {alpha_137_count} out of {len(all_results)} cases")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    
    N_values = [4, 6, 8, 10, 12]
    
    print("Starting Fixed D Test...")
    print(f"N values: {N_values}")
    print(f"D = N² (fixed, no fitting)")
    print()
    
    all_results = run_test(N_values, steps=200, num_runs=10)
    
    summary(all_results)
    
    plot_results(all_results, 
                 r"C:\Users\admin\Documents\ontological-resolution-theory\figures\fixed_D_test.png")
    
    print("\n" + "=" * 80)
    print("FIXED D TEST COMPLETE")
    print("=" * 80)