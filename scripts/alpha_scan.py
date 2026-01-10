"""
Alpha Scan v1.0
Поиск оптимального alpha, при котором defect/D максимально близок к (pi-3)
"""

import numpy as np
import matplotlib.pyplot as plt
from oscillator_sim import Membrane, PI_MINUS_3, PI

# ============================================================
# СКАНИРОВАНИЕ ALPHA
# ============================================================

def scan_alpha(N=8, steps=100, num_runs=5, alpha_range=None):
    """
    Сканируем alpha, для каждого находим лучший делитель D.
    Усредняем по num_runs запусков для стабильности.
    """
    
    if alpha_range is None:
        # От 1/200 до 1/50, включая 1/137
        alpha_range = [1/x for x in [200, 180, 160, 150, 140, 137, 135, 130, 
                                       120, 110, 100, 90, 80, 70, 60, 50]]
    
    print("=" * 70)
    print(f"ALPHA SCAN: N={N}, steps={steps}, runs={num_runs}")
    print(f"pi - 3 = {PI_MINUS_3:.8f}")
    print("=" * 70)
    print()
    print(f"{'alpha':>12}  {'1/alpha':>8}  {'defect':>10}  {'D':>6}  {'def/D':>10}  {'match%':>8}")
    print("-" * 70)
    
    results = []
    
    for alpha in sorted(alpha_range, reverse=True):
        defects = []
        
        # Несколько запусков для усреднения
        for run in range(num_runs):
            np.random.seed(run * 1000)  # Фиксированные, но разные seeds
            membrane = Membrane(N, alpha)
            
            for _ in range(steps):
                membrane.step()
            
            defect = membrane.compute_spectral_defect()
            defects.append(defect)
        
        avg_defect = np.mean(defects)
        std_defect = np.std(defects)
        
        # Ищем лучший делитель
        best_div = 1
        best_match = 0
        for div in range(1, 500):
            val = avg_defect / div
            match = 1 - abs(val - PI_MINUS_3) / PI_MINUS_3
            if match > best_match:
                best_match = match
                best_div = div
        
        inv_alpha = 1 / alpha
        def_over_d = avg_defect / best_div
        
        results.append({
            'alpha': alpha,
            'inv_alpha': inv_alpha,
            'defect': avg_defect,
            'std': std_defect,
            'D': best_div,
            'def_over_d': def_over_d,
            'match': best_match
        })
        
        marker = " <-- BEST" if best_match > 0.995 else ""
        print(f"{alpha:>12.6f}  {inv_alpha:>8.1f}  {avg_defect:>10.4f}  {best_div:>6}  "
              f"{def_over_d:>10.6f}  {best_match*100:>7.2f}%{marker}")
    
    return results


def find_optimum(results):
    """Находим alpha с максимальным match."""
    best = max(results, key=lambda x: x['match'])
    return best


def plot_results(results, filename):
    """График зависимости match от alpha."""
    
    alphas = [r['alpha'] for r in results]
    inv_alphas = [r['inv_alpha'] for r in results]
    matches = [r['match'] * 100 for r in results]
    defects = [r['defect'] for r in results]
    divs = [r['D'] for r in results]
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Match vs 1/alpha
    ax1 = axes[0, 0]
    ax1.plot(inv_alphas, matches, 'bo-', markersize=8)
    ax1.axvline(x=137, color='r', linestyle='--', label='α = 1/137')
    ax1.axhline(y=100, color='g', linestyle=':', alpha=0.5)
    ax1.set_xlabel('1/α')
    ax1.set_ylabel('Match to (π-3), %')
    ax1.set_title('Match vs α')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Defect vs 1/alpha
    ax2 = axes[0, 1]
    ax2.plot(inv_alphas, defects, 'go-', markersize=8)
    ax2.axvline(x=137, color='r', linestyle='--', label='α = 1/137')
    ax2.set_xlabel('1/α')
    ax2.set_ylabel('Spectral Defect')
    ax2.set_title('Defect vs α')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Divisor D vs 1/alpha
    ax3 = axes[1, 0]
    ax3.plot(inv_alphas, divs, 'mo-', markersize=8)
    ax3.axvline(x=137, color='r', linestyle='--', label='α = 1/137')
    ax3.set_xlabel('1/α')
    ax3.set_ylabel('Best Divisor D')
    ax3.set_title('D vs α')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Deviation from (pi-3)
    ax4 = axes[1, 1]
    deviations = [abs(r['def_over_d'] - PI_MINUS_3) for r in results]
    ax4.semilogy(inv_alphas, deviations, 'ro-', markersize=8)
    ax4.axvline(x=137, color='r', linestyle='--', label='α = 1/137')
    ax4.set_xlabel('1/α')
    ax4.set_ylabel('|defect/D - (π-3)|')
    ax4.set_title('Deviation from (π-3)')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    print(f"\nSaved: {filename}")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    
    # Основной скан
    results = scan_alpha(N=8, steps=100, num_runs=5)
    
    print()
    print("=" * 70)
    print("OPTIMUM")
    print("=" * 70)
    
    best = find_optimum(results)
    print(f"Best alpha: {best['alpha']:.6f} (1/{best['inv_alpha']:.1f})")
    print(f"Defect: {best['defect']:.4f}")
    print(f"Divisor D: {best['D']}")
    print(f"Defect/D: {best['def_over_d']:.6f}")
    print(f"Match: {best['match']*100:.2f}%")
    print(f"(π-3) = {PI_MINUS_3:.6f}")
    
    # График
    plot_results(results, 
                 r"C:\Users\admin\Documents\ontological-resolution-theory\figures\alpha_scan.png")
    
    print()
    print("=" * 70)
    print("FINE SCAN AROUND 1/137")
    print("=" * 70)
    
    # Тонкий скан вокруг 1/137
    fine_alphas = [1/x for x in range(130, 145)]
    fine_results = scan_alpha(N=8, steps=100, num_runs=10, alpha_range=fine_alphas)
    
    fine_best = find_optimum(fine_results)
    print()
    print(f"Fine scan best: alpha = 1/{fine_best['inv_alpha']:.0f}, match = {fine_best['match']*100:.2f}%")
    
    plot_results(fine_results,
                 r"C:\Users\admin\Documents\ontological-resolution-theory\figures\alpha_scan_fine.png")
    
    print()
    print("=== ALPHA SCAN COMPLETE ===")