"""
N Systematic Scan v1.1
Систематический анализ: как α_opt и D_opt зависят от N
+ бэйзлайн для случайных "дефектов".
"""

import numpy as np
import matplotlib.pyplot as plt
from oscillator_sim import Membrane, PI_MINUS_3, PI


# ============================================================
# СКАНИРОВАНИЕ ПО α ДЛЯ ДАННОГО N
# ============================================================

def scan_for_N(N, steps=200, num_runs=10):
    """
    Для данного N сканируем α и находим оптимум по совпадению
    avg_defect / D ≈ (pi - 3), подбирая целый D от 1 до max_div.
    """
    # Диапазон α: от 1/300 до 1/30 (можно менять при необходимости)
    alpha_values = [1/x for x in [
        300, 250, 200, 180, 160, 150, 140, 137,
        135, 130, 120, 110, 100, 80, 60, 50, 40, 30
    ]]

    best_match = 0.0
    best_alpha = None
    best_D = None
    best_defect = None
    results = []

    for alpha in alpha_values:
        defects = []

        for run in range(num_runs):
            # Детерминированные seeds для воспроизводимости
            np.random.seed(run * 1000 + N)
            membrane = Membrane(N, alpha)

            for _ in range(steps):
                membrane.step()

            defect = membrane.compute_spectral_defect()
            defects.append(defect)

        avg_defect = np.mean(defects)

        # Ищем лучший делитель D, приближающий avg_defect/D к (pi - 3)
        best_div_local = 1
        best_match_local = 0.0
        for div in range(1, 1000):
            val = avg_defect / div
            match = 1.0 - abs(val - PI_MINUS_3) / PI_MINUS_3
            if match > best_match_local:
                best_match_local = match
                best_div_local = div

        results.append({
            "alpha": alpha,
            "inv_alpha": 1.0 / alpha,
            "defect": avg_defect,
            "D": best_div_local,
            "match": best_match_local,
        })

        if best_match_local > best_match:
            best_match = best_match_local
            best_alpha = alpha
            best_D = best_div_local
            best_defect = avg_defect

    return {
        "N": N,
        "alpha_opt": best_alpha,
        "inv_alpha_opt": 1.0 / best_alpha if best_alpha else None,
        "D_opt": best_D,
        "defect_opt": best_defect,
        "match_opt": best_match,
        "all_results": results,
    }


def run_systematic_scan(N_values, steps=200, num_runs=10):
    """
    Запуск систематического скана по всем N.
    """
    print("=" * 80)
    print(f"SYSTEMATIC N SCAN: steps={steps}, runs={num_runs}")
    print(f"pi - 3 = {PI_MINUS_3:.8f}")
    print("=" * 80)
    print()

    all_data = []

    for N in N_values:
        print(f"\n{'=' * 60}")
        print(f"N = {N} (N² = {N * N} oscillators)")
        print("=" * 60)

        result = scan_for_N(N, steps, num_runs)
        all_data.append(result)

        print(f"\n  Optimal α = 1/{result['inv_alpha_opt']:.1f}")
        print(f"  Optimal D = {result['D_opt']}")
        print(f"  Defect = {result['defect_opt']:.4f}")
        print(f"  Defect/D = {result['defect_opt'] / result['D_opt']:.6f}")
        print(f"  Match = {result['match_opt'] * 100:.2f}%")

        # Таблица для этого N
        print(f"\n  {'1/α':>8}  {'defect':>10}  {'D':>6}  {'match%':>8}")
        print(f"  {'-' * 40}")
        for r in result["all_results"]:
            marker = " <--" if r["match"] == result["match_opt"] else ""
            print(
                f"  {r['inv_alpha']:>8.1f}  {r['defect']:>10.4f}  {r['D']:>6}  "
                f"{r['match'] * 100:>7.2f}%{marker}"
            )

    return all_data


# ============================================================
# АНАЛИЗ ТРЕНДОВ ПО N
# ============================================================

def analyze_trends(all_data):
    """
    Анализ трендов: как α_opt и D_opt зависят от N.
    """
    print("\n" + "=" * 80)
    print("TREND ANALYSIS")
    print("=" * 80)

    Ns = [d["N"] for d in all_data]
    N2s = [d["N"] ** 2 for d in all_data]
    inv_alphas = [d["inv_alpha_opt"] for d in all_data]
    Ds = [d["D_opt"] for d in all_data]
    matches = [d["match_opt"] for d in all_data]

    print(f"\n{'N':>6}  {'N²':>8}  {'1/α_opt':>10}  {'D_opt':>8}  {'D/N²':>10}  {'match%':>8}")
    print("-" * 70)

    for d in all_data:
        D_over_N2 = d["D_opt"] / (d["N"] ** 2)
        print(
            f"{d['N']:>6}  {d['N'] ** 2:>8}  {d['inv_alpha_opt']:>10.1f}  "
            f"{d['D_opt']:>8}  {D_over_N2:>10.4f}  {d['match_opt'] * 100:>7.2f}%"
        )

    # Проверяем гипотезы
    print("\n" + "-" * 70)
    print("HYPOTHESIS TESTS:")

    # H1: D ~ N²
    D_over_N2_values = [d["D_opt"] / (d["N"] ** 2) for d in all_data]
    avg_ratio = np.mean(D_over_N2_values)
    std_ratio = np.std(D_over_N2_values)
    print(f"\n  H1: D ∝ N²")
    print(f"      D/N² = {avg_ratio:.4f} ± {std_ratio:.4f}")

    # H2: D ~ N
    D_over_N_values = [d["D_opt"] / d["N"] for d in all_data]
    avg_ratio_N = np.mean(D_over_N_values)
    std_ratio_N = np.std(D_over_N_values)
    print(f"\n  H2: D ∝ N")
    print(f"      D/N = {avg_ratio_N:.4f} ± {std_ratio_N:.4f}")

    # H3: α_opt ~ 1/N
    alpha_times_N = [d["N"] / d["inv_alpha_opt"] for d in all_data]
    avg_aN = np.mean(alpha_times_N)
    std_aN = np.std(alpha_times_N)
    print(f"\n  H3: α_opt ∝ 1/N")
    print(f"      α·N = {avg_aN:.4f} ± {std_aN:.4f}")

    # H4: α_opt → 1/137 при больших N?
    print(f"\n  H4: α_opt → 1/137 at large N?")
    print(f"      Current trend: 1/α_opt = {inv_alphas}")

    # H5: D_opt → 43 для какого-то N?
    print(f"\n  H5: D_opt → 43 at some N?")
    for d in all_data:
        if abs(d["D_opt"] - 43) < 5:
            print(f"      N={d['N']}: D_opt = {d['D_opt']} (close to 43!)")

    return {
        "Ns": Ns,
        "N2s": N2s,
        "inv_alphas": inv_alphas,
        "Ds": Ds,
        "matches": matches,
    }


# ============================================================
# БЭЙЗЛАЙН: СЛУЧАЙНЫЕ DEFECT
# ============================================================

def random_defect_baseline(num_samples=2000, defect_range=(0.5, 8.0), max_div=1000):
    """
    Бэйзлайн: берём случайные значения defect из равномерного диапазона
    и смотрим, какие match% даёт тот же алгоритм подбора D.
    Это показывает, насколько "удивителен" высокий match сам по себе.
    """
    matches = []

    for i in range(num_samples):
        defect = np.random.uniform(*defect_range)

        best_match_local = 0.0
        best_div_local = 1

        for div in range(1, max_div):
            val = defect / div
            match = 1.0 - abs(val - PI_MINUS_3) / PI_MINUS_3
            if match > best_match_local:
                best_match_local = match
                best_div_local = div

        matches.append(best_match_local)

    matches = np.array(matches)

    print("\n" + "=" * 80)
    print("RANDOM DEFECT BASELINE")
    print("=" * 80)
    print(f"Samples: {num_samples}, defect_range={defect_range}, max_div={max_div}")
    print(f"match (mean ± std): {matches.mean() * 100:.2f}% ± {matches.std() * 100:.2f}%")
    print("match 5/50/95 percentiles: "
          f"{np.percentile(matches, 5) * 100:.2f}%, "
          f"{np.percentile(matches, 50) * 100:.2f}%, "
          f"{np.percentile(matches, 95) * 100:.2f}%")

    return matches


# ============================================================
# ГРАФИКИ ТРЕНДОВ
# ============================================================

def plot_trends(all_data, filename):
    """
    Графики трендов.
    """
    Ns = [d["N"] for d in all_data]
    N2s = [d["N"] ** 2 for d in all_data]
    inv_alphas = [d["inv_alpha_opt"] for d in all_data]
    Ds = [d["D_opt"] for d in all_data]
    matches = [d["match_opt"] * 100 for d in all_data]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. D_opt vs N
    ax1 = axes[0, 0]
    ax1.plot(Ns, Ds, "bo-", markersize=10, linewidth=2)
    ax1.axhline(y=43, color="r", linestyle="--", label="D = 43 (Q_crit)")
    ax1.set_xlabel("N")
    ax1.set_ylabel("D_opt")
    ax1.set_title("Optimal Divisor vs Grid Size")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 2. D_opt vs N²
    ax2 = axes[0, 1]
    ax2.plot(N2s, Ds, "go-", markersize=10, linewidth=2)
    ax2.axhline(y=43, color="r", linestyle="--", label="D = 43")
    if len(N2s) > 2:
        z = np.polyfit(N2s, Ds, 1)
        p = np.poly1d(z)
        ax2.plot(N2s, p(N2s), "g--", alpha=0.5,
                 label=f"Linear fit: D ≈ {z[0]:.4f}·N² + {z[1]:.1f}")
    ax2.set_xlabel("N²")
    ax2.set_ylabel("D_opt")
    ax2.set_title("D_opt vs N² (checking D ∝ N²)")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # 3. 1/α_opt vs N
    ax3 = axes[1, 0]
    ax3.plot(Ns, inv_alphas, "mo-", markersize=10, linewidth=2)
    ax3.axhline(y=137, color="r", linestyle="--", label="1/α = 137")
    ax3.set_xlabel("N")
    ax3.set_ylabel("1/α_opt")
    ax3.set_title("Optimal 1/α vs Grid Size")
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # 4. Match vs N
    ax4 = axes[1, 1]
    ax4.plot(Ns, matches, "ro-", markersize=10, linewidth=2)
    ax4.axhline(y=100, color="g", linestyle=":", alpha=0.5)
    ax4.set_xlabel("N")
    ax4.set_ylabel("Match to (π-3), %")
    ax4.set_title("Best Match vs Grid Size")
    ax4.set_ylim([95, 101])
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    print(f"\nSaved: {filename}")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    # Значения N для тестирования
    N_values = [4, 6, 8, 10, 12, 16]

    print("Starting systematic scan...")
    print(f"N values: {N_values}")
    print()

    all_data = run_systematic_scan(N_values, steps=200, num_runs=10)

    trends = analyze_trends(all_data)

    # Бэйзлайн по случайным "дефектам"
    random_defect_baseline(
        num_samples=2000,
        defect_range=(0.5, 8.0),
        max_div=1000
    )

    plot_trends(
        all_data,
        r"C:\Users\admin\Documents\ontological-resolution-theory\figures\N_systematic_scan.png"
    )

    print("\n" + "=" * 80)
    print("SYSTEMATIC SCAN COMPLETE")
    print("=" * 80)