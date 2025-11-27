#!/usr/bin/env python3
"""
09_dwarf_analysis.py - Ultra-faint dwarf galaxy M/L analysis
ORT prediction: high M/L independent of star formation history
"""

import numpy as np
import matplotlib.pyplot as plt
import os

print("="*60)
print("DWARF GALAXY ANALYSIS - ORT PREDICTION TEST")
print("="*60)

# Load data
data = []
with open('../data/dwarf_galaxies.csv', 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split(',')
        data.append({
            'name': parts[0],
            'L': float(parts[1]),
            'sigma': float(parts[2]),
            'r_half': float(parts[3]),
            'ML_obs': float(parts[4])
        })

print(f"Loaded {len(data)} dwarf galaxies\n")

# Calculate dynamical mass: M_dyn = 4 * sigma^2 * r_half / G
G = 4.302e-3  # pc (km/s)^2 / Msun

results = []
for d in data:
    M_dyn = 4 * d['sigma']**2 * d['r_half'] / G
    ML_calc = M_dyn / d['L']
    results.append({
        'name': d['name'],
        'L': d['L'],
        'sigma': d['sigma'],
        'r_half': d['r_half'],
        'M_dyn': M_dyn,
        'ML_calc': ML_calc,
        'ML_obs': d['ML_obs']
    })

# Print results
print(f"{'Galaxy':<15} {'L (Lsun)':<12} {'σ (km/s)':<10} {'M/L_calc':<10} {'M/L_obs':<10}")
print("-"*60)
for r in results:
    print(f"{r['name']:<15} {r['L']:<12.0f} {r['sigma']:<10.1f} {r['ML_calc']:<10.0f} {r['ML_obs']:<10.0f}")

# ORT Prediction test
print(f"\n{'='*60}")
print("ORT PREDICTION TEST")
print("="*60)
print("Prediction: High M/L (>>1) even for ultra-faint dwarfs")
print("           because DM is primordial, not from local processes\n")

ultra_faint = [r for r in results if r['L'] < 10000]
classical = [r for r in results if r['L'] >= 100000]

ML_uf = np.mean([r['ML_calc'] for r in ultra_faint])
ML_cl = np.mean([r['ML_calc'] for r in classical])

print(f"Ultra-faint dwarfs (L < 10^4 Lsun): <M/L> = {ML_uf:.0f}")
print(f"Classical dwarfs (L > 10^5 Lsun):   <M/L> = {ML_cl:.0f}")
print(f"\nRatio: {ML_uf/ML_cl:.1f}x")

if ML_uf > 100:
    print("\n✓ CONFIRMED: Ultra-faint dwarfs are DM-dominated")
    print("  This supports primordial DM origin (ORT hypothesis)")
else:
    print("\n✗ NOT CONFIRMED")

# Plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left: M/L vs Luminosity
ax = axes[0]
L_arr = np.array([r['L'] for r in results])
ML_arr = np.array([r['ML_calc'] for r in results])
names = [r['name'] for r in results]

ax.scatter(L_arr, ML_arr, s=100, c='blue', edgecolor='k', zorder=5)
for i, name in enumerate(names):
    ax.annotate(name, (L_arr[i], ML_arr[i]), fontsize=8, 
                xytext=(5,5), textcoords='offset points')

ax.axhline(1, color='green', ls='--', label='Stellar only (M/L=1)')
ax.axhline(5.4, color='red', ls=':', lw=2, label='Cosmic DM/baryon=5.4')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Luminosity (L☉)', fontsize=12)
ax.set_ylabel('M/L (M☉/L☉)', fontsize=12)
ax.set_title('Dwarf Galaxy Mass-to-Light Ratios', fontsize=14)
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xlim(100, 1e8)
ax.set_ylim(1, 1e4)

# Right: DM fraction
ax = axes[1]
f_DM = 1 - 1/ML_arr  # approximate DM fraction
ax.bar(range(len(names)), f_DM*100, color='darkblue', edgecolor='k')
ax.set_xticks(range(len(names)))
ax.set_xticklabels(names, rotation=45, ha='right', fontsize=9)
ax.set_ylabel('Dark Matter Fraction (%)', fontsize=12)
ax.set_title('DM Dominance in Dwarf Galaxies', fontsize=14)
ax.axhline(84, color='red', ls=':', lw=2, label='Cosmic: 84%')
ax.legend()
ax.set_ylim(0, 100)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
os.makedirs('../results/figures', exist_ok=True)
plt.savefig('../results/figures/dwarf_galaxies.png', dpi=150)
print(f"\nSaved: ../results/figures/dwarf_galaxies.png")
plt.show()