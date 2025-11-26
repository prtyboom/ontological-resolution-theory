#!/usr/bin/env python3
"""
03_dm_baryon_ratio.py - Dark matter to baryon ratio calculation
"""

import numpy as np
import matplotlib.pyplot as plt
import json
import os

# Planck 2018 observations
OMEGA_B_H2 = 0.02237
OMEGA_CDM_H2 = 0.1200
OMEGA_RATIO = OMEGA_CDM_H2 / OMEGA_B_H2  # = 5.36

G_SM = 106.75

print("="*60)
print("DARK MATTER TO BARYON RATIO")
print("="*60)
print(f"\nObserved (Planck 2018):")
print(f"  Omega_b h^2   = {OMEGA_B_H2}")
print(f"  Omega_CDM h^2 = {OMEGA_CDM_H2}")
print(f"  Ratio         = {OMEGA_RATIO:.2f}")
print()

def calculate_ratio(g_total, g_active=G_SM):
    """
    Simple model: ratio ~ (g_arch / g_active) * entropy_factor
    """
    if g_total <= g_active:
        return 0
    g_arch = g_total - g_active
    entropy_factor = (11/4)**(1/3)  # ~1.4, from e+e- annihilation
    return (g_arch / g_active) * entropy_factor

# Find required g_total
g_range = np.linspace(G_SM + 1, 1500, 1000)
ratios = [calculate_ratio(g) for g in g_range]
idx = np.argmin(np.abs(np.array(ratios) - OMEGA_RATIO))
g_required = g_range[idx]
f_required = (g_required - G_SM) / g_required

print("="*60)
print("MODEL RESULT")
print("="*60)
print(f"\nTo get Omega_DM/Omega_b = {OMEGA_RATIO:.2f}:")
print(f"  Required g_total = {g_required:.0f}")
print(f"  Decoupling fraction f = {f_required:.1%}")
print(f"  Archival DOF g_arch = {g_required - G_SM:.0f}")
print()

# BSM scenarios table
print("="*60)
print("BSM SCENARIOS")
print("="*60)
print(f"\n{'Model':<25} {'g_total':<10} {'Ratio':<10}")
print("-"*45)

scenarios = [
    ("SM only", G_SM),
    ("MSSM", 228.75),
    ("SO(10) GUT", 400),
    ("E6 GUT", 500),
    ("String (moderate)", 700),
    ("Required", g_required),
]

for name, g in scenarios:
    r = calculate_ratio(g)
    print(f"{name:<25} {g:<10.0f} {r:<10.2f}")

# Plot
plt.figure(figsize=(10, 6))
plt.plot(g_range, ratios, 'b-', linewidth=2, label='Model')
plt.axhline(OMEGA_RATIO, color='r', linestyle='--', linewidth=2, label=f'Planck: {OMEGA_RATIO:.2f}')
plt.axvline(g_required, color='g', linestyle=':', linewidth=2, label=f'Required: g_*={g_required:.0f}')
plt.xlabel('Total degrees of freedom g_*', fontsize=12)
plt.ylabel('Omega_DM / Omega_b', fontsize=12)
plt.title('Dark-to-Baryon Ratio from Holographic Decoupling', fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)
plt.xlim(100, 1500)
plt.ylim(0, 15)

os.makedirs('../results/figures', exist_ok=True)
plt.savefig('../results/figures/dm_baryon_ratio.png', dpi=150)
plt.savefig('../results/figures/dm_baryon_ratio.pdf', dpi=300)
print(f"\nSaved: ../results/figures/dm_baryon_ratio.png")

plt.show()

# Save
results = {
    'OMEGA_RATIO_observed': OMEGA_RATIO,
    'g_SM': G_SM,
    'g_required': float(g_required),
    'f_decoupled': float(f_required),
}
with open('../results/03_dm_baryon_ratio.json', 'w') as f:
    json.dump(results, f, indent=2)

print("Saved: ../results/03_dm_baryon_ratio.json")