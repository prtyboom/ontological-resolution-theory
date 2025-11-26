#!/usr/bin/env python3
"""
02_dof_evolution.py - Effective degrees of freedom g_*(T)
"""

import numpy as np
import matplotlib.pyplot as plt
import json
import os

# SM particles: (name, mass_GeV, bosonic_dof, fermionic_dof)
SM_PARTICLES = [
    ('photon', 0, 2, 0),
    ('gluon', 0, 16, 0),
    ('W+', 80.4, 3, 0),
    ('W-', 80.4, 3, 0),
    ('Z', 91.2, 3, 0),
    ('Higgs', 125, 1, 0),
    ('electron', 0.000511, 0, 4),
    ('muon', 0.106, 0, 4),
    ('tau', 1.78, 0, 4),
    ('nu_e', 0, 0, 2),
    ('nu_mu', 0, 0, 2),
    ('nu_tau', 0, 0, 2),
    ('up', 0.002, 0, 12),
    ('down', 0.005, 0, 12),
    ('strange', 0.095, 0, 12),
    ('charm', 1.27, 0, 12),
    ('bottom', 4.18, 0, 12),
    ('top', 172.5, 0, 12),
]

def g_star(T_GeV):
    """Effective relativistic degrees of freedom."""
    if T_GeV <= 0:
        return 0
    g = 0
    for name, mass, g_b, g_f in SM_PARTICLES:
        x = mass / T_GeV
        if x < 0.1:
            weight = 1.0
        elif x > 10:
            weight = 0.0
        else:
            weight = np.exp(-x)
        g += g_b * weight + (7/8) * g_f * weight
    return g

# Verify
g_high = g_star(1e6)
g_expected = sum(g_b + (7/8)*g_f for _, _, g_b, g_f in SM_PARTICLES)

print("="*60)
print("DEGREES OF FREEDOM g_*(T)")
print("="*60)
print(f"\nHigh-T limit: g_* = {g_high:.2f}")
print(f"Expected:     g_* = {g_expected:.2f}")
print()

# Table
print(f"{'T (GeV)':<12} {'g_*':<10}")
print("-"*22)
for T in [1e6, 1000, 100, 10, 1, 0.1, 0.01, 0.001]:
    print(f"{T:<12.0e} {g_star(T):<10.2f}")

# Plot
T_range = np.logspace(-4, 6, 500)
g_values = [g_star(T) for T in T_range]

plt.figure(figsize=(10, 6))
plt.semilogx(T_range, g_values, 'b-', linewidth=2)
plt.axhline(106.75, color='r', linestyle='--', label='Full SM = 106.75')
plt.xlabel('Temperature T [GeV]', fontsize=12)
plt.ylabel('g_*', fontsize=12)
plt.title('Standard Model Degrees of Freedom', fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)
plt.xlim(1e-4, 1e6)
plt.ylim(0, 120)

os.makedirs('../results/figures', exist_ok=True)
plt.savefig('../results/figures/g_star_evolution.png', dpi=150)
plt.savefig('../results/figures/g_star_evolution.pdf', dpi=300)
print("\nSaved: ../results/figures/g_star_evolution.png")

plt.show()

# Save data
with open('../results/02_dof_evolution.json', 'w') as f:
    json.dump({'g_SM': 106.75, 'g_high_T': g_high}, f, indent=2)