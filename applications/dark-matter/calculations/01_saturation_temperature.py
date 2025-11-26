#!/usr/bin/env python3
"""
01_saturation_temperature.py (FIXED)

Holographic saturation temperature calculation.
"""

import numpy as np
from scipy import constants
import json
import os

# Constants
c = constants.c
hbar = constants.hbar
G = constants.G
k_B = constants.k

l_P = np.sqrt(hbar * G / c**3)
M_P = np.sqrt(hbar * c / G)
E_P_GeV = M_P * c**2 / constants.eV / 1e9

print("="*60)
print("SATURATION TEMPERATURE")
print("="*60)
print(f"\nPlanck energy: E_P = {E_P_GeV:.3e} GeV\n")

def T_sat(g_star):
    """T_sat ~ M_P / sqrt(g_*) from holographic bound."""
    return E_P_GeV / np.sqrt(g_star)

g_SM = 106.75
T_sat_SM = T_sat(g_SM)

print(f"Standard Model (g_* = {g_SM}):")
print(f"  T_sat = {T_sat_SM:.2e} GeV")
print()

print("NOTE: T_sat ~ 10^18 GeV is ABOVE GUT scale.")
print("The 'GUT scale' (10^16 GeV) is the REHEATING temperature,")
print("which we take as INPUT for our model.")
print()

# Save
results = {
    'E_P_GeV': E_P_GeV,
    'g_SM': g_SM,
    'T_sat_GeV': T_sat_SM,
    'T_reheat_assumed_GeV': 1e16,
}

os.makedirs('../results', exist_ok=True)
with open('../results/01_saturation_temperature.json', 'w') as f:
    json.dump(results, f, indent=2)

print("Saved: ../results/01_saturation_temperature.json")