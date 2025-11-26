#!/usr/bin/env python3
"""
run_all.py - Run all Dark Matter calculations
"""

import subprocess
import sys
import os

scripts = [
    '01_saturation_temperature.py',
    '02_dof_evolution.py',
    '03_dm_baryon_ratio.py',
    '04_sparc_analysis.py',
    '05_entropy_budget.py',
]

def main():
    print("="*60)
    print("DARK MATTER FROM HOLOGRAPHIC SATURATION")
    print("Running all calculations...")
    print("="*60)
    
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    for script in scripts:
        print(f"\n{'='*60}")
        print(f">>> {script}")
        print('='*60)
        subprocess.run([sys.executable, script])
    
    print("\n" + "="*60)
    print("ALL COMPLETE")
    print("="*60)
    print("\nResults: ../results/")
    print("Figures: ../results/figures/")

if __name__ == '__main__':
    main()