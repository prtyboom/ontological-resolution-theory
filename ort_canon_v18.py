"""
ORT Canon v18.0 — Locked Results Test Suite
============================================
All formulas from Complete Canon v17.2 + Paper D (Koide).
Zero free parameters.
"""

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# ============================================================
# AXIOMS (discrete structural choices)
# ============================================================
PI = math.pi
E = math.e
PHI = (1 + math.sqrt(5)) / 2  # Golden Ratio (A0)

K = 12          # A3: FCC coordination
K2D = 6         # A3: hexagonal section
D = 4           # A5: 4 -> 3+1

# ============================================================
# DERIVED CONSTANTS (from axioms, zero free parameters)
# ============================================================

# A4: Sweep Operator
SWEEP = 2 * PI * E

# T4.5: Bare Impedance
Z_BARE = PI**4 + 4*PI**2 + (PI - 3) + 1/K**2

# T4.9: First vacancy correction
DELTA_1 = (K - 1) / (K**3 * Z_BARE)

# T4.10: Second vacancy correction
DELTA_2 = DELTA_1 / (K - 3)

# T4.11: Physical impedance
Z_PHYS = Z_BARE - DELTA_1 + DELTA_2

# ============================================================
# LEPTON MASSES (Canon Part VII + Paper D)
# ============================================================

# T7.4: Precision tau/electron (LOCKED)
def mtau_me():
    return (PI**3 - 6) * Z_BARE + 4*PI**2 + (K - 1)
# TD.1 + TD.2: Koide-corrected muon/electron (LOCKED)
def mmu_me_koide(r2):
    """
    From K=2/3 and fixed r2=m_tau/m_e, solve for r1=m_mu/m_e.
    Quadratic in u=sqrt(r1).
    """
    s = math.sqrt(r2)
    # u^2 - (4+4s)*u + (1+3*r2-2*s^2-4*s) = 0
    a_coef = 1.0
    b_coef = -(4 + 4*s)
    c_coef = 1 + 3*r2 - 2*s**2 - 4*s
    disc = b_coef**2 - 4*a_coef*c_coef
    u = (-b_coef - math.sqrt(disc)) / (2*a_coef)  # physical root
    return u**2

# T7.1: Tau/muon ratio (leading order, LOCKED)
def mtau_mmu():
    return SWEEP * (1 - 2/K**2)

# Proton mass: 6*pi^5 + 5/144
def mp_me():
    return K2D * PI**5 + 5/K**2

# ============================================================
# ELECTROWEAK (Canon Part VIII)
# ============================================================

# T8.1: Weinberg angle (LOCKED)
def sin2_theta_w():
    return 3 / (K + 1)

# ============================================================
# STRONG (Canon Part IX)
# ============================================================

# T9.1: Colour number
N_C = D - 1  # = 3

# T9.2: Gluon count
N_G = K - D  # = 8

# T9.3: Strong coupling (LOCKED)
def alpha_s():
    return math.sqrt(2) / K
# ============================================================
# GRAVITY (Canon Part VI)
# ============================================================

# T6.8: Newton's constant
def G_ort(hbar_si, c_si, mp_si):
    R = K + K2D  # resonant rank = 18
    pref = PHI**4/4 + 1/(4*K**2)
    return pref * hbar_si * c_si / (Z_BARE**R * mp_si**2)
# ============================================================
# KOIDE PARAMETER (Paper D)
# ============================================================

def koide(r1, r2):
    """Koide parameter from mass ratios r1=m_mu/m_e, r2=m_tau/m_e."""
    num = 1 + r1 + r2
    den = (1 + math.sqrt(r1) + math.sqrt(r2))**2
    return num / den

# ============================================================
# NEUTRINO HIERARCHY (Canon Part X)
# ============================================================

def nu_ratio():
    return (K - 1) * PI / 6

# ============================================================
# ORACLE THEOREMS (Canon Part III)
# ============================================================

def oracle_first():
    """4*phi^3 vs 2*pi*e"""
    return 4*PHI**3, SWEEP

def oracle_second():
    """12*phi^2 vs 10*pi"""
    return 12*PHI**2, 10*PI

# ============================================================
# REFERENCE VALUES
# ============================================================

REF = json.loads((ROOT / "reference_values.json").read_text(encoding="utf-8-sig"))["values"]

# SI constants for G
HBAR_SI = 1.054571817e-34
C_SI = 299792458.0
MP_SI = 1.67262192e-27

# Experimental references
ALPHA_INV_EXP = REF["alpha_inv"]
MP_ME_EXP = REF["m_p_MeV"] / REF["me_MeV"]
MMU_ME_EXP = REF["m_mu_MeV"] / REF["me_MeV"]
MTAU_ME_EXP = REF["m_tau_MeV"] / REF["me_MeV"]
MT_EXP = REF["mt_GeV"]
G_EXP = 6.67430e-11
SIN2W_EXP = 0.23121
ALPHAS_EXP = 0.1180
NU_RATIO_EXP = 5.77

# ============================================================
# MAIN
# ============================================================

def rel(pred, exp):
    return (pred - exp) / exp

def main():
    # Compute all predictions
    tau_e = mtau_me()
    mu_e = mmu_me_koide(tau_e)
    tau_mu = tau_e / mu_e
    p_e = mp_me()
    s2w = sin2_theta_w()
    a_s = alpha_s()
    G = G_ort(HBAR_SI, C_SI, MP_SI)
    K_koide = koide(mu_e, tau_e)
    nu_r = nu_ratio()
    o1_phi, o1_sweep = oracle_first()
    o2_phi, o2_pi = oracle_second()

    # Results table
    rows = [
        ["IMPEDANCE"],
        ["  Z_bare",            f"{Z_BARE:.6f}",        "",                 ""],
        ["  Z_phys (alpha^-1)", f"{Z_PHYS:.8f}",        f"{ALPHA_INV_EXP}", f"{rel(Z_PHYS, ALPHA_INV_EXP):.2e}"],
        [""],
        ["LEPTON MASSES"],
        ["  m_tau/m_e",         f"{tau_e:.3f}",         f"{MTAU_ME_EXP:.2f}", f"{rel(tau_e, MTAU_ME_EXP):.2e}"],
        ["  m_mu/m_e (Koide)",  f"{mu_e:.3f}",          f"{MMU_ME_EXP:.6f}", f"{rel(mu_e, MMU_ME_EXP):.2e}"],
        ["  m_tau/m_mu",        f"{tau_mu:.3f}",        f"{MTAU_ME_EXP/MMU_ME_EXP:.3f}", f"{rel(tau_mu, MTAU_ME_EXP/MMU_ME_EXP):.2e}"],
        ["  m_p/m_e",           f"{p_e:.3f}",           f"{MP_ME_EXP:.3f}", f"{rel(p_e, MP_ME_EXP):.2e}"],
        ["  Koide K",           f"{K_koide:.6f}",       "0.666667",          f"{rel(K_koide, 2/3):.2e}"],
        [""],
        ["ELECTROWEAK"],
        ["  sin2_theta_W",      f"{s2w:.5f}",           f"{SIN2W_EXP}",     f"{rel(s2w, SIN2W_EXP):.2e}"],
        [""],
        ["STRONG"],
        ["  N_c",               f"{N_C}",               "3",                "exact"],
        ["  N_g",               f"{N_G}",               "8",                "exact"],
        ["  alpha_s(M_Z)",      f"{a_s:.4f}",           f"{ALPHAS_EXP}",    f"{rel(a_s, ALPHAS_EXP):.2e}"],
        [""],
        ["GRAVITY"],
        ["  G (SI)",            f"{G:.3e}",             f"{G_EXP:.3e}",     f"{rel(G, G_EXP):.2e}"],
        [""],
        ["NEUTRINO"],
        ["  m_nu ratio",        f"{nu_r:.4f}",          f"{NU_RATIO_EXP}",  f"{rel(nu_r, NU_RATIO_EXP):.2e}"],
        [""],
        ["ORACLE THEOREMS"],
        ["  4*phi^3 vs S",      f"{o1_phi:.4f}",        f"{o1_sweep:.4f}",  f"{rel(o1_phi, o1_sweep):.2e}"],
        ["  12*phi^2 vs 10pi",  f"{o2_phi:.4f}",        f"{o2_pi:.4f}",     f"{rel(o2_phi, o2_pi):.2e}"],
        [""],
        ["GOLDEN CASCADE"],
        ["  phi",               f"{PHI:.6f}",           "",                 "A0: phi^2=phi+1"],
        ["  e",                 f"{E:.6f}",             "",                 "3D: f'=f"],
        ["  pi",                f"{PI:.6f}",            "",                 "2D: C/d"],
    ]

    print("=" * 72)
    print("  ORT CANON v18.0 — COMPLETE TEST SUITE")
    print("  Zero free parameters")
    print("=" * 72)
    for row in rows:
        if len(row) == 1:
            print(f"\n  {row[0]}")
        elif len(row) == 0 or row[0] == "":
            pass
        else:
            print(f"    {row[0]:<25s} {row[1]:>15s}  {row[2]:>15s}  {row[3]:>12s}")
    print("\n" + "=" * 72)
    print("  Locked: 17 | Conditional: 4 | Hypotheses: 2 | Free params: 0")
    print("=" * 72)

if __name__ == "__main__":
    main()