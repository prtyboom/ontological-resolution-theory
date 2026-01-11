import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

from tabulate import tabulate

ROOT = Path(__file__).resolve().parent

# =========================
# Constants (MeV-based)
# =========================
U_MEV = 931.49410242          # 1 u in MeV/c^2 (CODATA-style)
M_P = 938.27208816            # proton mass (MeV/c^2)
M_N = 939.56542052            # neutron mass (MeV/c^2)

@dataclass(frozen=True)
class ORTConstants:
    k2D: int = 6
    k3D: int = 12

def Z_value() -> float:
    pi = math.pi
    return (pi**4) + 4*(pi**2) + (pi - 3.0) + (1.0/144.0)

def alpha_inv_ort(Z: float, k3D: int) -> float:
    # vacancy correction: (k3D-1)/(k3D^3 * Z)
    delta = (k3D - 1) / (k3D**3 * Z)
    return Z - delta

def mt_ort_GeV(me_MeV: float, Z: float, k2D: int, k3D: int) -> float:
    mt_MeV = me_MeV * (k2D + k3D) * (Z**2)
    return mt_MeV / 1e3  # MeV -> GeV

def v_exp_from_GF(GF: float) -> float:
    # v = (sqrt(2) * GF)^(-1/2)
    return 1.0 / math.sqrt(math.sqrt(2.0) * GF)

def v_ort_from_mt(mt_GeV: float, k3D: int) -> float:
    Z0D = 1.0 / (k3D**2)  # 1/144
    return math.sqrt(2.0) * mt_GeV * (1.0 + Z0D)

# =========================
# ORT mass-ratio grammar (legacy)
# =========================
def mp_me_ort() -> float:
    pi = math.pi
    return 6*(pi**5) + (5.0/144.0)

def mmu_me_ort(Z: float) -> float:
    pi = math.pi
    return 1.5*Z + 9.0*(pi - 3.0) - (9.0/144.0)

def mtau_me_ort(Z: float, k3D: int) -> float:
    pi = math.pi
    return (pi**3 - 6.0)*Z + 4.0*(pi**2) + (k3D - 1)

# =========================
# Nuclear sector helpers
# =========================
def eb_per_a_from_atomic_mass_u(Z: int, A: int, atomic_mass_u: float, me_MeV: float) -> float:
    """
    Eb/A (MeV) from atomic mass (u):
      M_atom = m_u * U_MEV
      M_nucleus ≈ M_atom - Z*m_e  (electron binding energy ignored; keV-scale)
      Eb = Z*M_P + (A-Z)*M_N - M_nucleus
    """
    M_atom = atomic_mass_u * U_MEV
    M_nucleus = M_atom - Z * me_MeV
    Eb = Z * M_P + (A - Z) * M_N - M_nucleus
    return Eb / A

def ebind_peak_from_dataset(path: Path, me_MeV: float):
    """
    Reads data/nubase_stable.csv and returns:
      (best_nuclide_str, best_EbA_MeV)
    """
    if not path.exists():
        return None, None

    best_nuclide = None
    best_val = None

    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            try:
                if int(r["is_stable"]) != 1:
                    continue
                Z = int(r["Z"])
                A = int(r["A"])
                sym = r["symbol"].strip()
                m_u = float(r["atomic_mass_u"])
            except Exception:
                continue

            eb_a = eb_per_a_from_atomic_mass_u(Z, A, m_u, me_MeV)
            nucl = f"{A}{sym}"

            if (best_val is None) or (eb_a > best_val):
                best_val = eb_a
                best_nuclide = nucl

    return best_nuclide, best_val

# =========================
# Utility
# =========================
def rel_err(pred: float, ref: float) -> float:
    return (pred - ref) / ref

def sigma_residual(abs_res: float, sigma):
    if sigma is None:
        return None
    try:
        s = float(sigma)
        if s == 0:
            return None
        return abs_res / s
    except Exception:
        return None

# =========================
# Main
# =========================
def main():
    # BOM-safe read
    ref = json.loads((ROOT / "reference_values.json").read_text(encoding="utf-8-sig"))
    V = ref["values"]
    notes = ref.get("source_notes", {})

    const = ORTConstants()
    pi = math.pi

    Z = Z_value()

    # Core tests
    alpha_ort = alpha_inv_ort(Z, const.k3D)

    mt_ort = mt_ort_GeV(V["me_MeV"], Z, const.k2D, const.k3D)
    v_exp = v_exp_from_GF(V["GF_GeVminus2"])
    v_ort = v_ort_from_mt(mt_ort, const.k3D)

    # Reference ratios from pinned masses
    mp_me_ref = V["m_p_MeV"] / V["me_MeV"]
    mmu_me_ref = V["m_mu_MeV"] / V["me_MeV"]
    mtau_me_ref = V["m_tau_MeV"] / V["me_MeV"]

    # ORT ratios
    mp_me_pred = mp_me_ort()
    mmu_me_pred = mmu_me_ort(Z)
    mtau_me_pred = mtau_me_ort(Z, const.k3D)

    # Nuclear peak (ORT prediction + pinned reference)
    ebA_ort = (const.k3D - 1) * (pi / 2.0) * V["me_MeV"]
    ebA_ref = V.get("Ebind_peak_MeV_perA", None)
    ebA_nucl_ref = V.get("Ebind_peak_nuclide", None)

    # Dataset diagnostic (max among stable rows)
    dataset_path = ROOT / "data" / "nubase_stable.csv"
    best_nuclide_ds, ebA_ds = ebind_peak_from_dataset(dataset_path, V["me_MeV"])

    # Table rows
    rows = []
    rows.append(["Z", Z, None, None, None])

    rows.append([
        "alpha_inv (ORT)",
        alpha_ort,
        V["alpha_inv"],
        alpha_ort - V["alpha_inv"],
        rel_err(alpha_ort, V["alpha_inv"])
    ])

    rows.append([
        "mp/me (ORT)",
        mp_me_pred,
        mp_me_ref,
        mp_me_pred - mp_me_ref,
        rel_err(mp_me_pred, mp_me_ref)
    ])

    rows.append([
        "mmu/me (ORT)",
        mmu_me_pred,
        mmu_me_ref,
        mmu_me_pred - mmu_me_ref,
        rel_err(mmu_me_pred, mmu_me_ref)
    ])

    rows.append([
        "mtau/me (ORT)",
        mtau_me_pred,
        mtau_me_ref,
        mtau_me_pred - mtau_me_ref,
        rel_err(mtau_me_pred, mtau_me_ref)
    ])

    rows.append([
        "mt_GeV (ORT)",
        mt_ort,
        V["mt_GeV"],
        mt_ort - V["mt_GeV"],
        rel_err(mt_ort, V["mt_GeV"])
    ])

    rows.append(["v_GeV (exp from GF)", v_exp, None, None, None])

    rows.append([
        "v_GeV (ORT from mt)",
        v_ort,
        v_exp,
        v_ort - v_exp,
        rel_err(v_ort, v_exp)
    ])

    # Nuclear line only if pinned reference exists
    if ebA_ref is not None:
        rows.append([
            "Ebind/A peak (ORT) [MeV]",
            ebA_ort,
            ebA_ref,
            ebA_ort - ebA_ref,
            rel_err(ebA_ort, ebA_ref)
        ])

    table = tabulate(
        rows,
        headers=["quantity", "predicted", "reference", "abs_residual", "relative_residual"],
        floatfmt=".10g"
    )
    print(table)

    # Additional diagnostics (sigma residuals where meaningful)
    mt_sigma = V.get("mt_GeV_sigma", None)
    mt_sig = sigma_residual(mt_ort - V["mt_GeV"], mt_sigma)

    out = {
        "ort": {
            "Z": Z,
            "alpha_inv_ort": alpha_ort,
            "mp_me_ort": mp_me_pred,
            "mmu_me_ort": mmu_me_pred,
            "mtau_me_ort": mtau_me_pred,
            "mt_ort_GeV": mt_ort,
            "v_exp_from_GF_GeV": v_exp,
            "v_ort_from_mt_GeV": v_ort,
            "EbA_peak_ort_MeV": ebA_ort
        },
        "reference_values": V,
        "source_notes": notes,
        "reference_derived": {
            "mp_me_ref": mp_me_ref,
            "mmu_me_ref": mmu_me_ref,
            "mtau_me_ref": mtau_me_ref
        },
        "nuclear_diagnostic": {
            "dataset_path": str(dataset_path),
            "dataset_best_nuclide": best_nuclide_ds,
            "dataset_best_Eb_per_A_MeV": ebA_ds,
            "pinned_ref_nuclide": ebA_nucl_ref,
            "pinned_ref_Eb_per_A_MeV": ebA_ref
        },
        "sigma_residuals": {
            "mt_sigma": mt_sigma,
            "mt_sigma_residual": mt_sig
        }
    }

    # Save JSON
    (ROOT / "results.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

    # Save Markdown report
    md = []
    md.append("# ORT Immediate Tests — Results\n")
    md.append("Reference values are pinned in `reference_values.json`.\n")
    md.append("```\n" + table + "\n```\n")
    (ROOT / "results.md").write_text("\n".join(md), encoding="utf-8")

    # Save CSV
    with (ROOT / "results.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "predicted", "reference", "abs_residual", "relative_residual"])
        for r in rows:
            w.writerow(r)

if __name__ == "__main__":
    main()