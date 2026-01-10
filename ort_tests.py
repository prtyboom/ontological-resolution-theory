import json
import math
from dataclasses import dataclass
from pathlib import Path

from tabulate import tabulate

ROOT = Path(__file__).resolve().parent

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

# ORT mass-ratio grammar (legacy relations)
def mp_me_ort() -> float:
    pi = math.pi
    return 6*(pi**5) + (5.0/144.0)

def mmu_me_ort(Z: float) -> float:
    pi = math.pi
    return 1.5*Z + 9.0*(pi - 3.0) - (9.0/144.0)

def mtau_me_ort(Z: float, k3D: int) -> float:
    pi = math.pi
    return (pi**3 - 6.0)*Z + 4.0*(pi**2) + (k3D - 1)

def rel_err(pred: float, ref: float) -> float:
    return (pred - ref) / ref

def main():
    # BOM-safe read
    ref = json.loads((ROOT / "reference_values.json").read_text(encoding="utf-8-sig"))
    V = ref["values"]
    notes = ref.get("source_notes", {})

    const = ORTConstants()
    pi = math.pi

    Z = Z_value()

    # alpha
    alpha_ort = alpha_inv_ort(Z, const.k3D)

    # top + VEV
    mt_ort = mt_ort_GeV(V["me_MeV"], Z, const.k2D, const.k3D)
    v_exp = v_exp_from_GF(V["GF_GeVminus2"])
    v_ort = v_ort_from_mt(mt_ort, const.k3D)

    # reference ratios from masses
    mp_me_ref = V["m_p_MeV"] / V["me_MeV"]
    mmu_me_ref = V["m_mu_MeV"] / V["me_MeV"]
    mtau_me_ref = V["m_tau_MeV"] / V["me_MeV"]

    # ORT ratios
    mp_me_pred = mp_me_ort()
    mmu_me_pred = mmu_me_ort(Z)
    mtau_me_pred = mtau_me_ort(Z, const.k3D)

    rows = [
        ["Z", Z, None, None, None],

        ["alpha_inv (ORT)", alpha_ort, V["alpha_inv"],
         alpha_ort - V["alpha_inv"], rel_err(alpha_ort, V["alpha_inv"])],

        ["mp/me (ORT)", mp_me_pred, mp_me_ref,
         mp_me_pred - mp_me_ref, rel_err(mp_me_pred, mp_me_ref)],

        ["mmu/me (ORT)", mmu_me_pred, mmu_me_ref,
         mmu_me_pred - mmu_me_ref, rel_err(mmu_me_pred, mmu_me_ref)],

        ["mtau/me (ORT)", mtau_me_pred, mtau_me_ref,
         mtau_me_pred - mtau_me_ref, rel_err(mtau_me_pred, mtau_me_ref)],

        ["mt_GeV (ORT)", mt_ort, V["mt_GeV"],
         mt_ort - V["mt_GeV"], rel_err(mt_ort, V["mt_GeV"])],

        ["v_GeV (exp from GF)", v_exp, None, None, None],

        ["v_GeV (ORT from mt)", v_ort, v_exp,
         v_ort - v_exp, rel_err(v_ort, v_exp)],
    ]

    table = tabulate(
        rows,
        headers=["quantity", "predicted", "reference", "abs_residual", "relative_residual"],
        floatfmt=".10g"
    )
    print(table)

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
            "Delta_step_MeV": 1000.0*(pi - 3.0)/144.0
        },
        "reference_values": V,
        "source_notes": notes,
        "reference_derived": {
            "mp_me_ref": mp_me_ref,
            "mmu_me_ref": mmu_me_ref,
            "mtau_me_ref": mtau_me_ref
        },
        "residuals": {
            "alpha_inv_abs": alpha_ort - V["alpha_inv"],
            "alpha_inv_rel": rel_err(alpha_ort, V["alpha_inv"]),
            "mp_me_abs": mp_me_pred - mp_me_ref,
            "mp_me_rel": rel_err(mp_me_pred, mp_me_ref),
            "mmu_me_abs": mmu_me_pred - mmu_me_ref,
            "mmu_me_rel": rel_err(mmu_me_pred, mmu_me_ref),
            "mtau_me_abs": mtau_me_pred - mtau_me_ref,
            "mtau_me_rel": rel_err(mtau_me_pred, mtau_me_ref),
            "mt_abs": mt_ort - V["mt_GeV"],
            "mt_rel": rel_err(mt_ort, V["mt_GeV"]),
            "v_abs": v_ort - v_exp,
            "v_rel": rel_err(v_ort, v_exp)
        }
    }

    (ROOT / "results.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

    md = []
    md.append("# ORT Immediate Tests — Results\n")
    md.append("Reference values are pinned in `reference_values.json`.\n")
    md.append("```\n" + table + "\n```\n")
    (ROOT / "results.md").write_text("\n".join(md), encoding="utf-8")

if __name__ == "__main__":
    main()