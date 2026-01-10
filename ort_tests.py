import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
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
    # Standard electroweak relation: v = (sqrt(2) * GF)^(-1/2)
    # in natural units (GeV)
    return 1.0 / math.sqrt(math.sqrt(2.0) * GF)

def v_ort_from_mt(mt_GeV: float, k3D: int) -> float:
    Z0D = 1.0 / (k3D**2)  # 1/144
    return math.sqrt(2.0) * mt_GeV * (1.0 + Z0D)

def rel_err(pred: float, ref: float) -> float:
    return (pred - ref) / ref

def main():
    ref = json.loads((ROOT / "reference_values.json").read_text(encoding="utf-8-sig"))
    V = ref["values"]

    const = ORTConstants()

    Z = Z_value()
    alpha_ort = alpha_inv_ort(Z, const.k3D)
    mt_ort = mt_ort_GeV(V["me_MeV"], Z, const.k2D, const.k3D)

    v_exp = v_exp_from_GF(V["GF_GeVminus2"])
    v_ort = v_ort_from_mt(mt_ort, const.k3D)

    rows = [
        ["Z", Z, None, None, None],
        ["alpha_inv (ORT)", alpha_ort, V["alpha_inv"], alpha_ort - V["alpha_inv"], rel_err(alpha_ort, V["alpha_inv"])],
        ["mt_GeV (ORT)", mt_ort, V["mt_GeV"], mt_ort - V["mt_GeV"], rel_err(mt_ort, V["mt_GeV"])],
        ["v_GeV (exp from GF)", v_exp, None, None, None],
        ["v_GeV (ORT from mt)", v_ort, v_exp, v_ort - v_exp, rel_err(v_ort, v_exp)],
    ]

    table = tabulate(
        rows,
        headers=["quantity", "predicted", "reference", "abs_residual", "relative_residual"],
        floatfmt=".10g"
    )
    print(table)

    # Save JSON
    out = {
        "ort": {
            "Z": Z,
            "alpha_inv_ort": alpha_ort,
            "mt_ort_GeV": mt_ort,
            "v_exp_from_GF_GeV": v_exp,
            "v_ort_from_mt_GeV": v_ort
        },
        "reference": V,
        "notes": ref.get("source_notes", {}),
        "residuals": {
            "alpha_inv_abs": alpha_ort - V["alpha_inv"],
            "alpha_inv_rel": rel_err(alpha_ort, V["alpha_inv"]),
            "mt_abs": mt_ort - V["mt_GeV"],
            "mt_rel": rel_err(mt_ort, V["mt_GeV"]),
            "v_abs": v_ort - v_exp,
            "v_rel": rel_err(v_ort, v_exp)
        }
    }
    (ROOT / "results.json").write_text(json.dumps(out, indent=2), encoding="utf-8-sig")

    # Save Markdown report
    md = []
    md.append("# ORT Immediate Tests вЂ” Results\n")
    md.append("Reference values are pinned in `reference_values.json`.\n")
    md.append("```\n" + table + "\n```\n")
    (ROOT / "results.md").write_text("\n".join(md), encoding="utf-8-sig")

if __name__ == "__main__":
    main()
