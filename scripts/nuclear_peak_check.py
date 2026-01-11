import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "nubase_stable.csv"
OUT = ROOT / "results_nuclear_peak.json"

# constants (MeV)
U_MEV = 931.49410242
M_P = 938.27208816
M_N = 939.56542052
M_E = 0.51099895

K3D = 12
EB_PER_A_ORT = (K3D - 1) * (math.pi / 2.0) * M_E  # MeV

def eb_per_a_from_atomic_mass_u(Z: int, A: int, atomic_mass_u: float) -> float:
    M_atom = atomic_mass_u * U_MEV
    M_nucleus = M_atom - Z * M_E  # ignore electron binding energy (keV-scale)
    Eb = Z * M_P + (A - Z) * M_N - M_nucleus
    return Eb / A

def read_rows(path: Path):
    rows = []
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            Z = int(r["Z"])
            A = int(r["A"])
            sym = r["symbol"].strip()
            is_stable = int(r["is_stable"])
            m_u = float(r["atomic_mass_u"])
            rows.append((Z, A, sym, is_stable, m_u))
    return rows

def main():
    if not DATA.exists():
        raise FileNotFoundError(f"Missing dataset: {DATA}")

    rows = read_rows(DATA)
    stable = [r for r in rows if r[3] == 1]
    if not stable:
        raise RuntimeError("No stable nuclides in dataset (is_stable=1).")

    scored = []
    for (Z, A, sym, is_stable, m_u) in stable:
        eb_a = eb_per_a_from_atomic_mass_u(Z, A, m_u)
        scored.append({
            "Z": Z,
            "A": A,
            "symbol": sym,
            "nuclide": f"{A}{sym}",
            "Eb_per_A_MeV": eb_a
        })

    best = max(scored, key=lambda x: x["Eb_per_A_MeV"])
    abs_res = best["Eb_per_A_MeV"] - EB_PER_A_ORT
    rel_res = abs_res / EB_PER_A_ORT

    report = {
        "ort_prediction": {
            "Eb_per_A_peak_MeV": EB_PER_A_ORT,
            "formula": "(k3D-1)*(pi/2)*m_e",
            "k3D": K3D,
            "m_e_MeV": M_E
        },
        "dataset": {
            "path": str(DATA),
            "row_count": len(rows),
            "stable_count": len(stable),
            "note": "Template dataset; replace with pinned AME/NUBASE stable nuclides list."
        },
        "best_stable_nuclide_in_dataset": best,
        "residuals": {
            "abs_MeV": abs_res,
            "rel": rel_res
        }
    }

    OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print("=== ORT OP.4 Nuclear Peak Check ===")
    print(f"Dataset: {DATA}  (stable rows: {len(stable)})")
    print(f"ORT prediction (Eb/A)_max = {EB_PER_A_ORT:.6f} MeV")
    print(f"Best in dataset: {best['nuclide']}  Eb/A = {best['Eb_per_A_MeV']:.6f} MeV")
    print(f"Residual: {abs_res:+.6f} MeV  (rel {rel_res:+.3%})")
    print(f"Wrote: {OUT}")

if __name__ == "__main__":
    main()