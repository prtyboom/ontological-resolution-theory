import csv
import argparse
from pathlib import Path

def read_csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(line for line in f if not line.strip().startswith("#"))
        return list(reader)

def key_ZA(row):
    return (int(row["Z"]), int(row["A"]))

def main():
    ap = argparse.ArgumentParser(description="Join AME atomic masses with NUBASE stability into ORT nubase_stable.csv")
    ap.add_argument("--ame", required=True, help="Path to data/ame_atomic_masses.csv")
    ap.add_argument("--nubase", required=True, help="Path to data/nubase_stability.csv")
    ap.add_argument("--out", required=True, help="Path to output data/nubase_stable.csv")
    args = ap.parse_args()

    ame_path = Path(args.ame)
    nubase_path = Path(args.nubase)
    out_path = Path(args.out)

    ame = read_csv(ame_path)
    nubase = read_csv(nubase_path)

    ame_map = {key_ZA(r): r for r in ame}
    nub_map = {key_ZA(r): r for r in nubase}

    common = sorted(set(ame_map.keys()) & set(nub_map.keys()))
    missing_in_ame = sorted(set(nub_map.keys()) - set(ame_map.keys()))
    missing_in_nub = sorted(set(ame_map.keys()) - set(nub_map.keys()))

    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Z", "A", "symbol", "is_stable", "atomic_mass_u", "source_ame", "source_nubase"])
        for (Z, A) in common:
            a = ame_map[(Z, A)]
            n = nub_map[(Z, A)]
            w.writerow([
                Z,
                A,
                a.get("symbol", "").strip(),
                int(n["is_stable"]),
                float(a["atomic_mass_u"]),
                a.get("source", "").strip(),
                n.get("source", "").strip()
            ])

    print("=== ORT build_nubase_stable_csv ===")
    print(f"AME rows: {len(ame)}   NUBASE rows: {len(nubase)}")
    print(f"Joined rows: {len(common)}")
    print(f"Missing in AME: {len(missing_in_ame)}")
    print(f"Missing in NUBASE: {len(missing_in_nub)}")
    print(f"Wrote: {out_path}")

if __name__ == "__main__":
    main()