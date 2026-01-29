# scripts/lgf/compute_zlgf_fcc.py
"""
FCC lattice Green function in ORT normalization.

Symbol of transition operator P on FCC (12-neighbor, normalized):
  P̂(k) = (1/3) [cos k1 cos k2 + cos k2 cos k3 + cos k3 cos k1],  k in [-pi, pi]^3

LGF:
  G(z) = (1/(2π)^3) ∫_{[-π,π]^3} 1 / (1 - z P̂(k)) d^3k

Using even symmetry in each variable:
  G(z) = (1/π^3) ∫_{[0,π]^3} 1 / (1 - z P̂(k)) d^3k

Target:
  Z_LGF = lim_{z->1-} G(z)   (finite for 3D FCC)

Near z=1-, we fit:
  G(1-eps) = Z0 + a1 eps^(1/2) + a2 eps + a3 eps^(3/2) + ...
"""

from __future__ import annotations

import argparse
import json
import time
import mpmath as mp


def phat(k1, k2, k3):
    c1 = mp.cos(k1)
    c2 = mp.cos(k2)
    c3 = mp.cos(k3)
    return (c1 * c2 + c2 * c3 + c3 * c1) / mp.mpf(3)


def G_fcc(z, dps=50, split="0.15"):
    """
    Compute G(z) by nested adaptive quadrature on [0, pi]^3.

    split: each axis is split into [0, split] and [split, pi] to stabilize
           integration as z -> 1- (integrable near-singularity at k=0).
    """
    mp.mp.dps = int(dps)

    one = mp.mpf(1)
    pi = mp.pi

    z = mp.mpf(z)
    if z < 0:
        raise ValueError("z must be >= 0")

    # Exact shortcut: G(0) = 1
    if z == 0:
        return one

    # Do not allow exact z=1 in numeric integration (k=0 would hit division by zero)
    if z >= 1:
        z = one - mp.mpf("1e-30")

    split = mp.mpf(split)
    if not (mp.mpf("0") < split < pi):
        raise ValueError("split must satisfy 0 < split < pi")

    intervals = [mp.mpf("0.0"), split, pi]

    def integrand(k1, k2, k3):
        return one / (one - z * phat(k1, k2, k3))

    def inner_k3(k1, k2):
        return mp.quadts(lambda k3: integrand(k1, k2, k3), intervals)

    def inner_k2(k1):
        return mp.quadts(lambda k2: inner_k3(k1, k2), intervals)

    I = mp.quadts(lambda k1: inner_k2(k1), intervals)

    # (1/(2π)^3)*8 over [0,π]^3 equals 1/π^3
    return I / (pi ** 3)


def solve_fit_half_integer(eps_vals, G_vals, halfint_order):
    """
    Fit:
      G(eps) = Z0 + a1 eps^(1/2) + ... + aN eps^(N/2)
    where N = halfint_order.

    If number of points equals N+1 -> solve square system via LU.
    Else -> least squares via normal equations.
    """
    mmax = int(halfint_order)
    if mmax < 0:
        raise ValueError("halfint_order must be >= 0")

    X_rows = []
    for eps in eps_vals:
        row = [eps ** (mp.mpf(m) / 2) for m in range(0, mmax + 1)]
        X_rows.append(row)

    X = mp.matrix(X_rows)
    y = mp.matrix(G_vals)

    if X.rows == X.cols:
        beta = mp.lu_solve(X, y)
    else:
        beta = (X.T * X) ** -1 * (X.T * y)

    return beta


def extrapolate_Zlgf(eps_list, dps=80, split="0.15", halfint_order=4):
    """
    Compute G(1-eps) for eps in eps_list and fit half-integer model to estimate Z_LGF.
    Returns: (Z0, coeffs, table)
    """
    mp.mp.dps = int(dps)
    one = mp.mpf(1)

    eps_vals = [mp.mpf(e) for e in eps_list]
    if any(eps <= 0 for eps in eps_vals):
        raise ValueError("All eps must be > 0")

    G_vals = []
    table = []

    for eps in eps_vals:
        z = one - eps
        t0 = time.time()
        val = G_fcc(z, dps=dps, split=split)
        dt = time.time() - t0

        G_vals.append(val)
        table.append(
            {
                "eps": mp.nstr(eps, 20),
                "z": mp.nstr(z, 20),
                "G": mp.nstr(val, 40),
                "seconds": dt,
            }
        )

    beta = solve_fit_half_integer(eps_vals, G_vals, halfint_order=halfint_order)
    Z0 = beta[0]
    coeffs = [beta[i] for i in range(beta.rows)]
    return Z0, coeffs, table


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dps", type=int, default=60, help="mpmath decimal digits")
    p.add_argument("--split", type=str, default="0.15", help="split point near 0")
    p.add_argument("--z", type=str, default=None, help="compute single G(z)")
    p.add_argument("--eps", nargs="*", default=None, help="eps list for extrapolation (z=1-eps)")
    p.add_argument("--halfint_order", type=int, default=4, help="max half-integer power index")
    p.add_argument("--out_json", type=str, default=None, help="optional output JSON path")
    p.add_argument("--quick", action="store_true", help="coarse eps list for quick run")
    args = p.parse_args()

    result = {
        "method": "fcc_lgf_kspace_integral",
        "dps": args.dps,
        "split": args.split,
        "halfint_order": args.halfint_order,
    }

    if args.z is not None:
        z = mp.mpf(args.z)
        t0 = time.time()
        val = G_fcc(z, dps=args.dps, split=args.split)
        dt = time.time() - t0
        result["single"] = {"z": str(z), "G": mp.nstr(val, 60), "seconds": dt}
        print("G(z) =", result["single"]["G"])
        print("seconds =", dt)

    else:
        if args.quick:
            eps_list = ["1e-3", "3e-4", "1e-4", "3e-5"]
        else:
            eps_list = args.eps or ["1e-3", "1e-4", "1e-5", "1e-6"]

        Z0, coeffs, table = extrapolate_Zlgf(
            eps_list=eps_list,
            dps=args.dps,
            split=args.split,
            halfint_order=args.halfint_order,
        )

        result["fit"] = {
            "eps_list": [str(e) for e in eps_list],
            "Z_LGF_estimate": mp.nstr(Z0, 60),
            "coeffs": [mp.nstr(c, 40) for c in coeffs],
            "table": table,
        }

        print("Z_LGF estimate =", result["fit"]["Z_LGF_estimate"])
        print("coeffs =", result["fit"]["coeffs"])

    if args.out_json:
        with open(args.out_json, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2)
        print("wrote", args.out_json)


if __name__ == "__main__":
    main()