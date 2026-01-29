# scripts/lgf/compute_zlgf_fcc_grid.py
from __future__ import annotations

import argparse
import json
import time
import numpy as np


def phat_from_cos(c1, c2, c3):
    # P̂(k) = (1/3)(c1 c2 + c2 c3 + c3 c1)
    return (c1 * c2 + c2 * c3 + c3 * c1) / 3.0


def make_grid(N: int, power: int = 2):
    """
    Integrate over k in [0, pi]^3 using midpoint product rule with variable transform:
        k = pi * u^power,  u in [0,1]
    Then dk = pi * power * u^(power-1) du
    and with ORT symmetry-reduced form:
        G(z) = (1/pi^3) ∫_{[0,pi]^3} f(k) d^3k
    we get:
        G(z) = power^3 ∫_{[0,1]^3} f(pi u^power) u1^(power-1)u2^(power-1)u3^(power-1) du^3
    Midpoint discretization:
        G(z) ≈ (power^3 / N^3) Σ f(pi u^power) * w1*w2*w3
    where w = u^(power-1).
    """
    N = int(N)
    power = int(power)
    if N <= 0:
        raise ValueError("N must be positive")
    if power <= 0:
        raise ValueError("power must be positive integer")

    u = (np.arange(N, dtype=np.float64) + 0.5) / N
    k = np.pi * (u ** power)
    c = np.cos(k)
    w = u ** (power - 1)  # jacobian weight per axis (without constants)
    pref = (power ** 3) / (N ** 3)
    return c, w, pref


def G_fcc_grid(z: float, N: int = 140, power: int = 2) -> float:
    """
    FCC lattice Green function (ORT normalization), grid quadrature.

    P̂(k) = (1/3)(cos k1 cos k2 + cos k2 cos k3 + cos k3 cos k1)
    G(z)  = (1/pi^3) ∫_{[0,pi]^3} 1/(1 - z P̂(k)) d^3k

    Uses midpoint product quadrature on transformed grid to emphasize k≈0.
    Complexity ~ O(N^3) (but inner axis is vectorized).
    """
    z = float(z)
    if z < 0.0:
        raise ValueError("z must be >= 0")
    if z == 0.0:
        return 1.0
    if z >= 1.0:
        z = 1.0 - 1e-15  # keep strictly < 1

    c, w, pref = make_grid(N, power=power)

    s = 0.0
    c3 = c
    w3 = w

    # double loop + vectorized sum over k3
    for i in range(N):
        ci = c[i]
        wi = w[i]
        for j in range(N):
            cj = c[j]
            wij = wi * w[j]
            denom = 1.0 - z * phat_from_cos(ci, cj, c3)
            s += wij * np.sum(w3 / denom)

    return pref * s


def fit_half_integer(eps_list, G_list, halfint_order: int = 4):
    """
    Fit:
        G(1-eps) = Z0 + a1 eps^(1/2) + a2 eps + ... + aN eps^(N/2)
    by least squares.
    """
    eps = np.array(eps_list, dtype=np.float64)
    y = np.array(G_list, dtype=np.float64)

    cols = []
    for m in range(0, halfint_order + 1):
        cols.append(eps ** (0.5 * m))
    X = np.stack(cols, axis=1)

    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta


def estimate_Z_lgf(eps_list, N: int, power: int, halfint_order: int):
    table = []
    G_list = []

    for eps in eps_list:
        z = 1.0 - float(eps)
        t0 = time.time()
        Gv = G_fcc_grid(z, N=N, power=power)
        dt = time.time() - t0
        G_list.append(Gv)
        table.append({"eps": float(eps), "z": float(z), "G": float(Gv), "seconds": dt})

    beta = fit_half_integer(eps_list, G_list, halfint_order=halfint_order)
    Z0 = float(beta[0])
    return Z0, beta, table


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--z", type=float, default=None, help="compute single G(z)")
    ap.add_argument("--N", type=int, default=140, help="grid size per dimension")
    ap.add_argument("--power", type=int, default=2, help="k = pi*u^power (power>=1)")
    ap.add_argument("--eps", nargs="*", default=None, help="eps list for z=1-eps (extrapolation)")
    ap.add_argument("--halfint_order", type=int, default=4, help="max half-integer power index")
    ap.add_argument("--quick", action="store_true", help="use coarse eps list")
    ap.add_argument("--out_json", type=str, default=None, help="write JSON report")
    args = ap.parse_args()

    result = {
        "engine": "numpy_grid",
        "N": args.N,
        "power": args.power,
        "halfint_order": args.halfint_order,
    }

    if args.z is not None:
        t0 = time.time()
        val = G_fcc_grid(args.z, N=args.N, power=args.power)
        dt = time.time() - t0
        result["single"] = {"z": float(args.z), "G": float(val), "seconds": dt}
        print("G(z) =", result["single"]["G"])
        print("seconds =", dt)

    else:
        if args.quick:
            eps_list = [1e-3, 3e-4, 1e-4, 3e-5]
        else:
            eps_list = [float(e) for e in (args.eps or ["1e-3", "1e-4", "1e-5", "1e-6"])]

        Z0, beta, table = estimate_Z_lgf(
            eps_list=eps_list,
            N=args.N,
            power=args.power,
            halfint_order=args.halfint_order,
        )

        result["fit"] = {
            "eps_list": eps_list,
            "Z_LGF_estimate": Z0,
            "coeffs": [float(x) for x in beta],
            "table": table,
        }

        print("Z_LGF estimate =", result["fit"]["Z_LGF_estimate"])
        print("coeffs =", result["fit"]["coeffs"])
        for row in table:
            print("eps=", row["eps"], "G=", row["G"], "sec=", row["seconds"])

    if args.out_json:
        with open(args.out_json, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2)
        print("wrote", args.out_json)


if __name__ == "__main__":
    main()