# FCC lattice Green function (Watson integral) in ORT normalization

This note documents a reproducible numerical estimate of the FCC lattice Green function
constant used as a micro-bridge target in ORT v9.0.

Status: numerical / computational note (not a closed-form derivation).

---

## 1. Definitions (normalization)

We use the FCC step set with 12 neighbors and the normalized transition operator
P = (1/12)A, whose Brillouin-zone symbol is

P_hat(k) = (1/3)(cos k1 cos k2 + cos k2 cos k3 + cos k3 cos k1),  k in [-pi,pi]^3.

The lattice Green function (LGF) is defined by

G(z) = tr_bar (I - zP)^(-1) = (1/(2pi)^3) integral_{[-pi,pi]^3} d^3k / (1 - z P_hat(k)).

By even symmetry in each variable,

G(z) = (1/pi^3) integral_{[0,pi]^3} d^3k / (1 - z P_hat(k)).

The FCC Watson/LGF constant in ORT normalization is

Z_LGF := lim_{z -> 1^-} G(z).

For 3D FCC (transient random walk), this limit is finite.

---

## 2. Numerical method (fast grid quadrature)

We compute G(z) using a midpoint product quadrature on [0,pi]^3 with a variable
transformation that concentrates grid points near k=0 (where the integrand is most
sensitive when z -> 1^-).

We parameterize each axis by

k = pi u^p,  u in [0,1],  p in {2,3}.

Then

G(z) = p^3 integral_{[0,1]^3} (u1^(p-1) u2^(p-1) u3^(p-1)) / (1 - z P_hat(pi u1^p, pi u2^p, pi u3^p)) du1 du2 du3,

and we approximate this integral by the midpoint rule on an N×N×N grid in u.

Implementation lives in:

- scripts/lgf/compute_zlgf_fcc_grid.py

This implementation is double precision (NumPy) and is intended as a fast, reproducible baseline.
A high-precision (mpmath) version exists in development but is significantly slower.

---

## 3. Extrapolation to z -> 1^- (half-integer fit)

We estimate Z_LGF by evaluating G(1-epsilon) at a list of small epsilon and fitting:

G(1-epsilon) ≈ Z0 + a1 epsilon^(1/2) + a2 epsilon + a3 epsilon^(3/2).

Rationale: near k=0, P_hat(k) ≈ 1 - ||k||^2/3,
leading to a characteristic sqrt(epsilon) correction structure in 3D.

We used the 6-point eps list:

epsilon in {1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6}.

---

## 4. Results

### 4.1 Raw fitted estimates (eps6, order=3)

All runs used:
- halfint_order = 3 (basis 1, epsilon^(1/2), epsilon, epsilon^(3/2))
- eps list above (6 points)
- grid transform power p in {2,3} (two independent quadrature families)

#### power = 2 (k = pi u^2)
N     Z_LGF estimate
200   1.3427206628084154
240   1.3431341507189185
280   1.3434399659789387

#### power = 3 (k = pi u^3)
N     Z_LGF estimate
200   1.3415057095094620
240   1.3420837457892207
280   1.3425122772663063

The monotone drift with N indicates a remaining grid (quadrature) systematic.

---

### 4.2 Grid-limit extrapolation N -> infinity

We use a simple leading-order model:

Z(N) ≈ Z_infinity + a/N.

Fitting separately for each power family (N = 200, 240, 280) yields:

- Z_infinity (power=2) ≈ 1.34523
- Z_infinity (power=3) ≈ 1.34502

A conservative combined estimate (grid/fit systematic dominated) is:

Z_LGF ≈ 1.3451 ± 0.0003

where the uncertainty reflects finite N, dependence on power p, simple 1/N model, and double precision.

---

## 5. Reproducibility (Windows CMD)

Environment: .venv on Python 3.11.9; no extra packages beyond requirements.txt needed.

### Single G(z) evaluation:

.\.venv\Scripts\python scripts\lgf\compute_zlgf_fcc_grid.py --z 0 --N 80 --power 2

(Sanity check: should return G(0) = 1.0 exactly.)

.\.venv\Scripts\python scripts\lgf\compute_zlgf_fcc_grid.py --z 0.5 --N 120 --power 2

(Example mid-point: G(0.5) ≈ 1.0268...)

### Full Z_LGF estimate (eps6, halfint order 3):

power = 2:
.\.venv\Scripts\python scripts\lgf\compute_zlgf_fcc_grid.py --N 280 --power 2 --halfint_order 3 ^
  --eps 1e-3 3e-4 1e-4 3e-5 1e-5 3e-6 --out_json results_lgf_N280_p2_eps6_h3.json

power = 3:
.\.venv\Scripts\python scripts\lgf\compute_zlgf_fcc_grid.py --N 280 --power 3 --halfint_order 3 ^
  --eps 1e-3 3e-4 1e-4 3e-5 1e-5 3e-6 --out_json results_lgf_N280_p3_eps6_h3.json

JSON output files are local only (ignored by git via .gitignore).

---

## 6. Limitations and caveats

- This is not a high-precision integral and not publication-level accuracy yet.
- Half-integer fit on powers epsilon^(m/2) can be ill-conditioned at large order; order=3 chosen as stable compromise.
- Difference between p=2 vs p=3 used as quadrature systematic estimate, but does not replace independent validation (literature / mpmath / specialized cubature).
- Extrapolation model Z(N)=Z_infinity+a/N is minimal; requires check against a/N+b/N^2, Richardson extrapolation, etc.

---

## 7. Next steps (for external review readiness)

To make this block peer-reviewable:

1. Convergence: extend table to N=320, 400; assess stable digits in Z_infinity.
2. Systematic model: compare Z_infinity under a/N vs a/N+b/N^2; make error budget.
3. Independent validation:
   - (A) high-precision (mpmath, possibly with domain splitting / analytic singularity subtraction), or
   - (B) literature (Joyce; Glasser-Zucker) with explicit normalization cross-check.
4. Documentation: short note in docs/ + reproducibility commands + tables.
5. ORT mapping: after stabilizing Z_LGF, tackle Z=F(Z_LGF) with no tunable parameters (ORT v9.0 Open Problem).

---

## 8. Embedding in ORT v9.0 canon

This computational block should be cited in the canon as:

- Computation / Evidence for FCC spectral-invariant bridge;
- not as "proof" of formula for Z, but as a step toward parameter-free mapping Z=F(Z_LGF);
- with explicit statement of achieved precision and sources of systematic uncertainty.

---

End of note.