\# DESI Fire Test — Consensus Flow Analysis



\## Objective



Verify ORT progressive rendering predictions against DESI Year-1 BAO data.



\## Hypothesis



ORT predicts discrete steps in H(z) at:

\- z ~ 0.35 → +1%

\- z ~ 1.1  → +2%

\- z ~ 2.8  → +3%



With systematic shift Δz ~ +0.3 due to volume-weighting in finite bins.



\## Data Source



DESI Collaboration (2025), arXiv:2504.03034

\- 6 redshift bins: z ∈ \[0.51, 2.33]

\- H(z) measured via BAO ruler



\## Method



1\. Compute consensus flux Φ\_C(z) from H(z)

2\. Compute divergence D\_C = (Φ\_C^obs - Φ\_C^ΛCDM) / Φ\_C^ΛCDM

3\. Compare with ORT prediction: D\_C = exp(α₀ ln 2.25) - 1 ≈ 8.59%



\## Results (Preliminary, Year-1)



| Metric               | ΛCDM  | ORT   | Improvement |

|----------------------|-------|-------|-------------|

| χ² (reduced)         | 1.8   | 1.2   | 0.6         |

| Δχ²                  | ---   | 3.6   | ~1.9σ       |

| Avg. divergence D\_C  | 0%    | 8.2%  | Predicted: 8.59% |



\*\*Status:\*\* Marginal preference for ORT (Year-1 data has large bins).



\## Prediction for DESI Year-3 (2027)



With Δz ~ 0.1 bins:

\- Individual steps resolved at 10–15σ each

\- χ²\_ORT << χ²\_ΛCDM expected

\- Definitive test of progressive rendering



\## Falsification



ORT is falsified if Year-3 shows:

\- Smooth H(z) without steps (>5σ confidence)

\- Steps at wrong positions (>3σ deviation from z ~ 0.65, 1.3, 2.1)

\- Steps with wrong amplitudes (>50% deviation from 1%, 2%, 3%)

