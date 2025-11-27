\# Dark Matter from Holographic Saturation



\[!\[DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17721415.svg)](https://doi.org/10.5281/zenodo.17721415)



Application of ORT framework to the dark matter problem.



\## Paper



\*\*Dark Matter from Holographic Saturation: A Conceptual Framework with Preliminary Estimates\*\*



\- PDF: \[Zenodo](https://doi.org/10.5281/zenodo.17721415)

\- Status: Conceptual proposal with preliminary estimates



\## Key Idea



When comoving regions approach the Bekenstein-Hawking entropy bound during the radiation epoch, a fraction of degrees of freedom decouple from Standard Model gauge interactions while remaining gravitationally coupled.



\## Results Summary



| Quantity | Value | Notes |

|----------|-------|-------|

| Saturation temperature | T\_sat ~ 10^18 GeV | From holographic bound |

| Decoupling temperature | T\_dec ~ 10^16 GeV | Input (reheating scale) |

| Required g\_total | ~516 | To match Ω\_DM/Ω\_b = 5.36 |

| Decoupling fraction | ~79% | g\_arch / g\_total |

| Compatible BSM | E6 GUT, String | g\_\* ~ 500-700 |

| Late-time contribution | < 1% | DM is primordial |



\## Calculations



| Script | Description | Status |

|--------|-------------|--------|

| `01\\\_saturation\\\_temperature.py` | T\_sat from holographic bound | ✅ Done |

| `02\\\_dof\\\_evolution.py` | g\_\*(T) for Standard Model | ✅ Done |

| `03\\\_dm\\\_baryon\\\_ratio.py` | Required g\_total for Ω ratio | ✅ Done |

| `04\\\_sparc\\\_analysis.py` | NFW rotation curve fit | ✅ Done |

| `05\\\_entropy\\\_budget.py` | Entropy budget analysis | ✅ Done |



\## Run Calculations



```bash

cd calculations

pip install numpy scipy matplotlib

python run\\\_all.py

Figures

results/figures/g\\\_star\\\_evolution.png - g\\\_\\\*(T) evolution

results/figures/dm\\\_baryon\\\_ratio.png - Ω\\\_DM/Ω\\\_b vs g\\\_total

results/figures/rotation\\\_curve.png - NFW fit example

Next Steps

Simulation: "Entropy logs → halo formation" N-body code

JWST data: Test against "Little Red Dots" observations

Full SPARC fit: Fit all 175 galaxies

Citation

bibtex



@misc{kapitanov2025darkmatter,

\&nbsp; author = {Kapitanov, Fedor},

\&nbsp; title = {Dark Matter from Holographic Saturation: A Conceptual Framework with Preliminary Estimates},

\&nbsp; year = {2025},

\&nbsp; publisher = {Zenodo},

\&nbsp; doi = {10.5281/zenodo.17721415},

\&nbsp; url = {https://doi.org/10.5281/zenodo.17721415}

}

License

MIT License




