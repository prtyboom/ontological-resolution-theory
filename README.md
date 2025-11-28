\# Ontological Resolution Theory (ORT)



\*\*Gravity as memory, dark matter as archive, cognition as indexing\*\*



This repository contains code, data, and LaTeX sources for \*\*Ontological Resolution Theory (ORT)\*\* — an ontological framework where:



\- \*\*gravity\*\* is interpreted as the distributed memory of irreversible events (distinctions),

\- \*\*dark matter\*\* is an \*\*archival sector\*\* — highly scrambled information that no longer participates in Standard Model interactions but still gravitates,

\- \*\*consciousness\*\* is a \*\*cognitive operator\*\* that indexes this physical log into structured histories and meaning.



The project currently has three main layers:



1\. \*\*Emergent spacetime and gravity\*\* (consensus field, holographic emergence, GR as rendering) — see the legacy `ORT/` subproject.

2\. \*\*Dark matter as archival sector\*\* (galaxy rotation curves, dwarfs, entropy-halo toy model) — in `applications/dark-matter/`.

3\. \*\*Ontological closure\*\* (gravity as memory, cognition as indexing) — in `papers/ORT\_closing\_paper.tex`.



---



\## Repository structure



```text

ontological-resolution-theory/

├── ORT/                       # Original ORT v2.0 project (quantum uncertainty \& holography)

│   ├── v2.0/                  # PDF + LaTeX for ORT v2.0

│   └── README.md              # Detailed description of that project

├── applications/

│   └── dark-matter/           # Dark matter as archival sector

│       ├── data/              # CSV data: SPARC subsets, dwarfs, etc.

│       ├── calculations/      # Scripts: NFW fits, MOND comparison, dwarfs

│       ├── simulations/       # Toy entropy-halo simulation

│       └── results/           # Generated figures and JSON summaries (gitignored)

├── papers/

│   ├── ORT\_closing\_paper.tex  # Ontological closure: gravity as memory, cognition as indexing

│   └── DM\_holographic\_paper.tex (planned) # Dark Matter from Holographic Saturation

├── requirements.txt           # Python dependencies for dark-matter applications

└── README.md                  # This file

Dark matter application

All code related to the archival interpretation of dark matter lives in:



text



applications/dark-matter/

Environment

Install Python dependencies:



Bash



pip install -r requirements.txt

Minimal stack:



numpy

scipy

matplotlib

pandas

astropy

requests

Reproducing key results

From the project root:



Bash



cd applications/dark-matter

Synthetic SPARC-like test:



Bash



cd calculations

python 04\_sparc\_analysis.py

NGC 2403: real SPARC data, NFW fit



Bash



python 05\_sparc\_real\_fit.py

Produces:



Best-fit NFW parameters (rho0, rs)

Reduced chi-squared χ²/dof ≈ 0.52

Plot in ../results/figures/ngc2403\_nfw\_fit.png

NGC 3198: NFW fit and ORT vs MOND



Bash



python 06\_sparc\_multi\_fit.py       # NGC 3198 NFW fit

python 07\_ort\_vs\_nfw\_final.py      # ORT entropy-halo vs NFW comparison

python 08\_mond\_comparison.py       # NFW vs MOND χ² comparison

Entropy-halo toy simulation (DM/baryon ≈ 5.36)



Bash



cd ../simulations

python entropy\_halo\_sim.py

Dwarf galaxies: ultra-faint M/L



Bash



cd ../calculations

python 09\_dwarf\_analysis.py

Generated plots live in applications/dark-matter/results/figures/ (ignored by git).



Papers

LaTeX sources for conceptual and technical writeups live in papers/:



ORT\_closing\_paper.tex

Ontological Resolution Theory: Gravity as Memory, Cognition as Indexing



DM\_holographic\_paper.tex (to be added)

Dark Matter from Holographic Saturation: A Conceptual Framework with Preliminary Estimates and Galaxy-Scale Tests



The legacy ORT v2.0 paper on information-theoretic derivation of uncertainty and holographic bounds is in ORT/v2.0/ and documented in ORT/README.md.



Zenodo records / Публикации

Dark Matter from Holographic Saturation: A Conceptual Framework with Preliminary Estimates and Galaxy-Scale Tests

(RU: Тёмная материя из голографического насыщения: концептуальная основа с предварительными оценками)

Zenodo: https://zenodo.org/records/17721415

DOI: 10.5281/zenodo.17721415



Ontological Resolution Theory: Gravity as Memory, Cognition as Indexing

(RU: Теория онтологического разрешения: гравитация как память, познание как индексация)

Zenodo: https://zenodo.org/records/17741356

DOI: 10.5281/zenodo.17741356



Status

This is a conceptual and exploratory project with preliminary numerical tests.

It is not a finished replacement for ΛCDM or standard quantum gravity.



The goals are to:



explore an archival interpretation of dark matter compatible with galaxy-scale data,

reinterpret gravity as distributed memory of distinctions,

model consciousness as a cognitive indexing layer over this physical log,

and formulate concrete research directions for further development and testing.

Author

Fedor Kapitanov

Independent Researcher, Moscow

ORCID: 0009-0009-6438-8730

Contact: prtyboom@gmail.com

