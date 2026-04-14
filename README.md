# Protein-structural-alignment
Implementation of τ-equivalence algorithm from Subbarao &amp; Haneef (1991) for protein structural alignment
# Protein Structural Alignment using τ-Equivalence

Implementation of the graph-theoretical structural alignment method 
from Subbarao & Haneef (1991) - "Defining Topological Equivalences 
in Macromolecules"

## What this does
Finds topologically equivalent residues between two protein structures
using interatomic distance matching — without requiring prior 
superposition.

## Protein pairs studied
- L. casei DHFR vs E. coli DHFR (closely related)
- Haemoglobin-α vs Myoglobin (homologous)
- T4 Lysozyme vs Hen Egg White Lysozyme (distantly related)
- Cytochrome-c vs Cytochrome-c2 (homologous)
- Azurin vs Plastocyanin (distantly related)

## Requirements
pip install biopython numpy scipy matplotlib pandas

## How to run
python 01_download.py    # downloads PDB files
python 02_tau_equivalence.py   # runs the algorithm
python 03_plots.py       # generates figures

## Key Results
| Protein Pair | δ-equivalences | RMSD (Å) |
|---|---|---|
| DHFR | 137 | 1.09 |
| Globins | 115 | 1.19 |
| Lysozymes | 23 | 2.17 |
| Cytochromes | 77 | 0.85 |
| Azurin/Plastocyanin | 43 | 1.68 |

## Reference
Subbarao, N. and Haneef, I. (1991) Defining topological equivalences
in macromolecules. Protein Engineering, 4(8), 877–884.

## Tools used
Python, BioPython, PyMOL, VS Code
