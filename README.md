# agonist-finder
ABA-AgonistFinder: A structure-based tool for the discovery of novel ABA receptor agonists
ABA receptors (PYR/PYL/RCAR family) only close their gate and latch once an agonist is bound, and only then do they recruit a PP2C phosphatase. So when you want to screen for agonists computationally, you have to pick a structure to dock against. I couldn't find a straight answer in the literature, so I ran the comparison. Four receptor structures, five scoring functions, 15 known agonists and 660 property-matched decoys. This repo is the code for that.


RESULTS

In general, the ternary complex gave the best results, and it's not because it separates actives better overall, it's because it puts them higher up.

Receptor qvina2 smina vinardo gnina-Vina gnina-CNN
apo (5MMQ) 0.907 / 0.426 0.921 / 0.443 0.774 / 0.455 0.910 / 0.393 0.833 / 0.198
holo (5MMX) 0.791 / 0.173 0.796 / 0.110 0.746 / 0.172 — —
ternary, no PP2C (5MN0) 0.894 / 0.510 0.925 / 0.572 0.844 / 0.326 — —
ternary, full (5MN0) 0.912 / 0.555 0.932 / 0.576 0.893 / 0.392 0.920 / 0.483 0.779 / 0.381

AUC / BEDROC(α=20).


PROTEINS

All four share a coordinate frame, 5MMQ and 5MMX are superposed onto 5MN0 by Kabsch on Cα atoms with a 10 Å margin:

center  2.321   6.301  32.385
size    25.7 × 26.3 × 29.0 Å

Waters were removed everywhere. Protonation and Gasteiger charges from prepare_receptor4.py -A checkhydrogens -U nphs_lps_waters.


LIGANDS

The actives come from the PDB's Chemical Component Dictionary, not PubChem or ChEMBL, because most of them aren't in either. 22 PDB entries collapsed to 15 unique molecules. Decoys were DUDE-Z generated through TLDR. They arrived as pdbqt with several protomers per molecule. The analysis collapses protomers and keeps the best score per molecule.


DOCKING SETUP

qvina2 2.1, smina with and without --scoring vinardo, and GNINA 1.3.3 for the CNN rescoring. Everything at --exhaustiveness 8, --num_modes 9, --seed 42.
GNINA only ran on apo and the full ternary, and with the ligands rebuilt as SDF from SMILES.

METRICS

AUC via the Mann-Whitney U equivalence. EF at 1% and 5%. BEDROC with α = 20, which puts 80% of the weight in the top 8% of the ranking.

RUNNING IT

#bash
pip install -r requirements.txt
cd scripts
python 01_ligandos_smiles.py
python 02_activos_ism.py
python 03_preparar_receptores.py
python 04_preparar_ligandos.py decoys.zip ligandos/
python 05_correr_docking.py EXE REC.pdbqt ligandos/pdbqt poses/ out.csv NJ 8
python 08_analisis.py

#only for GNINA
python 06_ligandos_sdf.py
python 07_correr_gnina.py REC.pdbqt ligandos/sdf poses/ out.csv NJ 8
python 09_analisis_gnina.py

Both scritps are resumable. They can skip what's already on the CSV.


LAYOUT

scripts/
  01_ligandos_smiles.py      actives' SMILES from the RCSB API
  02_activos_ism.py          input file for DUDE-Z decoy generation
  03_preparar_receptores.py  Kabsch superposition, box, four receptors
  04_preparar_ligandos.py    unpacks decoys, builds the manifest
  05_correr_docking.py       parallel resumable docking (qvina2 / smina)
  06_ligandos_sdf.py         ligands as SDF for GNINA
  07_correr_gnina.py         GNINA runs, also resumable
  08_analisis.py             metrics, bootstrap, permutations, figures
  09_analisis_gnina.py       same for the GNINA score variants
notebooks/
