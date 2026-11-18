#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regenera los ligandos en SDF desde los SMILES para GNINA"""

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

MAN = "/content/ligandos/manifiesto.csv"
OUT = "/content/ligandos/sdf"
os.makedirs(OUT, exist_ok=True)

man = pd.read_csv(MAN).drop_duplicates(subset="archivo")
unicas = man.drop_duplicates(subset="id")
print("moleculas unicas en el manifiesto:", len(unicas))

ok = 0
ya = 0
fallos = []

for _, r in unicas.iterrows():
    nombre = str(r["id"])
    dest = os.path.join(OUT, nombre + ".sdf")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        ya += 1
        continue
    smi = r.get("smiles")
    if not isinstance(smi, str) or not smi.strip():
        fallos.append((nombre, "sin smiles"))
        continue
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        fallos.append((nombre, "smiles invalido"))
        continue
    mol = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 42
    if AllChem.EmbedMolecule(mol, ps) != 0:
        ps.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, ps) != 0:
            fallos.append((nombre, "sin conformacion 3D"))
            continue
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    except Exception:
        pass
    mol.SetProp("_Name", nombre)
    w = Chem.SDWriter(dest)
    w.write(mol)
    w.close()
    ok += 1

print("SDF nuevos:", ok, "| ya existian:", ya, "| fallos:", len(fallos))
for f in fallos[:15]:
    print("   ", f[0], "->", f[1])
print("TOTAL EN", OUT, ":", len([f for f in os.listdir(OUT) if f.endswith(".sdf")]))
