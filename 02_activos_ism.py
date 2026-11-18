#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepara el fichero de activos que se sube a TLDR para generar los decoys."""

import pandas as pd

df = pd.read_csv("ligandos_pyrpyl_smiles.csv").dropna(subset=["smiles"])
unicas = df.drop_duplicates(subset=["smiles"]).copy()

with open("activos_dude_input.ism", "w") as f:
    for _, fila in unicas.iterrows():
        nombre = fila["ligand_label"].replace(" ", "_")
        f.write(f"{fila['smiles']} {nombre}\n")

print(f"Escritas {len(unicas)} moleculas unicas en activos_dude_input.ism")
print(unicas[["pdb_id", "ligand_label", "smiles"]])
