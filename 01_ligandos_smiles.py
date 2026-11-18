#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import time

import pandas as pd
import requests

LIGANDOS = {
    "3NJO": ("A8S", "ABA"),
    "3QN1": ("A8S", "ABA"),
    "3KDI": ("A8S", "ABA"),
    "3OQU": ("A8S", "ABA"),
    "4JDA": ("A9S", "ABA (-)"),
    "3KB0": ("A8S", "ABA"),
    "3KDJ": ("A8S", "ABA"),
    "5MN0": ("A8S", "ABA"),
    "5UR5": ("8KP", "4m"),
    "4LGB": ("1W7", "AM1"),
    "3WG8": ("6AS", "AS6"),
    "6NWB": ("L6V", "hexabactin"),
    "6NWC": ("L6P", "3CB"),
    "3NJ0": ("PYV", "pyrabactin"),
    "5UR6": ("8KM", "cyanabactin"),
    "7MLC": ("ZLA", "4a"),
    "7MLD": ("ZLG", "antabactin"),
    "6ZUC": ("QPZ", "isb09"),
    "7AVW": ("S5N", "isb07"),
    "8AY7": ("S5N", "isb09"),
    "9H0H": ("A1IRJ", "opabactin"),
    "9H0I": ("A1IRN", "icb"),
}

CHEMCOMP_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{}"


def get_smiles(comp_id):
    resp = requests.get(CHEMCOMP_URL.format(comp_id), timeout=20)
    if resp.status_code != 200:
        print(f"   [HTTP {resp.status_code}] {comp_id}")
        return None, None
    data = resp.json()

    rcsb_desc = data.get("rcsb_chem_comp_descriptor", {}) or {}
    smiles = rcsb_desc.get("smiles_stereo") or rcsb_desc.get("smiles")

    if not smiles:
        descriptors = data.get("pdbx_chem_comp_descriptor", []) or []
        smiles_canon, smiles_openeye = None, None
        for d in descriptors:
            if d["type"] == "SMILES_CANONICAL":
                smiles_canon = d["descriptor"]
            elif d["type"] == "SMILES":
                smiles_openeye = d["descriptor"]
        smiles = smiles_canon or smiles_openeye

    if not smiles:
        print(f"   [SIN SMILES EN LA API] {comp_id} -> revisar a mano en "
              f"https://www.rcsb.org/ligand/{comp_id}")

    chem = data.get("chem_comp", {})
    return smiles, {"name": chem.get("name"), "formula": chem.get("formula")}


def main():
    filas = []
    for pdb_id, (comp_id, etiqueta) in LIGANDOS.items():
        smiles, meta = get_smiles(comp_id)
        filas.append({
            "pdb_id": pdb_id,
            "ligand_label": etiqueta,
            "ccd_code": comp_id,
            "smiles": smiles,
            "chem_name": meta["name"] if meta else None,
            "formula": meta["formula"] if meta else None,
        })
        time.sleep(0.3)   # cortesia con la API

    df = pd.DataFrame(filas)
    df.to_csv("ligandos_pyrpyl_smiles.csv", index=False)
    print(df)
    print("\nGuardado en ligandos_pyrpyl_smiles.csv")


if __name__ == "__main__":
    main()
