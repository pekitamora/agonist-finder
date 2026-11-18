import csv, io, os, re, sys, tarfile, zipfile

Z = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else "ligandos"
P = re.compile(r"^(.+?)\.(\d+\.[NM])\.pdbqt$")
d = OUT + "/pdbqt"
os.makedirs(d, exist_ok=True)
z = zipfile.ZipFile(Z)

act = {}
for l in z.read("decoys/ligands.smi").decode().splitlines():
    p = l.split()
    if len(p) >= 2:
        act[p[1]] = p[0]

dec = {}
for n in z.namelist():
    if "/maps/" in n and n.endswith(".txt"):
        pa = None
        for l in z.read(n).decode().splitlines():
            l = l.strip()
            if l.startswith("LIGAND:"):
                pa = l.split()[2].split("_")[0]
            elif l.startswith("DECOY"):
                c = l.split()
                dec[c[3]] = (c[2], pa)

rows = []
tg = [n for n in z.namelist() if n.endswith(".tar.gz") and "build3d" in n][0]
T = tarfile.open(fileobj=io.BytesIO(z.read(tg)), mode="r:gz")
for m in T.getmembers():
    if not m.name.endswith(".tar.gz"):
        continue
    try:
        T2 = tarfile.open(fileobj=io.BytesIO(T.extractfile(m).read()), mode="r:gz")
        ms = [x for x in T2.getmembers() if x.isfile()]
    except Exception:
        ms = []
    for x in ms:
        mm = P.match(x.name.split("/")[-1])
        if not mm:
            continue
        mid, pr = mm.group(1), mm.group(2)
        fn = mid + "__" + pr.replace(".", "_") + ".pdbqt"
        open(d + "/" + fn, "wb").write(T2.extractfile(x).read())
        s, pa = dec.get(mid, ("", "?"))
        rows.append([fn, mid, pr, "decoy", pa, s])
print("decoys pdbqt:", len(rows), "| moleculas unicas:", len({r[1] for r in rows}))

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    ok = 0
    for i, smi in sorted(act.items(), key=lambda kv: int(kv[0])):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print("ERROR smiles activo", i)
            continue
        mol = Chem.AddHs(mol)
        ps = AllChem.ETKDGv3()
        ps.randomSeed = 42
        if AllChem.EmbedMolecule(mol, ps) != 0:
            print("ERROR 3D activo", i)
            continue
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
        for cfg in MoleculePreparation().prepare(mol):
            t, good, msg = PDBQTWriterLegacy.write_string(cfg)
            if not good:
                print("ERROR meeko activo", i, msg)
                break
            fn = "activo_" + i + "__0_N.pdbqt"
            open(d + "/" + fn, "w").write(t)
            rows.append([fn, "activo_" + i, "0.N", "activo", i, smi])
            ok += 1
            break
    print("activos construidos:", ok, "de", len(act))
except ImportError as e:
    print("FALTA rdkit o meeko:", e)

w = csv.writer(open(OUT + "/manifiesto.csv", "w", newline=""))
w.writerow(["archivo", "id", "protomero", "tipo", "activo_padre", "smiles"])
w.writerows(rows)
print("total pdbqt:", len(rows), "->", OUT)
