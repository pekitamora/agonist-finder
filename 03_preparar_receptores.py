import numpy as np

REF, LIG, CH = "5MN0.pdb", "A8S", "A"
JUNK = {"HOH","WAT","GOL","EDO","PEG","MPD","DMS","SO4","PO4","CL","NA","K",
        "BR","IOD","ACT","MN","MG","ZN","CA","FE","CU","NI","CO","CD"}

def rd(p):
    return [l.rstrip("\n") for l in open(p)
            if l[:6] in ("ATOM  ","HETATM") and len(l) >= 54 and l[16] in " A"]

def xyz(l):
    return np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])

def cas(ls, ch):
    return {int(l[22:26]): xyz(l) for l in ls
            if l[:6] == "ATOM  " and l[21] == ch and l[12:16].strip() == "CA"}

def kabsch(P, Q):
    cP, cQ = P.mean(0), Q.mean(0)
    P0, Q0 = P - cP, Q - cQ
    V, _, Wt = np.linalg.svd(P0.T @ Q0)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    return R, cQ - cP @ R, float(np.sqrt((((P0 @ R) - Q0) ** 2).sum() / len(P)))

def mv(l, R, t):
    x, y, z = xyz(l) @ R + t
    return f"{l[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{l[54:]}"

def wr(p, ls):
    open(p, "w").write("\n".join(ls) + "\nEND\n")

ref = rd(REF)
cref = cas(ref, CH)
for src, out in [("5MMQ.pdb","5MMQ_alineado.pdb"), ("5MMX.pdb","5MMX_alineado.pdb")]:
    mob = rd(src)
    cmob = cas(mob, CH)
    k = sorted(set(cref) & set(cmob))
    R, t, r = kabsch(np.array([cmob[i] for i in k]), np.array([cref[i] for i in k]))
    wr(out, [mv(l, R, t) for l in mob])
    print(f"{src}: {len(k)} residuos emparejados, RMSD {r:.3f} A -> {out}")

lig = [l for l in ref if l[17:20].strip() == LIG]
claves = sorted({l[21] + l[22:27] for l in lig})
lig = [l for l in lig if l[21] + l[22:27] == claves[0]]
wr("ligando_ref.pdb", lig)
c = np.array([xyz(l) for l in lig])
cen, ext = c.mean(0), c.max(0) - c.min(0)
caja = np.maximum(25.0, np.round(ext + 20.0, 1))
print(f"\nligando {LIG}: {len(lig)} atomos, copias presentes {claves}")
print("elementos:", sorted({l[76:78].strip() for l in lig}))
print(f"CENTRO {cen[0]:.3f} {cen[1]:.3f} {cen[2]:.3f}")
print(f"CAJA   {caja[0]} {caja[1]} {caja[2]}\n")

def prot(ls, chs):
    return [l for l in ls if l[:6] == "ATOM  " and l[21] in chs
            and l[17:20].strip() not in JUNK]

wr("tern_conPP2C_receptor.pdb", prot(ref, "AB"))
wr("tern_sinPP2C_receptor.pdb", prot(ref, "A"))
wr("holo_receptor.pdb", prot(rd("5MMX_alineado.pdb"), "A"))
wr("apo_receptor.pdb", prot(rd("5MMQ_alineado.pdb"), "A"))
for n in ["tern_conPP2C","tern_sinPP2C","holo","apo"]:
    print(n, sum(1 for _ in open(n + "_receptor.pdb")) - 1, "atomos")
