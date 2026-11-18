#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, glob, math, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")

DRIVE = "/content/drive/MyDrive/tfm docking"
RES = os.path.join(DRIVE, "resultados")
OUT = os.path.join(DRIVE, "analisis")
MAN = "/content/ligandos/manifiesto.csv"
os.makedirs(OUT, exist_ok=True)

B_BOOT = 2000      # remuestreos bootstrap
B_PERM = 2000      # permutaciones para el p-valor
ALPHA_BEDROC = 20.0
rng = np.random.default_rng(20261114)

# 1. Colapso de protomeros

man = pd.read_csv(MAN).drop_duplicates(subset="archivo")
print("MANIFIESTO:", man.shape[0], "ficheros |", man["id"].nunique(), "moleculas unicas")
print(man["conjunto"].value_counts().to_string())
print()

CLASES = [c for c in man["conjunto"].unique() if c != "decoy"]


def cargar(path):
    """Devuelve un DataFrame con una fila por molecula (mejor protomero)."""
    d = pd.read_csv(path)
    n_bruto = len(d)
    d = d.drop_duplicates(subset="ligando", keep="first")
    n_fallos = int((d["estado"] != "ok").sum())
    d = d[d["estado"] == "ok"].copy()
    d["score"] = pd.to_numeric(d["score"], errors="coerce")
    d = d.dropna(subset=["score"])
    d = d.merge(man, left_on="ligando", right_on="archivo", how="inner")
    # colapso de protomeros: se queda el mas negativo (mejor) por molécula
    d = d.sort_values("score").groupby("id", as_index=False).first()
    d["mejor"] = -d["score"]
    return d, n_bruto, n_fallos


# 2. Metricas

def auc_roc(y, s):
    y = np.asarray(y); s = np.asarray(s, float)
    n1 = int(y.sum()); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    r = pd.Series(s).rank().to_numpy()
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def enrich(y, s, frac):
    """Factor de enriquecimiento al frac (0.01 = 1%)."""
    y = np.asarray(y); s = np.asarray(s, float)
    N = len(y); n = int(y.sum())
    if n == 0:
        return np.nan
    k = max(1, int(round(N * frac)))
    idx = np.argsort(-s, kind="mergesort")[:k]
    hits = int(y[idx].sum())
    return (hits / n) / (k / N)


def ef1(y, s):
    return enrich(y, s, 0.01)


def ef5(y, s):
    return enrich(y, s, 0.05)


def bedroc(y, s, alpha=ALPHA_BEDROC):
    """BEDROC de Truchon & Bayly (2007). s: mayor = mejor."""
    y = np.asarray(y); s = np.asarray(s, float)
    N = len(y); n = int(y.sum())
    if n == 0 or n == N:
        return np.nan
    Ra = n / N
    orden = np.argsort(-s, kind="mergesort")
    rangos = np.where(y[orden] == 1)[0] + 1
    rie_num = np.sum(np.exp(-alpha * rangos / N))
    rie_den = Ra * (1 - math.exp(-alpha)) / (math.exp(alpha / N) - 1)
    rie = rie_num / rie_den
    f1 = Ra * math.sinh(alpha / 2) / (math.cosh(alpha / 2) - math.cosh(alpha / 2 - alpha * Ra))
    f2 = 1.0 / (1 - math.exp(alpha * (1 - Ra)))
    return rie * f1 + f2


METRICAS = [("AUC", auc_roc), ("EF1", ef1), ("EF5", ef5), ("BEDROC", bedroc)]


def boot_ic(y, s, fn, B=B_BOOT):
    """IC95% percentil, remuestreando activos y decoys por separado."""
    y = np.asarray(y); s = np.asarray(s, float)
    ia = np.where(y == 1)[0]; idec = np.where(y == 0)[0]
    if len(ia) < 2 or len(idec) < 2:
        return np.nan, np.nan
    vals = np.empty(B)
    for b in range(B):
        k = np.concatenate([rng.choice(ia, len(ia), replace=True),
                            rng.choice(idec, len(idec), replace=True)])
        vals[b] = fn(y[k], s[k])
    return float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))


def perm_p(y, s, fn, B=B_PERM):
    """p-valor de una cola por permutacion de etiquetas + media de la nula."""
    y = np.asarray(y); s = np.asarray(s, float)
    obs = fn(y, s)
    if np.isnan(obs):
        return obs, np.nan, np.nan
    nulos = np.empty(B)
    for b in range(B):
        nulos[b] = fn(rng.permutation(y), s)
    p = (np.sum(nulos >= obs) + 1) / (B + 1)
    return float(obs), float(p), float(np.nanmean(nulos))

ficheros = sorted(glob.glob(os.path.join(RES, "*.csv")))
print("Ficheros de resultados encontrados:", len(ficheros))

filas = []
datos = {}          # (programa, receptor) -> DataFrame colapsado

for f in ficheros:
    base = os.path.basename(f)[:-4]
    if "_" not in base:
        continue
    programa, receptor = base.split("_", 1)
    d, n_bruto, n_fallos = cargar(f)
    if len(d) == 0:
        continue
    datos[(programa, receptor)] = d
    cobertura = man["archivo"].isin(d["ligando"]).mean()

    for clase in ["TODOS"] + CLASES:
        if clase == "TODOS":
            sub = d.copy()
            sub["y"] = (sub["conjunto"] != "decoy").astype(int)
        else:
            sub = d[(d["conjunto"] == "decoy") | (d["conjunto"] == clase)].copy()
            sub["y"] = (sub["conjunto"] == clase).astype(int)
        y = sub["y"].to_numpy(); s = sub["mejor"].to_numpy()
        n_act = int(y.sum()); n_dec = len(y) - n_act
        if n_act < 2 or n_dec < 10:
            continue

        fila = {"programa": programa, "receptor": receptor, "activos_usados": clase,
                "n_moleculas": len(y), "n_activos": n_act, "n_decoys": n_dec,
                "ligandos_evaluados": n_bruto, "fallos_docking": n_fallos,
                "cobertura_fichero": round(float(cobertura), 3)}

        for nombre, fn in METRICAS:
            obs, p, nulo = perm_p(y, s, fn)
            lo, hi = boot_ic(y, s, fn)
            fila[nombre] = round(obs, 4)
            fila[nombre + "_IC95_inf"] = round(lo, 4) if lo == lo else np.nan
            fila[nombre + "_IC95_sup"] = round(hi, 4) if hi == hi else np.nan
            fila[nombre + "_p_perm"] = round(p, 4)
            fila[nombre + "_nula"] = round(nulo, 4)

        u, p_mw = mannwhitneyu(s[y == 1], s[y == 0], alternative="greater")
        fila["MannWhitney_p"] = float(p_mw)
        fila["score_medio_activos"] = round(float(-s[y == 1].mean()), 3)
        fila["score_medio_decoys"] = round(float(-s[y == 0].mean()), 3)
        fila["mejor_score"] = round(float(-s.max()), 2)
        filas.append(fila)

met = pd.DataFrame(filas)
met.to_csv(os.path.join(OUT, "metricas.csv"), index=False)

print()
print("=" * 100)
print("METRICAS DE CRIBADO VIRTUAL")
print("=" * 100)
cols = ["programa", "receptor", "activos_usados", "n_activos", "n_decoys",
        "AUC", "AUC_IC95_inf", "AUC_IC95_sup", "AUC_p_perm",
        "EF1", "EF5", "BEDROC", "BEDROC_p_perm"]
if len(met):
    print(met[cols].to_string(index=False))
print()


def auc_pareado(dA, dB, B=B_BOOT):
    """IC95% y p bilateral de AUC(A) - AUC(B) sobre las moleculas comunes."""
    j = dA[["id", "conjunto", "mejor"]].merge(
        dB[["id", "mejor"]], on="id", suffixes=("_A", "_B"))
    y = (j["conjunto"] != "decoy").astype(int).to_numpy()
    sA = j["mejor_A"].to_numpy(); sB = j["mejor_B"].to_numpy()
    if y.sum() < 2 or (len(y) - y.sum()) < 10:
        return None
    obs = auc_roc(y, sA) - auc_roc(y, sB)
    ia = np.where(y == 1)[0]; idec = np.where(y == 0)[0]
    difs = np.empty(B)
    for b in range(B):
        k = np.concatenate([rng.choice(ia, len(ia), replace=True),
                            rng.choice(idec, len(idec), replace=True)])
        difs[b] = auc_roc(y[k], sA[k]) - auc_roc(y[k], sB[k])
    lo, hi = np.nanpercentile(difs, [2.5, 97.5])
    p = 2 * min((difs <= 0).mean(), (difs >= 0).mean())
    return dict(n_comun=len(y), dAUC=round(float(obs), 4),
                IC95_inf=round(float(lo), 4), IC95_sup=round(float(hi), 4),
                p=round(float(min(1.0, p)), 4))


comp = []
claves = sorted(datos.keys())
for i in range(len(claves)):
    for j in range(i + 1, len(claves)):
        (pA, rA), (pB, rB) = claves[i], claves[j]
        if pA != pB and rA != rB:
            continue          # solo comparaciones que cambian una cosa
        r = auc_pareado(datos[claves[i]], datos[claves[j]])
        if r:
            r.update({"A": f"{pA}_{rA}", "B": f"{pB}_{rB}",
                      "tipo": "receptores" if pA == pB else "programas"})
            comp.append(r)

if comp:
    cmp_df = pd.DataFrame(comp)[["tipo", "A", "B", "n_comun", "dAUC", "IC95_inf", "IC95_sup", "p"]]
    cmp_df.to_csv(os.path.join(OUT, "comparaciones.csv"), index=False)
    print("=" * 100)
    print("COMPARACIONES PAREADAS DE AUC  (dAUC = AUC(A) - AUC(B))")
    print("=" * 100)
    print(cmp_df.to_string(index=False))
    print()


for (prog, rec), d in datos.items():
    t = d.sort_values("score")[["id", "conjunto", "activo_padre", "score", "smiles"]].head(20)
    t.insert(0, "puesto", range(1, len(t) + 1))
    t.to_csv(os.path.join(OUT, f"top20_{prog}_{rec}.csv"), index=False)

print("Rankings top-20 guardados:", len(datos), "ficheros")

def curva_roc(y, s):
    o = np.argsort(-np.asarray(s, float), kind="mergesort")
    y = np.asarray(y)[o]
    tpr = np.cumsum(y) / max(1, y.sum())
    fpr = np.cumsum(1 - y) / max(1, (1 - y).sum())
    return np.r_[0, fpr], np.r_[0, tpr]


programas = sorted({p for p, _ in datos})
for prog in programas:
    recs = sorted([r for p, r in datos if p == prog])
    if not recs:
        continue
    plt.figure(figsize=(5.5, 5.5))
    for rec in recs:
        d = datos[(prog, rec)]
        y = (d["conjunto"] != "decoy").astype(int).to_numpy()
        s = d["mejor"].to_numpy()
        fpr, tpr = curva_roc(y, s)
        plt.plot(fpr, tpr, lw=1.8, label=f"{rec}  AUC={auc_roc(y, s):.3f}")
    plt.plot([0, 1], [0, 1], "k--", lw=0.8, label="azar")
    plt.xlabel("Tasa de falsos positivos")
    plt.ylabel("Tasa de verdaderos positivos")
    plt.title(f"Curvas ROC - {prog}")
    plt.legend(loc="lower right", fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, f"roc_{prog}.png"), dpi=200)
    plt.close()

n = len(datos)
if n:
    fig, axes = plt.subplots(1, n, figsize=(3.1 * n, 4.2), squeeze=False)
    for ax, k in zip(axes[0], claves):
        d = datos[k]
        act = d.loc[d["conjunto"] != "decoy", "score"].to_numpy()
        dec = d.loc[d["conjunto"] == "decoy", "score"].to_numpy()
        ax.boxplot([dec, act], widths=0.6, showfliers=False)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["decoys", "activos"])
        ax.scatter(np.full(len(act), 2) + rng.normal(0, 0.05, len(act)), act,
                   s=14, alpha=0.8, zorder=3)
        ax.set_title(f"{k[0]}\n{k[1]}", fontsize=9)
        ax.set_ylabel("score (kcal/mol)")
        ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "distribucion_scores.png"), dpi=200)
    plt.close()

print("Figuras guardadas en", OUT)
print("LISTO")
