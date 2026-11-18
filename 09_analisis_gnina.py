#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import math
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

DRIVE = "/content/drive/MyDrive/tfm docking"
RES = os.path.join(DRIVE, "resultados")
OUT = os.path.join(DRIVE, "analisis")
MAN = "/content/ligandos/manifiesto.csv"
os.makedirs(OUT, exist_ok=True)

B_BOOT = 2000
B_PERM = 2000
ALPHA = 20.0
rng = np.random.default_rng(20261116)

man = pd.read_csv(MAN).drop_duplicates(subset="archivo")
clase = man.drop_duplicates(subset="id").set_index("id")["conjunto"]
clase.index = clase.index.astype(str)
print("moleculas en el manifiesto:", len(clase))
print(clase.value_counts().to_string())
print()


def auc_roc(y, s):
    y = np.asarray(y); s = np.asarray(s, float)
    n1 = int(y.sum()); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    r = pd.Series(s).rank().to_numpy()
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def enrich(y, s, frac):
    y = np.asarray(y); s = np.asarray(s, float)
    N = len(y); n = int(y.sum())
    if n == 0:
        return np.nan
    k = max(1, int(round(N * frac)))
    idx = np.argsort(-s, kind="mergesort")[:k]
    return (int(y[idx].sum()) / n) / (k / N)


def ef1(y, s):
    return enrich(y, s, 0.01)


def ef5(y, s):
    return enrich(y, s, 0.05)


def bedroc(y, s, alpha=ALPHA):
    y = np.asarray(y); s = np.asarray(s, float)
    N = len(y); n = int(y.sum())
    if n == 0 or n == N:
        return np.nan
    Ra = n / N
    orden = np.argsort(-s, kind="mergesort")
    rangos = np.where(y[orden] == 1)[0] + 1
    rie = np.sum(np.exp(-alpha * rangos / N)) / (Ra * (1 - math.exp(-alpha)) / (math.exp(alpha / N) - 1))
    f1 = Ra * math.sinh(alpha / 2) / (math.cosh(alpha / 2) - math.cosh(alpha / 2 - alpha * Ra))
    f2 = 1.0 / (1 - math.exp(alpha * (1 - Ra)))
    return rie * f1 + f2


METRICAS = [("AUC", auc_roc), ("EF1", ef1), ("EF5", ef5), ("BEDROC", bedroc)]


def boot_ic(y, s, fn, B=B_BOOT):
    y = np.asarray(y); s = np.asarray(s, float)
    ia = np.where(y == 1)[0]; idec = np.where(y == 0)[0]
    if len(ia) < 2 or len(idec) < 2:
        return np.nan, np.nan
    v = np.empty(B)
    for b in range(B):
        k = np.concatenate([rng.choice(ia, len(ia), replace=True),
                            rng.choice(idec, len(idec), replace=True)])
        v[b] = fn(y[k], s[k])
    return float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))


def perm_p(y, s, fn, B=B_PERM):
    y = np.asarray(y); s = np.asarray(s, float)
    obs = fn(y, s)
    if np.isnan(obs):
        return obs, np.nan
    nulos = np.empty(B)
    for b in range(B):
        nulos[b] = fn(rng.permutation(y), s)
    return float(obs), float((np.sum(nulos >= obs) + 1) / (B + 1))


series = {}

for f in sorted(glob.glob(os.path.join(RES, "gnina_*.csv"))):
    receptor = os.path.basename(f)[len("gnina_"):-4]
    d = pd.read_csv(f).drop_duplicates(subset="ligando", keep="first")
    d = d[d["estado"] == "ok"]
    d["ligando"] = d["ligando"].astype(str)
    d = d.set_index("ligando")
    series[("gnina-CNNaffinity", receptor)] = pd.to_numeric(d["cnn_affinity"], errors="coerce").dropna()
    series[("gnina-CNNaffinity-mejor", receptor)] = pd.to_numeric(d["mejor_cnn_affinity"], errors="coerce").dropna()
    series[("gnina-CNNscore", receptor)] = pd.to_numeric(d["cnn_score"], errors="coerce").dropna()
    series[("gnina-Vina", receptor)] = -pd.to_numeric(d["mejor_affinity"], errors="coerce").dropna()
    print("gnina", receptor, ":", len(d), "moleculas ok de", len(pd.read_csv(f)))

for f in sorted(glob.glob(os.path.join(RES, "*.csv"))):
    base = os.path.basename(f)[:-4]
    if base.startswith("gnina") or "_" not in base:
        continue
    programa, receptor = base.split("_", 1)
    d = pd.read_csv(f).drop_duplicates(subset="ligando", keep="first")
    d = d[d["estado"] == "ok"].copy()
    d["score"] = pd.to_numeric(d["score"], errors="coerce")
    d = d.dropna(subset=["score"]).merge(man, left_on="ligando", right_on="archivo")
    g = d.sort_values("score").groupby("id", as_index=False).first()
    s = pd.Series(-g["score"].to_numpy(), index=g["id"].astype(str))
    series[(programa, receptor)] = s

print()
print("series cargadas:", len(series))
print()


filas = []
for (metodo, receptor), s in sorted(series.items()):
    comunes = s.index.intersection(clase.index)
    s2 = s.loc[comunes]
    c2 = clase.loc[comunes]
    for grupo in ["TODOS", "core12"]:
        if grupo == "TODOS":
            y = (c2 != "decoy").astype(int).to_numpy()
            v = s2.to_numpy()
        else:
            m = (c2 == "decoy") | (c2 == grupo)
            y = (c2[m] == grupo).astype(int).to_numpy()
            v = s2[m].to_numpy()
        if y.sum() < 2:
            continue
        fila = {"metodo": metodo, "receptor": receptor, "activos": grupo,
                "n": len(y), "n_act": int(y.sum())}
        for nombre, fn in METRICAS:
            obs, p = perm_p(y, v, fn)
            lo, hi = boot_ic(y, v, fn)
            fila[nombre] = round(obs, 4)
            fila[nombre + "_inf"] = round(lo, 4) if lo == lo else np.nan
            fila[nombre + "_sup"] = round(hi, 4) if hi == hi else np.nan
            fila[nombre + "_p"] = round(p, 4)
        filas.append(fila)

met = pd.DataFrame(filas)
met.to_csv(os.path.join(OUT, "metricas_gnina.csv"), index=False)

print("=" * 110)
print("METRICAS  (TODOS = 15 activos)")
print("=" * 110)
t = met[met.activos == "TODOS"]
print(t[["metodo", "receptor", "n", "n_act", "AUC", "AUC_inf", "AUC_sup", "AUC_p",
         "EF1", "EF5", "BEDROC", "BEDROC_p"]].to_string(index=False))
print()


def comparar(sa, sb, B=B_BOOT):
    comunes = sa.index.intersection(sb.index).intersection(clase.index)
    if len(comunes) < 50:
        return None
    y = (clase.loc[comunes] != "decoy").astype(int).to_numpy()
    a = sa.loc[comunes].to_numpy(); b = sb.loc[comunes].to_numpy()
    if y.sum() < 2:
        return None
    obs = auc_roc(y, a) - auc_roc(y, b)
    ia = np.where(y == 1)[0]; idec = np.where(y == 0)[0]
    difs = np.empty(B)
    for i in range(B):
        k = np.concatenate([rng.choice(ia, len(ia), replace=True),
                            rng.choice(idec, len(idec), replace=True)])
        difs[i] = auc_roc(y[k], a[k]) - auc_roc(y[k], b[k])
    lo, hi = np.nanpercentile(difs, [2.5, 97.5])
    p = 2 * min((difs <= 0).mean(), (difs >= 0).mean())
    return dict(n=len(comunes), dAUC=round(float(obs), 4),
                inf=round(float(lo), 4), sup=round(float(hi), 4),
                p=round(float(min(1.0, p)), 4))


receptores_gnina = sorted({r for m, r in series if m.startswith("gnina")})
comp = []

for rec in receptores_gnina:
    ref = ("gnina-CNNaffinity", rec)
    if ref not in series:
        continue
    otros = [k for k in series if k[1] == rec and k != ref]
    for k in sorted(otros):
        r = comparar(series[ref], series[k])
        if r:
            r.update({"A": "gnina-CNNaffinity", "B": k[0], "receptor": rec})
            comp.append(r)

for m in sorted({k[0] for k in series}):
    if len(receptores_gnina) == 2 and (m, receptores_gnina[0]) in series and (m, receptores_gnina[1]) in series:
        r = comparar(series[(m, receptores_gnina[1])], series[(m, receptores_gnina[0])])
        if r:
            r.update({"A": m + " / " + receptores_gnina[1], "B": m + " / " + receptores_gnina[0],
                      "receptor": "entre receptores"})
            comp.append(r)

if comp:
    cmp_df = pd.DataFrame(comp)[["receptor", "A", "B", "n", "dAUC", "inf", "sup", "p"]]
    cmp_df.to_csv(os.path.join(OUT, "comparaciones_gnina.csv"), index=False)
    print("=" * 110)
    print("COMPARACIONES PAREADAS DE AUC  (dAUC = AUC(A) - AUC(B))")
    print("=" * 110)
    print(cmp_df.to_string(index=False))
    print()


def curva(y, s):
    o = np.argsort(-np.asarray(s, float), kind="mergesort")
    y = np.asarray(y)[o]
    return np.r_[0, np.cumsum(1 - y) / max(1, (1 - y).sum())], np.r_[0, np.cumsum(y) / max(1, y.sum())]


n = len(receptores_gnina)
if n:
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 5.0), squeeze=False)
    for ax, rec in zip(axes[0], receptores_gnina):
        for m in ["gnina-CNNaffinity", "gnina-CNNscore", "gnina-Vina", "smina"]:
            if (m, rec) not in series:
                continue
            s = series[(m, rec)]
            com = s.index.intersection(clase.index)
            y = (clase.loc[com] != "decoy").astype(int).to_numpy()
            fpr, tpr = curva(y, s.loc[com].to_numpy())
            ax.plot(fpr, tpr, lw=1.7, label="%s  AUC=%.3f" % (m, auc_roc(y, s.loc[com].to_numpy())))
        ax.plot([0, 1], [0, 1], "k--", lw=0.8, label="azar")
        ax.set_title(rec)
        ax.set_xlabel("Tasa de falsos positivos")
        ax.set_ylabel("Tasa de verdaderos positivos")
        ax.legend(loc="lower right", fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "roc_gnina.png"), dpi=200)
    plt.close()

print("Guardado en", OUT)
print("LISTO")
