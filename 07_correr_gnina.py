#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

REC, LIGD, OUTD, CSVP = sys.argv[1:5]
NJ = int(sys.argv[5]) if len(sys.argv) > 5 else 4
EXH = sys.argv[6] if len(sys.argv) > 6 else "8"

GNINA = "/content/gnina"
REF = "/content/ligando_ref.pdb"
CABECERA = ["ligando", "estado", "affinity", "cnn_score", "cnn_affinity",
            "mejor_affinity", "mejor_cnn_affinity", "segundos"]

os.makedirs(OUTD, exist_ok=True)
os.makedirs(os.path.dirname(CSVP) or ".", exist_ok=True)

done = set()
if os.path.exists(CSVP):
    done = {r[0] for r in csv.reader(open(CSVP)) if r and r[0] != "ligando"}

todos = sorted(f[:-4] for f in os.listdir(LIGD) if f.endswith(".sdf"))
ligs = [n for n in todos if n not in done]
print("ligandos %d | ya hechos %d | pendientes %d" % (len(todos), len(done), len(ligs)),
      flush=True)

FILA = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")


def run(n):
    t0 = time.time()
    cmd = [GNINA, "-r", REC, "-l", os.path.join(LIGD, n + ".sdf"),
           "--autobox_ligand", REF, "--autobox_add", "10",
           "--exhaustiveness", EXH, "--num_modes", "9",
           "--seed", "42", "--cpu", "1",
           "-o", os.path.join(OUTD, n + ".sdf")]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    except Exception:
        return [n, "timeout", "", "", "", "", "", round(time.time() - t0, 1)]
    if p.returncode != 0:
        return [n, "error", "", "", "", "", "", round(time.time() - t0, 1)]

    modos = []
    for linea in p.stdout.splitlines():
        m = FILA.match(linea)
        if m:
            modos.append((int(m.group(1)), float(m.group(2)),
                          float(m.group(4)), float(m.group(5))))
    if not modos:
        return [n, "sin_score", "", "", "", "", "", round(time.time() - t0, 1)]

    modos.sort()
    _, aff1, cnn1, cnnaff1 = modos[0]
    mejor_aff = min(x[1] for x in modos)
    mejor_cnnaff = max(x[3] for x in modos)
    return [n, "ok", "%.2f" % aff1, "%.4f" % cnn1, "%.3f" % cnnaff1,
            "%.2f" % mejor_aff, "%.3f" % mejor_cnnaff, round(time.time() - t0, 1)]


fh = open(CSVP, "a", newline="")
w = csv.writer(fh)
if not done:
    w.writerow(CABECERA)
    fh.flush()

t0 = time.time()
n = 0
with ProcessPoolExecutor(max_workers=NJ) as ex:
    futs = [ex.submit(run, x) for x in ligs]
    for fut in as_completed(futs):
        w.writerow(fut.result())
        n += 1
        if n % 10 == 0 or n == len(ligs):
            fh.flush()
            el = time.time() - t0
            print("%d/%d  %.1f s/lig  faltan ~%.0f min"
                  % (n, len(ligs), el / n, (el / n) * (len(ligs) - n) / 60), flush=True)
fh.flush()
fh.close()
print("FIN %d moleculas en %.1f min -> %s" % (n, (time.time() - t0) / 60, CSVP), flush=True)
