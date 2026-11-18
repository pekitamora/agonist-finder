import csv, os, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed

EXE, REC, LIGD, POSED, CSVP = sys.argv[1:6]
NJ = int(sys.argv[6]) if len(sys.argv) > 6 else 2
EXH = sys.argv[7] if len(sys.argv) > 7 else "8"
EXTRA = sys.argv[8].split() if len(sys.argv) > 8 else []
CEN = ("2.321", "6.301", "32.385")
SZ = ("25.7", "26.3", "29.0")

os.makedirs(POSED, exist_ok=True)
os.makedirs(os.path.dirname(CSVP) or ".", exist_ok=True)

done = set()
if os.path.exists(CSVP):
    done = {r[0] for r in csv.reader(open(CSVP)) if r and r[0] != "ligando"}
todos = sorted(f for f in os.listdir(LIGD) if f.endswith(".pdbqt"))
ligs = [f for f in todos if f not in done]
print("ligandos %d | ya hechos %d | pendientes %d" % (len(todos), len(done), len(ligs)), flush=True)


def run(f):
    out = os.path.join(POSED, f.replace(".pdbqt", "_out.pdbqt"))
    cmd = [EXE, "--receptor", REC, "--ligand", os.path.join(LIGD, f),
           "--center_x", CEN[0], "--center_y", CEN[1], "--center_z", CEN[2],
           "--size_x", SZ[0], "--size_y", SZ[1], "--size_z", SZ[2],
           "--cpu", "1", "--exhaustiveness", EXH, "--num_modes", "9",
           "--seed", "42", "--out", out] + EXTRA
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    except Exception:
        return [f, "timeout", "", round(time.time() - t0, 1)]
    if p.returncode != 0:
        return [f, "error", "", round(time.time() - t0, 1)]
    s = ""
    for l in p.stdout.splitlines():
        q = l.split()
        if len(q) >= 2 and q[0] == "1":
            try:
                s = "%.3f" % float(q[1])
                break
            except ValueError:
                pass
    return [f, "ok" if s else "sin_score", s, round(time.time() - t0, 1)]


fh = open(CSVP, "a", newline="")
w = csv.writer(fh)
if not done:
    w.writerow(["ligando", "estado", "score", "segundos"])
    fh.flush()

t0 = time.time()
n = 0
with ProcessPoolExecutor(max_workers=NJ) as ex:
    futs = [ex.submit(run, f) for f in ligs]
    for fut in as_completed(futs):
        w.writerow(fut.result())
        n += 1
        if n % 25 == 0 or n == len(ligs):
            fh.flush()
            el = time.time() - t0
            print("%d/%d  %.1f s/lig  faltan ~%.0f min"
                  % (n, len(ligs), el / n, (el / n) * (len(ligs) - n) / 60), flush=True)
fh.flush()
fh.close()
print("FIN %d ligandos en %.1f min -> %s" % (n, (time.time() - t0) / 60, CSVP), flush=True)
