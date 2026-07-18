"""Collect the big sweep into a surrogate-training table.

For each solved (geometry, alpha) case, time-average the last N iterations of
the force coefficients and join with that geometry's shape parameters. Writes
dataset.csv (one row per geometry x alpha) and prints coverage.
"""
import csv, glob, json, os, re
import numpy as np

BS = os.path.dirname(os.path.abspath(__file__))
COL = {"Cd": 1, "Cl": 4, "CmPitch": 7}
NAVG = 300


def avg_case(path):
    dat = os.path.join(path, "postProcessing/forceCoeffs1/0/coefficient.dat")
    if not os.path.exists(dat):
        return None
    rows = [ln.split() for ln in open(dat) if ln.strip() and not ln.startswith("#")]
    n = len(rows)
    if n < 150:                       # too short to have cleared the startup transient
        return None
    # average a post-transient window: always skip >= the first half, and take
    # at most the last NAVG iters. (Averaging the raw "last 300" of a case that
    # converged in ~220 iters swept in the wild iter-1 transient -> bogus means.)
    tail = rows[n - min(NAVG, n // 2):]
    return {k: float(np.mean([float(r[c]) for r in tail])) for k, c in COL.items()} | \
           {"n_iter": n}


def main():
    man = json.load(open(os.path.join(BS, "geometries.json")))["geometries"]
    gp = {g["id"]: g for g in man}
    pcols = ["rim_width_mm", "rim_shoulder_mm", "parting_line_mm",
             "dome_mm", "plate_thick_mm", "nose_radius_mm"]

    rows = []
    for d in sorted(glob.glob(os.path.join(BS, "solve", "geom_*_a*"))):
        m = re.search(r"(geom_\d+)_a(-?\d+)$", d)
        if not m:
            continue
        gid, alpha = m.group(1), float(m.group(2))
        if gid not in gp:
            continue
        c = avg_case(d)
        if c is None:
            continue
        rows.append({"id": gid, **{k: gp[gid][k] for k in pcols}, "alpha_deg": alpha,
                     "Cl": round(c["Cl"], 5), "Cd": round(c["Cd"], 5),
                     "Cm": round(c["CmPitch"], 5), "n_iter": c["n_iter"]})

    out = os.path.join(BS, "dataset.csv")
    if rows:
        with open(out, "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=list(rows[0]))
            wr.writeheader(); wr.writerows(rows)

    geoms_done = len({r["id"] for r in rows})
    print(f"rows: {len(rows)}   geometries with >=1 solve: {geoms_done}/{len(man)}")
    print(f"wrote {out}")
    # per-geometry angle coverage (a full row set = 7 angles)
    from collections import Counter
    cov = Counter(r["id"] for r in rows)
    incomplete = {g: n for g, n in cov.items() if n < 7}
    if incomplete:
        print(f"incomplete geometries (<7 angles): {len(incomplete)}")


if __name__ == "__main__":
    main()
