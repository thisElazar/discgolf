"""Collect the sweep's converged coefficients, derive the validation metrics,
and compare against the Potts (Kamaruddin/Potts & Crowther 2018) Aviar row."""
import glob, json, os, re
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
SWEEP = os.path.join(BASE, "sweep")

# coefficient.dat columns (ESI forceCoeffs): Time Cd Cd(f) Cd(r) Cl Cl(f) Cl(r) CmPitch ...
COL = {"Cd": 1, "Cl": 4, "CmPitch": 7}

POTTS_AVIAR = {  # putter row, Re=3.78e5, non-spinning
    "ClCd_a0": 1.8, "ClCd_peak": 2.8, "dCm_dalpha_per_deg": 0.002,
    "trim_deg": 7.5, "xac_over_d": 0.05,
}


def read_case(d, navg=300):
    """Time-average the last `navg` iterations (steady RANS on a bluff wake
    oscillates; the mean is the physical answer). Returns means + scatter."""
    m = re.search(r"a_(-?\d+)$", d)
    if not m:
        return None
    a = float(m.group(1))
    dat = os.path.join(d, "postProcessing/forceCoeffs1/0/coefficient.dat")
    if not os.path.exists(dat):
        return None
    rows = [ln.split() for ln in open(dat) if ln.strip() and not ln.startswith("#")]
    if not rows:
        return None
    tail = rows[-navg:]
    cl = np.array([float(r[COL["Cl"]]) for r in tail])
    cd = np.array([float(r[COL["Cd"]]) for r in tail])
    cm = np.array([float(r[COL["CmPitch"]]) for r in tail])
    return dict(alpha=a, Cd=float(cd.mean()), Cl=float(cl.mean()), Cm=float(cm.mean()),
                Cl_sd=float(cl.std()), n_iter=len(rows), n_avg=len(tail))


def main():
    rows = sorted(filter(None, (read_case(d) for d in glob.glob(os.path.join(SWEEP, "a_*")))),
                  key=lambda r: r["alpha"])
    if not rows:
        print("no results yet"); return
    a = np.array([r["alpha"] for r in rows])
    Cl = np.array([r["Cl"] for r in rows])
    Cd = np.array([r["Cd"] for r in rows])
    Cm = np.array([r["Cm"] for r in rows])
    ClCd = Cl / Cd

    # linear region -5..+5 for gradients / aerodynamic centre
    lin = (a >= -5) & (a <= 5)
    dCl = np.polyfit(a[lin], Cl[lin], 1)[0]          # per degree
    dCm = np.polyfit(a[lin], Cm[lin], 1)[0]          # per degree
    xac_d = dCm / dCl if dCl else float("nan")
    # trim point: actual Cm=0 crossing (CM is non-linear across the full range,
    # so interpolate the real sign change rather than extrapolate a global fit)
    trim = float("nan")
    for i in range(len(a) - 1):
        if Cm[i] == 0:
            trim = a[i]; break
        if Cm[i] * Cm[i + 1] < 0:                    # sign change between i, i+1
            trim = a[i] - Cm[i] * (a[i + 1] - a[i]) / (Cm[i + 1] - Cm[i])
            break

    print(f"{'alpha':>6} {'Cl':>9} {'Cd':>9} {'Cl/Cd':>8} {'Cm':>10}")
    for r, lc in zip(rows, ClCd):
        print(f"{r['alpha']:6.0f} {r['Cl']:9.4f} {r['Cd']:9.4f} {lc:8.3f} {r['Cm']:10.4f}")

    i0 = int(np.argmin(np.abs(a)))
    metrics = {
        "ClCd_a0": round(float(ClCd[i0]), 2),
        "ClCd_peak": round(float(np.max(ClCd)), 2),
        "ClCd_peak_alpha": float(a[int(np.argmax(ClCd))]),
        "dCm_dalpha_per_deg": round(float(dCm), 4),
        "dCl_dalpha_per_deg": round(float(dCl), 4),
        "trim_deg": round(float(trim), 1),
        "xac_over_d": round(float(xac_d), 3),
        "n_cases": len(rows),
    }
    print("\n=== CFD metrics ===")
    print(json.dumps(metrics, indent=2))
    print("\n=== Potts Aviar target ===")
    print(json.dumps(POTTS_AVIAR, indent=2))
    print("\n=== comparison (CFD vs Potts) ===")
    for k in ("ClCd_a0", "ClCd_peak", "dCm_dalpha_per_deg", "trim_deg", "xac_over_d"):
        print(f"  {k:22} {metrics[k]:>8}  vs  {POTTS_AVIAR[k]:>6}")

    with open(os.path.join(BASE, "sweep_results.json"), "w") as f:
        json.dump({"rows": rows, "metrics": metrics, "potts": POTTS_AVIAR}, f, indent=2)


if __name__ == "__main__":
    main()
