"""Fit the anchor-disc CFD sweeps and derive the CL/CD bias correction.

Compares our coarse-mesh CFD (Aviar/Roc/Wraith, same pipeline as the big
sweep) against Kamaruddin's wind-tunnel measurements of the same discs
(data/kamaruddin_wind_tunnel.json) and derives multiplicative correction
factors kCL, kCD that map the CFD level onto the measured level:

    kCL = least-squares scale between CFD CL(alpha) and measured CL(alpha)
    kCD = same for CD(alpha), on the sweep's alpha grid, pooled over discs

CM is compared but not corrected unless it disagrees (the earlier validation
showed trim/gradient match). Writes anchor_bias.json.

Run locally after syncing anchors/solve back:  python3 cfd/anchors/fit_anchors.py
"""
import glob, json, math, os, re, sys

import numpy as np

AN = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(AN))
DEG = math.pi / 180.0
COL = {"Cd": 1, "Cl": 4, "CmPitch": 7}
NAVG = 300
ALPHAS = [-4.0, -1.0, 2.0, 5.0, 8.0, 11.0, 14.0]
NAME = {"anchor_aviar": "Aviar", "anchor_roc": "Roc", "anchor_wraith": "Wraith"}


def avg_case(path):
    """Same post-transient averaging as bigsweep/collect.py."""
    dat = os.path.join(path, "postProcessing/forceCoeffs1/0/coefficient.dat")
    if not os.path.exists(dat):
        return None
    rows = [ln.split() for ln in open(dat) if ln.strip() and not ln.startswith("#")]
    n = len(rows)
    if n < 100:        # residual-control runs can converge well before 1000 iters
        return None
    tail = rows[n - min(NAVG, n // 2):]
    return {k: float(np.mean([float(r[c]) for r in tail])) for k, c in COL.items()}


def fit_hummel(a_deg, Cl, Cd, Cm):
    """Same coefficient form as bigsweep/fit_surrogate.py (alpha in radians)."""
    a = np.radians(a_deg)
    CLa, CL0 = np.polyfit(a, Cl, 1)
    CMa, CM0 = np.polyfit(a, Cm, 1)
    qa, qb, qc = np.polyfit(a, Cd, 2)
    alpha0 = -qb / (2 * qa) if qa else 0.0
    CD0 = qc - qb * qb / (4 * qa) if qa else float(np.min(Cd))
    return dict(CL0=CL0, CLa=CLa, CD0=CD0, CDa=qa, alpha0=alpha0, CM0=CM0, CMa=CMa)


def main():
    meas = json.load(open(os.path.join(ROOT, "data", "kamaruddin_wind_tunnel.json")))["discs"]

    cfd = {}
    for d in sorted(glob.glob(os.path.join(AN, "solve", "anchor_*_a*"))):
        m = re.search(r"(anchor_\w+)_a(-?\d+)$", d)
        if not m:
            continue
        c = avg_case(d)
        if c is None:
            print(f"  (no converged data: {os.path.basename(d)})")
            continue
        cfd.setdefault(m.group(1), {})[float(m.group(2))] = c

    grid = np.array(ALPHAS)
    num_L = den_L = num_D = den_D = 0.0
    report = {}
    for gid, by_a in sorted(cfd.items()):
        name = NAME[gid]
        a = np.array(sorted(by_a))
        Cl = np.array([by_a[x]["Cl"] for x in sorted(by_a)])
        Cd = np.array([by_a[x]["Cd"] for x in sorted(by_a)])
        Cm = np.array([by_a[x]["CmPitch"] for x in sorted(by_a)])
        fit = fit_hummel(a, Cl, Cd, Cm)

        mr = meas[name]
        mCl = mr["CL0"] + mr["dCL_dalpha"] * grid
        mCd_0 = mr["CD0"]                      # measured CD at alpha=0 (table)
        cCl = np.interp(grid, a, Cl)
        cCd = np.interp(grid, a, Cd)

        num_L += float(np.dot(cCl, mCl)); den_L += float(np.dot(cCl, cCl))
        # CD anchored at the alpha=0 table value only (curvature not tabulated)
        cCd0 = float(np.interp(0.0, a, Cd))
        num_D += cCd0 * mCd_0; den_D += cCd0 * cCd0

        trim_cfd = -fit["CM0"] / fit["CMa"] / DEG if fit["CMa"] else float("nan")
        report[name] = {
            "cfd_fit": {k: round(float(v), 5) for k, v in fit.items()},
            "cfd_CL_at_0": round(float(np.interp(0.0, a, Cl)), 4),
            "meas_CL0": mr["CL0"],
            "cfd_CD_at_0": round(cCd0, 4),
            "meas_CD0": mr["CD0"],
            "cfd_CLa_per_deg": round(float(fit["CLa"]) * DEG, 4),
            "meas_CLa_per_deg": mr["dCL_dalpha"],
            "cfd_CMa_per_deg": round(float(fit["CMa"]) * DEG, 4),
            "meas_CMa_per_deg": mr["dCM_dalpha"],
            "cfd_trim_deg": round(trim_cfd, 1),
            "meas_trim_deg": mr["trim_deg"],
        }

    kCL = num_L / den_L if den_L else float("nan")
    kCD = num_D / den_D if den_D else float("nan")

    print(json.dumps(report, indent=2))
    print(f"\nbias correction (measured = k * CFD):  kCL = {kCL:.3f}   kCD = {kCD:.3f}")

    out = {"kCL": round(kCL, 4), "kCD": round(kCD, 4),
           "method": "least-squares scale CFD->measured; CL over alpha grid "
                     f"{ALPHAS}, CD at alpha=0 (Kamaruddin Table 5.6); "
                     "pooled over Aviar/Roc/Wraith",
           "per_disc": report}
    json.dump(out, open(os.path.join(AN, "anchor_bias.json"), "w"), indent=2)
    print(f"wrote {os.path.join(AN, 'anchor_bias.json')}")


if __name__ == "__main__":
    main()
