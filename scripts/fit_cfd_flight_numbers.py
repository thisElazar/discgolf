#!/usr/bin/env python3
"""
Flight numbers from CFD-surrogate coefficients.
===============================================

Replaces the Tier-2 geometry->rating ridge fit as the app's flight-number
source: features here are the CFD-surrogate aerodynamic coefficients (the
same cfdAero() polynomial the designer evaluates), so the four numbers are
driven by computed aerodynamics rather than raw geometry.

Same honest methodology as scripts/fit_flight_numbers.py:
- tiny ridge models (LOO-validated) against the manufacturer ratings of the
  pro discs in data/pro_discs.json;
- manufacturer ratings are a marketing scale — LOO-RMSE is the number to
  trust;
- the surrogate itself was trained at fixed 211 mm diameter; diameter
  effects on the ratings are absorbed by the fit, not modelled.

Feature glossary (alpha in radians inside the surrogate):
  CD0m  = minimum drag (drag-bucket floor)          -> how fast it cuts air
  CL0, CDa = lift offset, drag-bucket curvature
  CM0   = pitching moment at alpha=0 (high-speed attitude) -> turn
  CM5, CD5 = moment and drag at alpha=5deg (late-flight regime) -> fade
  trim  = -CM0/CMa in degrees (CM=0 crossing)        -> where fade takes over

Feature sets were chosen by LOO comparison among small physically-signed
candidates (same mild selection bias as the Tier-2 fit, documented there).

Outputs: data/cfd_flight_fit.json + cfd_flight_fit.js (CFD_FLIGHT_FIT).

Run:  python3 scripts/fit_cfd_flight_numbers.py
"""
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cfd"))
from geometry import params_from_certified  # noqa: E402  (app-identical mapping)

RIDGE_LAMBDA = 1.0
DEG = math.pi / 180.0

FEATURES = {
    "speed": ["CD0m", "CDa", "CL0"],
    "glide": ["CL0", "CD0m", "trim"],
    "turn":  ["CM0", "CL0", "CD0m"],
    "fade":  ["CM5", "CD5", "CL0"],
}
CLAMPS = {"speed": (1, 14.5), "glide": (1, 7), "turn": (-5, 1), "fade": (0, 5)}


def load_surrogate():
    txt = (ROOT / "cfd_surrogate.js").read_text()
    return json.loads(txt[txt.index("{"): txt.rindex("}") + 1])


def surrogate_coeffs(S, p):
    x = np.array([getattr(p, f) for f in S["features"]], float)
    xs = (x - np.array(S["mean"])) / np.array(S["std"])
    t = [1.0] + list(xs)
    for i in range(6):
        for j in range(i, 6):
            t.append(xs[i] * xs[j])
    t = np.array(t)
    return {k: float(np.dot(S["targets"][k], t)) for k in S["targets"]}


def aero_features(c):
    a5 = 5.0 * DEG
    trim = (-c["CM0"] / c["CMa"]) / DEG if c["CMa"] else 0.0
    return {"CD0m": c["CD0"], "CDa": c["CDa"], "CL0": c["CL0"],
            "CM0": c["CM0"], "CM5": c["CM0"] + c["CMa"] * a5,
            "CD5": c["CD0"] + c["CDa"] * (a5 - c["alpha0"]) ** 2,
            "trim": float(np.clip(trim, -20, 30))}


def ridge_loo(X, y, lam=RIDGE_LAMBDA):
    """Standardized ridge with leave-one-out CV. Returns weights on RAW
    features (w, b) plus train R2 and LOO RMSE."""
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd
    n, d = Xs.shape

    def fit(idx):
        A = Xs[idx].T @ Xs[idx] + lam * np.eye(d)
        w = np.linalg.solve(A, Xs[idx].T @ (y[idx] - y[idx].mean()))
        b = y[idx].mean()
        return w, b

    preds = np.empty(n)
    for i in range(n):
        idx = np.array([j for j in range(n) if j != i])
        w, b = fit(idx)
        preds[i] = Xs[i] @ w + b
    w, b = fit(np.arange(n))
    yhat = Xs @ w + b
    ss_res = np.sum((y - yhat) ** 2); ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot else float("nan")
    loo = float(np.sqrt(np.mean((y - preds) ** 2)))
    w_raw = w / sd
    b_raw = float(b - np.sum(w_raw * mu))
    return w_raw, b_raw, float(r2), loo


def main():
    S = load_surrogate()
    discs = json.load(open(ROOT / "data" / "pro_discs.json"))["discs"]
    rows = []
    for d in discs:
        if not d.get("flight"):
            continue
        p = params_from_certified(d)
        feats = aero_features(surrogate_coeffs(S, p))
        rows.append({"name": d["name"], **feats,
                     "speed": d["flight"][0], "glide": d["flight"][1],
                     "turn": d["flight"][2], "fade": d["flight"][3]})
    print(f"{len(rows)} pro discs with ratings")

    out = {"method": f"ridge (lambda={RIDGE_LAMBDA} on standardized features) on "
                     "CFD-surrogate aero coefficients vs manufacturer ratings",
           "caveats": [
               "targets are manufacturer marketing ratings, not measurements",
               "surrogate trained at fixed 211 mm diameter",
               "trust loo_rmse, not train_r2",
           ],
           "n_discs": len(rows), "fit": {}}
    print(f"{'target':6} {'train_r2':>8} {'loo_rmse':>8}   features")
    for tgt, fs in FEATURES.items():
        X = np.array([[r[f] for f in fs] for r in rows], float)
        y = np.array([r[tgt] for r in rows], float)
        w, b, r2, loo = ridge_loo(X, y)
        out["fit"][tgt] = {"features": {f: round(float(v), 6) for f, v in zip(fs, w)},
                           "intercept": round(b, 6), "clamp": list(CLAMPS[tgt]),
                           "train_r2": round(r2, 3), "loo_rmse": round(loo, 3)}
        print(f"{tgt:6} {r2:8.3f} {loo:8.3f}   {fs}")

    json.dump(out, open(ROOT / "data" / "cfd_flight_fit.json", "w"), indent=1)
    js = ("// Generated by scripts/fit_cfd_flight_numbers.py — do not edit by hand.\n"
          "// Flight numbers from CFD-surrogate aero coefficients, calibrated to the\n"
          f"// manufacturer ratings of {len(rows)} pro discs; see data/cfd_flight_fit.json.\n"
          "const CFD_FLIGHT_FIT = " + json.dumps(out["fit"]) + ";\n")
    (ROOT / "cfd_flight_fit.js").write_text(js)
    print(f"wrote cfd_flight_fit.js + data/cfd_flight_fit.json")


if __name__ == "__main__":
    main()
