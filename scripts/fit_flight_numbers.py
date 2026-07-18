#!/usr/bin/env python3
"""
Calibrate the app's live flight-number estimate against real discs.
================================================================

Fits one small ridge-regression model per flight number (speed, glide,
turn, fade) mapping geometry features -> the manufacturer's printed rating,
using the pro discs in data/pro_discs.json (PDGA-certified measurements
+ manufacturer flight numbers).

Methodology, honestly stated
----------------------------
- n is small (see data/pro_discs.json). Models are deliberately tiny (2-4 features + intercept)
  and ridge-regularized; anything bigger overfits.
- Validation is leave-one-out CV (LOO): each disc is predicted by a model
  trained on the other 15. LOO-RMSE is the number to trust, not train R².
- Targets are MANUFACTURER RATINGS — a marketing scale, not a measurement.
  The fit inherits that scale's noise and brand inconsistencies.
- The `pr` feature (parting-line height / shoulder height) is CIRCULAR for
  turn/fade on this dataset: preset parting lines were themselves derived
  from the printed stability. Including it encodes our stability assumption
  (documented in METHODS.md); it is NOT independent evidence. The measured
  features (rim width, height, cavity depth, diameter) are grounded.

Outputs
-------
- data/flight_fit.json  — coefficients + fit statistics (provenance record)
- flight_fit.js         — `const FLIGHT_FIT = {...}` consumed by the app

Run:  python3 scripts/fit_flight_numbers.py
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from disc_model import DiscParams, mass_properties  # noqa: E402

RIDGE_LAMBDA = 1.0  # on standardized features

# Per-class estimates for dimensions the PDGA doesn't publish.
# MUST match CLASS_EST / paramsFromCertified in disc_designer.html.
CLASS_EST = {
    "Putter":   dict(plate=3.2, nose=3.0),
    "Midrange": dict(plate=2.6, nose=2.5),
    "Fairway":  dict(plate=2.0, nose=2.0),
    "Distance": dict(plate=1.7, nose=1.7),
}

# Feature sets per target: small and physically signed.
# rw = rim width, dome = dome height, h = total height, cav = cavity depth,
# pr = parting-line/shoulder ratio (assumption-encoded — see docstring).
FEATURES = {
    "speed": ["rw", "dia"],
    "glide": ["rw", "dia", "pr", "dome"],
    "turn":  ["pr", "h", "rw"],
    "fade":  ["rw", "pr", "dome"],
}
CLAMPS = {"speed": (1, 14.5), "glide": (1, 7), "turn": (-5, 1), "fade": (0, 5)}


def params_from_certified(d: dict) -> DiscParams:
    """Mirror of paramsFromCertified() in disc_designer.html."""
    est = CLASS_EST[d["cls"]]
    shoulder = d["rd"] + est["plate"]
    dome = float(np.clip(d["h"] - shoulder, 0.2, 8))
    stab = d["flight"][3] + d["flight"][2]
    parting = float(np.clip(shoulder * (0.60 - 0.04 * stab), 5, shoulder - 2))
    p = DiscParams(D_out_mm=d["dia"], rim_width_mm=d["rw"],
                   rim_shoulder_mm=round(shoulder, 1),
                   parting_line_mm=round(parting, 1), dome_mm=round(dome, 1),
                   plate_thick_mm=est["plate"], nose_radius_mm=est["nose"],
                   density_gcc=1.0)
    target = min(d["wt"], 8.3 * d["dia"] / 10) - 0.5
    vol = mass_properties(p)["volume_cm3"]
    p.density_gcc = float(np.clip(target / max(vol, 1e-9), 0.85, 1.00))
    return p


def features_of(p: DiscParams) -> dict:
    mp = mass_properties(p)
    return dict(rw=p.rim_width_mm, dia=p.D_out_mm, dome=p.dome_mm,
                h=mp["total_height_mm"], cav=mp["cavity_depth_mm"],
                pr=p.parting_line_mm / p.rim_shoulder_mm)


def ridge_fit(X: np.ndarray, y: np.ndarray, lam: float):
    """Ridge on standardized features; returns raw-space (coefs, intercept)."""
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd
    A = Xs.T @ Xs + lam * np.eye(X.shape[1])
    b = Xs.T @ (y - y.mean())
    w_s = np.linalg.solve(A, b)
    w = w_s / sd
    return w, float(y.mean() - w @ mu)


def main() -> None:
    data = json.loads((ROOT / "data" / "pro_discs.json").read_text())["discs"]
    feats = [features_of(params_from_certified(d)) for d in data]
    names = [f'{d["mfr"]} {d["name"]}' for d in data]

    fit, report = {}, []
    for ti, target in enumerate(["speed", "glide", "turn", "fade"]):
        cols = FEATURES[target]
        X = np.array([[f[c] for c in cols] for f in feats])
        y = np.array([d["flight"][ti] for d in data], dtype=float)

        w, b = ridge_fit(X, y, RIDGE_LAMBDA)
        pred = X @ w + b
        ss_res = ((y - pred) ** 2).sum()
        ss_tot = ((y - y.mean()) ** 2).sum()
        r2 = 1 - ss_res / ss_tot

        # leave-one-out CV
        loo = np.empty(len(y))
        for i in range(len(y)):
            m = np.ones(len(y), bool); m[i] = False
            wi, bi = ridge_fit(X[m], y[m], RIDGE_LAMBDA)
            loo[i] = X[i] @ wi + bi
        loo_rmse = float(np.sqrt(((y - loo) ** 2).mean()))

        lo, hi = CLAMPS[target]
        fit[target] = dict(features={c: round(float(wc), 5) for c, wc in zip(cols, w)},
                           intercept=round(b, 5), clamp=[lo, hi],
                           train_r2=round(float(r2), 3), loo_rmse=round(loo_rmse, 3))
        report.append((target, r2, loo_rmse, y, pred, loo))

    # ---- provenance JSON ----
    out = dict(
        method=f"ridge regression (lambda={RIDGE_LAMBDA} on standardized features), "
               f"fit to manufacturer flight numbers of the {len(data)} discs in pro_discs.json",
        caveats=[
            f"n={len(data)}; trust loo_rmse, not train_r2",
            "targets are manufacturer marketing ratings, not measurements",
            "feature 'pr' is derived from printed stability for these presets (circular) — "
            "it encodes our parting-line assumption, not independent evidence",
            "feature sets were chosen partly by LOO comparison on the same 16 discs "
            "(mild selection bias); glide keeps a small dome term at slight LOO cost "
            "because the physical effect (dome -> glide) is expected and it keeps the "
            "designer's dome slider live",
        ],
        features_glossary=dict(
            rw="rim width (mm)", dia="outside diameter (mm)", dome="dome height (mm)",
            h="total height (mm)", cav="cavity depth (mm)",
            pr="parting-line height / shoulder height"),
        fit=fit,
    )
    (ROOT / "data" / "flight_fit.json").write_text(json.dumps(out, indent=2) + "\n")

    (ROOT / "flight_fit.js").write_text(
        "// Generated by scripts/fit_flight_numbers.py — do not edit by hand.\n"
        f"// Calibrated to manufacturer flight numbers of {len(data)} pro discs; "
        "see data/flight_fit.json for stats & caveats.\n"
        f"const FLIGHT_FIT = {json.dumps(fit, indent=1)};\n")

    # ---- console report ----
    print(f"{'':24s}  rated        fit          LOO")
    for i, n in enumerate(names):
        rated = "/".join(f"{report[t][3][i]:g}" for t in range(4))
        fitv = "/".join(f"{report[t][4][i]:.1f}" for t in range(4))
        loov = "/".join(f"{report[t][5][i]:.1f}" for t in range(4))
        print(f"{n:24s}  {rated:12s} {fitv:12s} {loov}")
    print()
    for target, r2, loo_rmse, *_ in report:
        print(f"{target:6s} train R2 = {r2:5.3f}   LOO-RMSE = {loo_rmse:.2f}  "
              f"(features: {', '.join(FEATURES[target])})")
    print("\nwrote data/flight_fit.json and flight_fit.js")


if __name__ == "__main__":
    main()
