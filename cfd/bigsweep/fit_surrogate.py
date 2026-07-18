"""Fit the geometry -> aero-coefficient surrogate from the big-sweep dataset.

Two stages:
  1. Per geometry, fit the CFD alpha-sweep to the Hummel coefficient FORM the
     6-DOF sim (flight_sim.py AeroCoeffs) expects (alpha in RADIANS):
         CL = CL0 + CLa*alpha
         CD = CD0 + CDa*(alpha - alpha0)^2
         CM = CM0 + CMa*alpha
     -> 7 static coefficients per geometry: CL0,CLa,CD0,CDa,alpha0,CM0,CMa.
  2. Train a regression-tree ensemble mapping the 6 shape parameters ->
     each coefficient (RTE beat polynomials in the 2021 Springer paper).
     Cross-validated (K-fold). Saves models + a JSON export for later JS port.

Note on absolute level: the validated coarse mesh carries a systematic
~15-20% CL/CD bias (wall-function regime, small cluster). The surrogate
learns geometry TRENDS well; anchor the absolute level against the discs
Potts measured (Aviar/Roc/Wraith) before trusting magnitudes downstream.
"""
import json, os
import numpy as np
import pandas as pd

BS = os.path.dirname(os.path.abspath(__file__))
PCOLS = ["rim_width_mm", "rim_shoulder_mm", "parting_line_mm",
         "dome_mm", "plate_thick_mm", "nose_radius_mm"]
TARGETS = ["CL0", "CLa", "CD0", "CDa", "alpha0", "CM0", "CMa"]


def fit_coeffs(g):
    """g: rows for one geometry (alpha_deg, Cl, Cd, Cm). Returns coefficient dict."""
    a = np.radians(g["alpha_deg"].to_numpy())
    Cl, Cd, Cm = g["Cl"].to_numpy(), g["Cd"].to_numpy(), g["Cm"].to_numpy()
    # linear CL, CM
    CLa, CL0 = np.polyfit(a, Cl, 1)
    CMa, CM0 = np.polyfit(a, Cm, 1)
    # quadratic CD -> CD0 + CDa*(a-a0)^2
    qa, qb, qc = np.polyfit(a, Cd, 2)
    CDa = qa
    alpha0 = -qb / (2 * qa) if qa else 0.0
    CD0 = qc - qb * qb / (4 * qa) if qa else float(np.min(Cd))
    return dict(CL0=CL0, CLa=CLa, CD0=CD0, CDa=CDa, alpha0=alpha0, CM0=CM0, CMa=CMa)


def main():
    df = pd.read_csv(os.path.join(BS, "dataset.csv"))
    # keep geometries with a full-enough alpha set to fit the curves
    rows = []
    for gid, g in df.groupby("id"):
        if len(g) < 5:
            continue
        c = fit_coeffs(g.sort_values("alpha_deg"))
        rec = {"id": gid, **{k: g[k].iloc[0] for k in PCOLS}, **c}
        rows.append(rec)
    fit = pd.DataFrame(rows)
    fit.to_csv(os.path.join(BS, "coeff_table.csv"), index=False)
    print(f"fitted coefficients for {len(fit)} geometries -> coeff_table.csv")

    try:
        from sklearn.ensemble import RandomForestRegressor
        from sklearn.model_selection import cross_val_predict, KFold
    except ImportError:
        print("scikit-learn not installed; run: micromamba install -p ~/cfd/env scikit-learn")
        return

    X = fit[PCOLS].to_numpy()
    kf = KFold(n_splits=8, shuffle=True, random_state=0)
    models, report = {}, {}
    for t in TARGETS:
        y = fit[t].to_numpy()
        rf = RandomForestRegressor(n_estimators=400, min_samples_leaf=2, random_state=0)
        yp = cross_val_predict(rf, X, y, cv=kf)
        ss_res = np.sum((y - yp) ** 2); ss_tot = np.sum((y - y.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot else float("nan")
        rmse = float(np.sqrt(np.mean((y - yp) ** 2)))
        rf.fit(X, y)
        models[t] = rf
        report[t] = {"cv_r2": round(float(r2), 3), "cv_rmse": round(rmse, 5)}

    print("\n=== surrogate cross-validation (8-fold) ===")
    for t in TARGETS:
        print(f"  {t:7} R2={report[t]['cv_r2']:>6}  RMSE={report[t]['cv_rmse']}")

    import joblib
    joblib.dump({"models": models, "features": PCOLS, "targets": TARGETS},
                os.path.join(BS, "surrogate.joblib"))
    json.dump({"features": PCOLS, "targets": TARGETS, "cv": report,
               "n_geometries": len(fit)},
              open(os.path.join(BS, "surrogate_report.json"), "w"), indent=2)
    print("\nsaved surrogate.joblib + surrogate_report.json")


if __name__ == "__main__":
    main()
