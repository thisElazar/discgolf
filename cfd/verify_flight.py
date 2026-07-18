"""Verify the JS flight port against flight_sim.py on identical coefficients.

Computes the CFD-surrogate coefficients for a test disc (same degree-2 poly
the JS uses), feeds them into flight_sim.py's 6-DOF integrator, and writes the
test case for the Node side to run through cfd_flight.js. Endpoints must match.
"""
import json, os, re, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import disc_model as dm
import flight_sim as fs

# load the exported surrogate (strip the JS wrapper -> JSON)
txt = open(os.path.join(ROOT, "cfd_surrogate.js")).read()
S = json.loads(txt[txt.index("{"):txt.rindex("}") + 1])
FEATS = S["features"]; mean = np.array(S["mean"]); std = np.array(S["std"])


def surrogate_coeffs(p):
    x = np.array([getattr(p, f) for f in FEATS])
    xs = (x - mean) / std
    t = [1.0] + list(xs)
    for i in range(6):
        for j in range(i, 6):
            t.append(xs[i] * xs[j])
    t = np.array(t)
    c = {k: float(np.dot(S["targets"][k], t)) for k in S["targets"]}
    return c


def main():
    p = dm.DiscParams(D_out_mm=211, rim_width_mm=18.0, rim_shoulder_mm=15.0,
                      parting_line_mm=9.0, dome_mm=2.0, plate_thick_mm=2.0,
                      nose_radius_mm=2.0, density_gcc=0.95)
    mp = dm.mass_properties(p)
    c = surrogate_coeffs(p)

    # inject surrogate coeffs into the sim in place of aero_from_geometry
    ac = fs.AeroCoeffs(CL0=c["CL0"], CLa=c["CLa"], CD0=c["CD0"], CDa=c["CDa"],
                       alpha0=c["alpha0"], CM0=c["CM0"], CMa=c["CMa"],
                       CMq=-5.0e-3, CRr=1.4e-2, CRp=-5.5e-3, CNr=-3.4e-5)
    fs.aero_from_geometry = lambda p, mp=None: ac

    opts = dict(v0=24.0, launch_deg=8.0, hyzer_deg=0.0, nose_deg=0.0,
                h0=1.4, dt=1e-3, handed="RHBH")
    res = fs.simulate(p, **opts)

    case = {"params": {k: getattr(p, k) for k in
                       ("D_out_mm", "rim_width_mm", "rim_shoulder_mm", "parting_line_mm",
                        "dome_mm", "plate_thick_mm", "nose_radius_mm")},
            "mp": {"Iz_gcm2": mp["Iz_gcm2"], "mass_g": mp["mass_g"]},
            "opts": opts, "coeffs": c}
    json.dump(case, open(os.path.join(ROOT, "cfd", "test_case.json"), "w"), indent=2)

    print("PYTHON (flight_sim.py):")
    print("  coeffs:", {k: round(v, 5) for k, v in c.items()})
    print(f"  distance_ft = {res['distance_ft']:.2f}")
    print(f"  lateral_m   = {res['lateral_max_m']:+.3f}")
    print(f"  apex_m      = {res['apex_m']:.3f}")
    print(f"  time_s      = {res['flight_time_s']:.3f}")


if __name__ == "__main__":
    main()
