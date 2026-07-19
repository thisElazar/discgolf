"""Validate the 6-DOF trajectory model against Kamaruddin's published flights.

Kamaruddin (2011) ch.6 / Kamaruddin, Potts & Crowther (2018) Fig 7 simulate
the Aviar/Roc/Wraith with their MEASURED wind-tunnel coefficients and publish
the resulting ranges (putter 49 m, mid-range 54 m, driver 63 m) for a fixed
launch: V=20 m/s, pitch 15 deg, roll 0, alpha 0, AdvR 0.5, ground-level
launch, m=0.177 kg, d=0.214 m, uniform-disc inertia. Driving OUR integrator
(flight_sim.py) with THEIR measured coefficients isolates the trajectory
model from our CFD: if ranges and flight shape match, the 6-DOF is sound.

Measured coefficient sources (data/kamaruddin_wind_tunnel.json):
  CL, CM linear fits + CD at alpha=0 from thesis Table 5.6. CD curvature is
  closed with the alpha=15 deg value read off thesis Fig 5.9(b)/5.11(b)
  (+-0.02 figure-read tolerance), vertex at the zero-lift angle.

Run:  python3 cfd/verify_trajectory.py
"""
import json, math, os, sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import flight_sim as fs
import disc_model as dm

DEG = math.pi / 180.0
# CD at alpha=15deg read from thesis Fig 5.9(b)/5.11(b)
CD15 = {"Aviar": 0.28, "Roc": 0.23, "Wraith": 0.23}
BENCH = {"Aviar": 49.0, "Roc": 54.0, "Wraith": 63.0}   # published ranges, m


def measured_aero(rec, cd15):
    """Kamaruddin Table 5.6 linear fits -> Hummel-form AeroCoeffs (radians)."""
    CL0 = rec["CL0"]
    CLa = rec["dCL_dalpha"] / DEG
    CM0 = rec["CM0"]
    CMa = rec["dCM_dalpha"] / DEG
    alpha0 = -CL0 / CLa                       # zero-lift angle = drag-bucket vertex
    # quadratic through (0, CD0_table) and (15deg, cd15) with vertex at alpha0
    a15 = 15.0 * DEG
    CDa = (cd15 - rec["CD0"]) / ((a15 - alpha0) ** 2 - alpha0 ** 2)
    CD0 = rec["CD0"] - CDa * alpha0 ** 2
    return fs.AeroCoeffs(CL0=CL0, CLa=CLa, CD0=CD0, CDa=CDa, alpha0=alpha0,
                         CM0=CM0, CMa=CMa, CMq=-5.0e-3, CRr=1.4e-2,
                         CRp=-5.5e-3, CNr=-3.4e-5)


def kamaruddin_phys():
    """Their sim: m=0.177 kg, d=0.214 m, uniform thin disc (I = m d^2 / 8)."""
    m, d = 0.177, 0.214
    Ia = m * d * d / 8.0
    return fs.DiscPhysical(m=m, d=d, A=math.pi * (d / 2) ** 2, Ia=Ia, Id=Ia / 2)


def run_disc(ac, fade_gain, matched=False):
    """matched=True strips the model terms Kamaruddin's sim does not have
    (aero damping moments, spin decay) for an assumptions-matched comparison."""
    if matched:
        ac = fs.AeroCoeffs(CL0=ac.CL0, CLa=ac.CLa, CD0=ac.CD0, CDa=ac.CDa,
                           alpha0=ac.alpha0, CM0=ac.CM0, CMa=ac.CMa,
                           CMq=0.0, CRr=0.0, CRp=0.0, CNr=0.0)
    fs.aero_from_geometry = lambda p, mp=None: ac
    fs.physical_from_geometry = lambda p, mp=None: kamaruddin_phys()
    fs.LOWSPEED_FADE_GAIN = fade_gain
    phys = kamaruddin_phys()
    v0 = 20.0
    omega = 0.5 * 2.0 * v0 / phys.d           # AdvR = omega*d/(2V) = 0.5
    p = dm.DiscParams()                       # placeholder; phys/aero are patched
    res = fs.simulate(p, v0=v0, spin_rev=omega / (2 * math.pi), launch_deg=15.0,
                      hyzer_deg=0.0, nose_deg=0.0, h0=0.01, dt=1e-3, handed="RHBH")
    x, y = res["x"], res["y"]
    rng = float(math.hypot(x[-1], y[-1]))
    # lateral reversal (S-path): does y change direction of travel mid-flight?
    dy = np.diff(y)
    reversed_ = bool(np.any(dy[:-1] * dy[1:] < 0))
    return dict(range_m=rng, apex_m=res["apex_m"], t_s=res["flight_time_s"],
                lateral_end_m=float(y[-1]), s_path=reversed_)


def main():
    data = json.load(open(os.path.join(ROOT, "data", "kamaruddin_wind_tunnel.json")))
    orig_aero, orig_phys = fs.aero_from_geometry, fs.physical_from_geometry
    orig_gain = fs.LOWSPEED_FADE_GAIN
    out = {}
    try:
        for fade_gain, matched, label in (
                (0.0, True,  "assumptions matched to Kamaruddin (no damping, constant spin)"),
                (0.0, False, "full model, no low-speed fade"),
                (0.060, False, "full model with app low-speed fade term")):
            print(f"\n=== {label} ===")
            print(f"{'disc':8} {'range m':>8} {'bench':>6} {'err%':>6} "
                  f"{'apex m':>7} {'t s':>5} {'lat m':>6}  S-path")
            out[label] = {}
            for name in ("Aviar", "Roc", "Wraith"):
                ac = measured_aero(data["discs"][name], CD15[name])
                r = run_disc(ac, fade_gain, matched)
                err = 100.0 * (r["range_m"] - BENCH[name]) / BENCH[name]
                out[label][name] = {**r, "bench_m": BENCH[name], "err_pct": round(err, 1)}
                print(f"{name:8} {r['range_m']:8.1f} {BENCH[name]:6.0f} {err:+6.1f} "
                      f"{r['apex_m']:7.2f} {r['t_s']:5.2f} {r['lateral_end_m']:+6.1f}  {r['s_path']}")
    finally:
        fs.aero_from_geometry, fs.physical_from_geometry = orig_aero, orig_phys
        fs.LOWSPEED_FADE_GAIN = orig_gain

    json.dump(out, open(os.path.join(ROOT, "cfd", "trajectory_validation.json"), "w"),
              indent=2)
    print("\nwrote cfd/trajectory_validation.json")


if __name__ == "__main__":
    main()
