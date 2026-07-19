"""Full-loop sanity check: disc geometry -> (bias-corrected) CFD surrogate ->
6-DOF trajectory, for the three Kamaruddin anchor discs.

This closes the whole pipeline the designer app uses: the same parametric
geometry, the same exported surrogate polynomial (cfd_surrogate.js), the same
integrator. Compares the surrogate's coefficients against Kamaruddin's
measured values (Table 5.6) and the resulting ranges against her published
6-DOF benchmark (49 / 54 / 63 m; assumptions-matched config).

Run:  python3 cfd/verify_anchors_flight.py
"""
import json, math, os, sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "cfd"))
sys.path.insert(0, os.path.join(ROOT, "cfd", "anchors"))
import flight_sim as fs
import disc_model as dm
from gen_anchors import ANCHORS
from geometry import params_from_certified
from verify_trajectory import kamaruddin_phys, BENCH

DEG = math.pi / 180.0
NAME = {"anchor_aviar": "Aviar", "anchor_roc": "Roc", "anchor_wraith": "Wraith"}


def surrogate():
    txt = open(os.path.join(ROOT, "cfd_surrogate.js")).read()
    return json.loads(txt[txt.index("{"): txt.rindex("}") + 1])


def surrogate_coeffs(S, p):
    x = np.array([getattr(p, f) for f in S["features"]], float)
    xs = (x - np.array(S["mean"])) / np.array(S["std"])
    t = [1.0] + list(xs)
    for i in range(6):
        for j in range(i, 6):
            t.append(xs[i] * xs[j])
    return {k: float(np.dot(S["targets"][k], np.array(t))) for k in S["targets"]}


def main():
    S = surrogate()
    meas = json.load(open(os.path.join(ROOT, "data", "kamaruddin_wind_tunnel.json")))["discs"]
    k = S.get("bias_correction", {})
    print(f"surrogate bias correction: kCL={k.get('kCL', 1.0)}  kCD={k.get('kCD', 1.0)}\n")

    print(f"{'disc':7} {'CL@0 sur/meas':>14} {'CD@0 sur/meas':>14} "
          f"{'trim sur/meas':>14} {'range m':>8} {'bench':>6}")
    orig_aero, orig_phys = fs.aero_from_geometry, fs.physical_from_geometry
    orig_gain = fs.LOWSPEED_FADE_GAIN
    try:
        for gid, cert in ANCHORS.items():
            name = NAME[gid]
            p = params_from_certified(cert)
            c = surrogate_coeffs(S, p)
            # assumptions-matched config (same as verify_trajectory.py)
            ac = fs.AeroCoeffs(CL0=c["CL0"], CLa=c["CLa"], CD0=c["CD0"],
                               CDa=c["CDa"], alpha0=c["alpha0"], CM0=c["CM0"],
                               CMa=c["CMa"], CMq=0.0, CRr=0.0, CRp=0.0, CNr=0.0)
            fs.aero_from_geometry = lambda pp, mp=None, ac=ac: ac
            fs.physical_from_geometry = lambda pp, mp=None: kamaruddin_phys()
            fs.LOWSPEED_FADE_GAIN = 0.0
            phys = kamaruddin_phys()
            v0 = 20.0
            res = fs.simulate(dm.DiscParams(), v0=v0, spin_rev=(v0 / phys.d) / (2 * math.pi),
                              launch_deg=15.0, hyzer_deg=0.0, nose_deg=0.0,
                              h0=0.01, dt=1e-3, handed="RHBH")
            rng = float(math.hypot(res["x"][-1], res["y"][-1]))
            m = meas[name]
            trim_sur = (-c["CM0"] / c["CMa"]) / DEG if c["CMa"] else float("nan")
            print(f"{name:7} {c['CL0']:6.3f}/{m['CL0']:5.3f} "
                  f"{c['CD0'] + c['CDa'] * c['alpha0'] ** 2:7.3f}/{m['CD0']:5.3f} "
                  f"{trim_sur:6.1f}/{m['trim_deg']:4.0f} {rng:8.1f} {BENCH[name]:6.0f}")
    finally:
        fs.aero_from_geometry, fs.physical_from_geometry = orig_aero, orig_phys
        fs.LOWSPEED_FADE_GAIN = orig_gain


if __name__ == "__main__":
    main()
