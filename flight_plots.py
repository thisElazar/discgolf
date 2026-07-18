"""Simulate the four disc presets and plot flight paths + derive flight numbers
from the actual trajectories."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from disc_model import DiscParams
import flight_sim as fs

PRESETS = {
    "Putter":   dict(D_out_mm=212, rim_width_mm=11.5, rim_shoulder_mm=16,  parting_line_mm=9,  dome_mm=5.5, plate_thick_mm=2.4, nose_radius_mm=2.0, density_gcc=.98),
    "Midrange": dict(D_out_mm=212, rim_width_mm=13.5, rim_shoulder_mm=14.5, parting_line_mm=10, dome_mm=3.5, plate_thick_mm=2.2, nose_radius_mm=2.2, density_gcc=.97),
    "Fairway":  dict(D_out_mm=212, rim_width_mm=17.5, rim_shoulder_mm=13.5, parting_line_mm=9,  dome_mm=2.0, plate_thick_mm=1.9, nose_radius_mm=1.8, density_gcc=.95),
    "Distance": dict(D_out_mm=211, rim_width_mm=22,   rim_shoulder_mm=13.5, parting_line_mm=8,  dome_mm=1.0, plate_thick_mm=1.6, nose_radius_mm=1.7, density_gcc=.90),
}
SPEEDS = {"Putter": 14, "Midrange": 17, "Fairway": 22, "Distance": 27}
COLORS = {"Putter": "#2563eb", "Midrange": "#16a34a", "Fairway": "#d97706", "Distance": "#dc2626"}


def derived_flight_numbers(res):
    """Turn a simulated trajectory into speed/glide/turn/fade-style descriptors."""
    x, y, z = res["x"], res["y"], res["z"]
    n = len(x)
    dist_ft = res["distance_ft"]
    # speed proxy: how far it went (normalized to a ~14-class scale)
    speed = np.clip(dist_ft / 34.0, 1, 14)
    # glide proxy: carry per unit drop -> hang time relative to distance
    glide = np.clip(res["flight_time_s"] / (dist_ft / 100.0) * 1.5, 1, 7)
    # turn: max RIGHT excursion (-y) during first 60% of flight (high-speed phase)
    i60 = int(0.6 * (n - 1))
    turn_right = -min(0.0, y[:i60].min())          # metres to the right
    turn = np.clip(-turn_right * 1.1, -5, 1)
    # fade: net LEFT pull from the turning point to landing (low-speed phase)
    fade_left = y[-1] - y[i60:].min()
    fade = np.clip(max(0.0, y[-1]) * 0.9 + max(0.0, fade_left) * 0.4, 0, 5)
    return dict(speed=round(float(speed), 1), glide=round(float(glide), 1),
                turn=round(float(turn), 1), fade=round(float(fade), 1))


fig, (axtop, axside) = plt.subplots(1, 2, figsize=(13, 6.2), dpi=130)
rows = []
for name, pd in PRESETS.items():
    p = DiscParams(**pd)
    res = fs.simulate(p, v0=SPEEDS[name], launch_deg=9)
    c = COLORS[name]
    # top-down: downrange (x) vs lateral (y). RHBH: +y left, -y right.
    axtop.plot(res["y"], res["x"], color=c, lw=2.4, label=name)
    axtop.plot(res["y"][-1], res["x"][-1], "o", color=c, ms=7)
    # side: downrange vs altitude
    axside.plot(res["x"], res["z"], color=c, lw=2.4, label=name)
    fn = derived_flight_numbers(res)
    rows.append((name, res["distance_ft"], res["apex_m"], res["flight_time_s"], fn, res["coeffs"]))

axtop.axvline(0, color="#bbb", ls=":", lw=1)
axtop.set_xlabel("lateral offset (m)   —   tee view:   ← left        right →")
axtop.set_ylabel("downrange (m)")
axtop.set_title("Top-down flight path (flat RHBH throw)")
axtop.invert_xaxis()   # so right is on the right visually
axtop.grid(alpha=.15); axtop.legend()

axside.set_xlabel("downrange (m)"); axside.set_ylabel("altitude (m)")
axside.set_title("Altitude profile")
axside.axhline(0, color="#8a5a3a", lw=1)
axside.grid(alpha=.15); axside.legend()
fig.tight_layout()
fig.savefig("flight_paths.png", dpi=130)

print(f"{'disc':9} {'dist(ft)':>8} {'apex(m)':>7} {'t(s)':>5}   flight (S/G/T/F)     CD0   CL0    CM0")
for name, dist, apex, t, fn, ac in rows:
    print(f"{name:9} {dist:8.0f} {apex:7.1f} {t:5.2f}   "
          f"{fn['speed']:>4}/{fn['glide']:>3}/{fn['turn']:>4}/{fn['fade']:>3}   "
          f"{ac.CD0:.3f} {ac.CL0:.3f} {ac.CM0:+.3f}")
print("saved flight_paths.png")
