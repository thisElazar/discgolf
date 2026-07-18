"""
6-DOF flying-disc flight simulator  --  physics tier
====================================================

Equations of motion and aerodynamic coefficient FORM follow Hummel (2003),
"Frisbee Flight Simulation and Throw Biomechanics" (UC Davis M.S. thesis):

    L = 1/2 rho v^2 A (CL0 + CLa*alpha)
    D = 1/2 rho v^2 A (CD0 + CDa*(alpha-alpha0)^2)
    M = 1/2 rho v^2 A d (CM0 + CMa*alpha + CMq*q)      (pitch)
    R = 1/2 rho v^2 A d (CRr*r + CRp*p)                (roll)
    N = 1/2 rho v^2 A d (CNr*r)                        (spin-down)

The disc is treated as an axisymmetric spinning rigid body; gyroscopic
precession (what turns a pitching moment into roll = "turn"/"fade")
emerges from the full Euler equations, not from any hand-coded rule.

WHAT IS GROUNDED vs CALIBRATED (see docs/METHODS.md):
  * EOM + coefficient functional forms + moment coeffs  -> Hummel (2003)  [grounded]
  * Disc mass & moment of inertia                        -> our geometry model [grounded]
  * Lift/drag magnitudes for the disc-golf regime and the
    geometry->coefficient scaling                        -> CALIBRATED / heuristic,
    pending a CFD sweep. Signs & trends are physical; absolute
    numbers are first-order.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from disc_model import DiscParams, mass_properties

RHO = 1.225          # air density kg/m^3
G = 9.81
DEG = np.pi / 180.0

# Low-speed fade model (CALIBRATED, not from Hummel): as airspeed drops the
# disc behaves more overstable and fades to the RHBH-left. This stands in for
# advance-ratio-dependent moments that need CFD/wind-tunnel data to pin down.
LOWSPEED_FADE_VREF = 14.0   # m/s, fade grows below this
LOWSPEED_FADE_GAIN = 0.060  # CM0 shift at zero speed


# --------------------------------------------------------------------------
# geometry -> aerodynamic coefficients
# --------------------------------------------------------------------------
@dataclass
class AeroCoeffs:
    CL0: float; CLa: float
    CD0: float; CDa: float; alpha0: float
    CM0: float; CMa: float; CMq: float
    CRr: float; CRp: float; CNr: float


def aero_from_geometry(p: DiscParams, mp: dict | None = None) -> AeroCoeffs:
    """Map parametric geometry -> aero coefficients.
    Anchored so a wide-rim low-dome driver is fast/low-drag and a
    narrow-rim high-dome putter is slow/high-lift/high-drag."""
    if mp is None:
        mp = mass_properties(p)
    rim = p.rim_width_mm
    dome = p.dome_mm
    h = mp["total_height_mm"]
    plr = p.parting_line_mm / max(p.rim_shoulder_mm, 1e-9)   # parting-line ratio

    # --- drag: falls with rim width (streamlining), rises with height/bluffness
    CD0 = float(np.clip(0.100 - 0.0030 * (rim - 9.0) + 0.0016 * (h - 15.0), 0.035, 0.140))
    CDa = 0.69                                              # Hummel

    # --- lift: rises with dome (camber) and overall height
    CL0 = float(np.clip(0.115 + 0.0130 * dome + 0.0040 * (h - 15.0), 0.05, 0.35))
    CLa = 2.0                                               # ~Hummel (1.9)
    alpha0 = -4.0 * DEG

    # --- pitching moment offset CM0 sets high-speed turn / low-speed fade.
    # "Understability" u rises with dome (camber), parting-line height, and rim
    # width. Baseline is biased overstable so neutral geometry fades left (RHBH);
    # only domey/high-parting shapes cross into high-speed turn.
    u = (plr - 0.60) * 1.5 + dome / 8.0 + (rim - 14.0) / 40.0
    CM0 = float(np.clip(-0.028 + 0.060 * u, -0.060, 0.055))
    CMa = 0.43                                              # Hummel (pitch stiffness)
    CMq = -5.0e-3                                           # Hummel
    CRr = 1.4e-2                                            # Hummel
    CRp = -5.5e-3                                           # Hummel
    CNr = -3.4e-5                                           # scaled from Hummel (spin decay)
    return AeroCoeffs(CL0, CLa, CD0, CDa, alpha0, CM0, CMa, CMq, CRr, CRp, CNr)


# --------------------------------------------------------------------------
# rigid-body helpers
# --------------------------------------------------------------------------
def rot_body_to_inertial(phi, theta, psi):
    """3-2-1 (yaw psi, pitch theta, roll phi) body->inertial rotation."""
    cph, sph = np.cos(phi), np.sin(phi)
    cth, sth = np.cos(theta), np.sin(theta)
    cps, sps = np.cos(psi), np.sin(psi)
    return np.array([
        [cth * cps, sph * sth * cps - cph * sps, cph * sth * cps + sph * sps],
        [cth * sps, sph * sth * sps + cph * cps, cph * sth * sps - sph * cps],
        [-sth,      sph * cth,                   cph * cth],
    ])


@dataclass
class DiscPhysical:
    m: float      # kg
    d: float      # m
    A: float      # m^2
    Ia: float     # axial MOI (spin axis) kg m^2
    Id: float     # diametral MOI kg m^2


def physical_from_geometry(p: DiscParams, mp: dict | None = None) -> DiscPhysical:
    if mp is None:
        mp = mass_properties(p)
    d = p.D_out_mm / 1000.0
    A = np.pi * (d / 2.0) ** 2
    Ia = mp["Iz_gcm2"] * 1e-7                 # g*cm^2 -> kg*m^2
    Id = 0.5 * Ia                             # thin-disc perpendicular-axis approx
    return DiscPhysical(mp["mass_g"] / 1000.0, d, A, Ia, Id)


# --------------------------------------------------------------------------
# equations of motion
# --------------------------------------------------------------------------
def _deriv(s, phys: DiscPhysical, ac: AeroCoeffs):
    x, y, z, vx, vy, vz, phi, theta, psi, p, q, r = s
    V = np.array([vx, vy, vz])
    vmag = np.linalg.norm(V)
    if vmag < 1e-3:
        return np.zeros(12)
    Vhat = V / vmag

    R = rot_body_to_inertial(phi, theta, psi)
    zb = R @ np.array([0.0, 0.0, 1.0])        # disc normal (inertial)

    # angle of attack: angle between velocity and disc plane (+ = underside into wind)
    alpha = -np.arcsin(np.clip(np.dot(Vhat, zb), -1.0, 1.0))

    qdyn = 0.5 * RHO * vmag * vmag
    Lf = qdyn * phys.A * (ac.CL0 + ac.CLa * alpha)
    Df = qdyn * phys.A * (ac.CD0 + ac.CDa * (alpha - ac.alpha0) ** 2)

    # force directions (inertial)
    drag_dir = -Vhat
    perp = zb - np.dot(zb, Vhat) * Vhat       # zb component perpendicular to V
    lift_dir = perp / (np.linalg.norm(perp) + 1e-12)
    F = Lf * lift_dir + Df * drag_dir + np.array([0.0, 0.0, -phys.m * G])
    acc = F / phys.m

    # aero moments (magnitudes), nondim rates by d/(2v)
    k = qdyn * phys.A * phys.d
    fade = LOWSPEED_FADE_GAIN * max(0.0, (LOWSPEED_FADE_VREF - vmag) / LOWSPEED_FADE_VREF)
    CM0_eff = ac.CM0 - fade                                                 # more overstable when slow
    Mp = k * (CM0_eff + ac.CMa * alpha + ac.CMq * q * phys.d / (2 * vmag))  # pitch
    Rl = k * (ac.CRr * r + ac.CRp * p) * phys.d / (2 * vmag)                # roll
    Nz = k * (ac.CNr * r) * phys.d / (2 * vmag)                            # spin-down

    # moment axes in inertial frame
    pitch_axis = np.cross(zb, Vhat)
    pitch_axis /= (np.linalg.norm(pitch_axis) + 1e-12)
    roll_axis = Vhat - np.dot(Vhat, zb) * zb
    roll_axis /= (np.linalg.norm(roll_axis) + 1e-12)
    tau_inertial = Mp * pitch_axis + Rl * roll_axis + Nz * zb
    tau_b = R.T @ tau_inertial

    # symmetric-body Euler eqns (Ixx=Iyy=Id, Izz=Ia)
    Id_, Ia_ = phys.Id, phys.Ia
    pdot = (tau_b[0] - (Ia_ - Id_) * q * r) / Id_
    qdot = (tau_b[1] - (Id_ - Ia_) * r * p) / Id_
    rdot = tau_b[2] / Ia_

    # Euler-angle kinematics
    phidot = p + (q * np.sin(phi) + r * np.cos(phi)) * np.tan(theta)
    thetadot = q * np.cos(phi) - r * np.sin(phi)
    psidot = (q * np.sin(phi) + r * np.cos(phi)) / np.cos(theta)

    return np.array([vx, vy, vz, acc[0], acc[1], acc[2],
                     phidot, thetadot, psidot, pdot, qdot, rdot])


def simulate(p: DiscParams, v0=24.0, spin_rev=None, launch_deg=8.0,
             hyzer_deg=0.0, nose_deg=0.0, h0=1.4, dt=5e-4, tmax=15.0,
             handed="RHBH"):
    """Throw the disc. Returns trajectory dict.
    v0: launch speed m/s. launch_deg: upward launch angle. hyzer_deg: bank
    (positive = hyzer for RHBH -> left-tilt). nose_deg: nose up(+)/down(-)."""
    mp = mass_properties(p)
    phys = physical_from_geometry(p, mp)
    ac = aero_from_geometry(p, mp)

    if spin_rev is None:
        spin_rev = v0 / (np.pi * phys.d)      # ~rolling-contact spin estimate
    spin = 2 * np.pi * spin_rev
    # RHBH spins clockwise seen from above -> spin angular velocity points DOWN
    r0 = -spin if handed == "RHBH" else spin

    la = launch_deg * DEG
    # z-up frame: +theta in the 3-2-1 matrix pitches nose DOWN, so negate.
    # Disc plane attitude = climb angle + nose angle -> initial alpha = nose_deg.
    theta0 = -(launch_deg + nose_deg) * DEG
    phi0 = -hyzer_deg * DEG if handed == "RHBH" else hyzer_deg * DEG

    s = np.array([0, 0, h0,
                  v0 * np.cos(la), 0.0, v0 * np.sin(la),
                  phi0, theta0, 0.0,
                  0.0, 0.0, r0], dtype=float)

    traj = [s.copy()]
    t = 0.0
    while t < tmax and s[2] > 0.0:
        k1 = _deriv(s, phys, ac)
        k2 = _deriv(s + 0.5 * dt * k1, phys, ac)
        k3 = _deriv(s + 0.5 * dt * k2, phys, ac)
        k4 = _deriv(s + dt * k3, phys, ac)
        s = s + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
        traj.append(s.copy())
    T = np.array(traj)

    x, y, z = T[:, 0], T[:, 1], T[:, 2]
    return dict(
        t=np.arange(len(T)) * dt, x=x, y=y, z=z, state=T,
        distance_m=float(x[-1]),
        distance_ft=float(x[-1] / 0.3048),
        lateral_max_m=float(y[np.argmax(np.abs(y))]),
        apex_m=float(z.max()),
        flight_time_s=float((len(T) - 1) * dt),
        coeffs=ac, phys=phys, mp=mp,
    )


if __name__ == "__main__":
    for name, kw in [("driver-ish", {}), ]:
        p = DiscParams(rim_width_mm=21, rim_shoulder_mm=13.5, parting_line_mm=8,
                       dome_mm=1.0, plate_thick_mm=1.6, density_gcc=0.90, D_out_mm=211)
        res = simulate(p, v0=27, launch_deg=10)
        print(f"{name}: dist {res['distance_ft']:.0f} ft  apex {res['apex_m']:.1f} m  "
              f"lateral {res['lateral_max_m']:+.1f} m  t {res['flight_time_s']:.2f} s  "
              f"CD0 {res['coeffs'].CD0:.3f} CL0 {res['coeffs'].CL0:.3f} CM0 {res['coeffs'].CM0:+.3f}")
