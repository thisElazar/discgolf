"""
Parametric disc golf disc model  --  v1 foundation
==================================================

A disc is modelled as a *solid of revolution*: we define a radial half
cross-section (a top surface z_top(r) and a material-bottom surface
z_bot(r), for r from the spin axis out to the rim edge) and revolve it
360 deg about the vertical axis.

From that single cross-section we derive everything downstream:
  - cross-sectional area, enclosed volume (Pappus / ring integration)
  - mass (volume x material density)
  - polar moment of inertia Iz about the spin axis (gyroscopic stability)
  - centre-of-gravity height
  - full PDGA legality check
  - a *placeholder* flight-number estimate (to be replaced by a
    CFD-trained surrogate later)

All linear dimensions are in MILLIMETRES internally. Mass in grams,
density in g/cm^3.

This is deliberately a readable reference implementation. The same
parameter set and profile logic is mirrored in the GhPython version so
Rhino and this script stay in agreement.
"""

from __future__ import annotations
from dataclasses import dataclass, asdict
import numpy as np

# --------------------------------------------------------------------------
# PDGA technical standards (source: pdga.com/technical-standards/guidelines)
# --------------------------------------------------------------------------
PDGA = dict(
    diameter_min_cm=21.0,
    diameter_max_cm=30.0,
    rim_depth_pct_min=0.05,        # >=5% of outside diameter
    rim_depth_pct_max=0.12,        # <=12% of outside diameter
    rim_depth_abs_min_cm=1.10,     # and at least 1.1 cm
    rim_width_max_cm=2.60,
    inner_rim_dia_min_cm=15.80,
    plate_thickness_max_cm=0.50,   # standard flight plate
    weight_per_cm_max=8.30,        # g per cm of outside diameter
    weight_abs_max_g=200.0,
    edge_radius_gauge_mm=1.60,     # 1/16" leading-edge radius gauge
)


@dataclass
class DiscParams:
    """The parametric 'sliders'. Defaults describe a stable putter."""
    D_out_mm:        float = 212.0   # outside diameter
    rim_width_mm:    float = 11.5    # radial width of the thick rim band
    rim_shoulder_mm: float = 14.5    # z-height of the rim top at its inner shoulder
    parting_line_mm: float = 8.5     # z-height of the outer edge (widest point)
    dome_mm:         float = 4.5     # extra dome height at centre, above the shoulder
    plate_thick_mm:  float = 1.9     # flight-plate thickness
    nose_radius_mm:  float = 2.0     # leading-edge (nose) radius
    top_rim_power:   float = 2.0     # >1 flatter shoulder, curls near edge (ellipse-like)
    dome_power:      float = 2.0     # 2 = parabolic dome; higher = flatter puddle-top
    density_gcc:     float = 0.905   # ~polypropylene; disc plastics ~0.90-0.95

    # ---- derived radii ----
    @property
    def R(self) -> float:
        return self.D_out_mm / 2.0

    @property
    def R_in(self) -> float:
        return self.R - self.rim_width_mm


# --------------------------------------------------------------------------
# Profile: the two surfaces of the cross-section as functions of radius r
# --------------------------------------------------------------------------
def _z_top(p: DiscParams, r: np.ndarray) -> np.ndarray:
    """Top surface height (mm) vs radius (mm)."""
    R, R_in = p.R, p.R_in
    z = np.empty_like(r, dtype=float)

    # Flight-plate dome:  shoulder -> shoulder+dome at centre
    plate = r <= R_in
    x = np.clip(r[plate] / max(R_in, 1e-9), 0.0, 1.0)
    z[plate] = p.rim_shoulder_mm + p.dome_mm * (1.0 - x ** p.dome_power)

    # Rim top:  shoulder at R_in  ->  parting line at R  (convex, ellipse-like)
    rim = ~plate
    u = np.clip((r[rim] - R_in) / max(p.rim_width_mm, 1e-9), 0.0, 1.0)
    drop = (p.rim_shoulder_mm - p.parting_line_mm)
    z[rim] = p.parting_line_mm + drop * (1.0 - u) ** (1.0 / p.top_rim_power)
    return z


def _z_bot(p: DiscParams, r: np.ndarray) -> np.ndarray:
    """Material bottom height (mm) vs radius (mm).

    Plate region: underside of the (thin) flight plate  -> air below it.
    Rim region:   solid rim sits on the ground (z=0) and rolls up at the
                  leading edge to meet the top at the parting line.
    """
    R, R_in = p.R, p.R_in
    z = np.empty_like(r, dtype=float)

    plate = r <= R_in
    z[plate] = _z_top(p, r[plate]) - p.plate_thick_mm

    rim = ~plate
    rr = r[rim]
    z[rim] = 0.0
    # leading-edge rollover: last nose_radius mm curve up to the parting line
    le = rr >= (R - p.nose_radius_mm)
    if np.any(le):
        # quarter-circle-ish rise so bottom meets top (thickness -> ~0 at edge)
        t = np.clip((rr[le] - (R - p.nose_radius_mm)) / max(p.nose_radius_mm, 1e-9), 0, 1)
        z_rim = z[rim]
        z_rim[le] = p.parting_line_mm * (1.0 - np.sqrt(np.clip(1.0 - t ** 2, 0, 1)))
        z[rim] = z_rim
    return z


# --------------------------------------------------------------------------
# Geometry / mass properties by ring integration over radius
# --------------------------------------------------------------------------
def mass_properties(p: DiscParams, n: int = 4000) -> dict:
    r = np.linspace(0.0, p.R, n)                 # mm
    top = _z_top(p, r)
    bot = _z_bot(p, r)
    h = np.clip(top - bot, 0.0, None)            # material thickness column (mm)
    zmid = 0.5 * (top + bot)

    # ring integration:  dV = 2*pi*r*h*dr
    dr = r[1] - r[0]
    ring_V = 2.0 * np.pi * r * h * dr            # mm^3 per ring
    V_mm3 = ring_V.sum()
    V_cm3 = V_mm3 / 1000.0
    mass_g = V_cm3 * p.density_gcc

    # polar moment of inertia about spin (z) axis:  Iz = integral r^2 dm
    #   dm = rho * dV ;  work in cm to get g*cm^2
    r_cm = r / 10.0
    h_cm = h / 10.0
    ring_dm = p.density_gcc * 2.0 * np.pi * r_cm * h_cm * (dr / 10.0)   # g
    Iz_gcm2 = np.sum(ring_dm * r_cm ** 2)        # g*cm^2
    Iz_kgm2 = Iz_gcm2 * 1e-7

    # CG height (material-weighted)
    z_cg_mm = float(np.sum(ring_V * zmid) / V_mm3) if V_mm3 > 0 else 0.0

    # rim mass fraction (rim region vs total) -- drives 'feel' & stability
    rim_mask = r > p.R_in
    rim_V = ring_V[rim_mask].sum()
    rim_mass_frac = float(rim_V / V_mm3) if V_mm3 > 0 else 0.0

    # cavity (underside inner rim) depth = plate underside height at R_in
    cavity_depth_mm = float(_z_top(p, np.array([p.R_in]))[0] - p.plate_thick_mm)

    total_height_mm = float(top.max())

    return dict(
        volume_cm3=V_cm3,
        mass_g=mass_g,
        Iz_gcm2=Iz_gcm2,
        Iz_kgm2=Iz_kgm2,
        cg_height_mm=z_cg_mm,
        rim_mass_fraction=rim_mass_frac,
        cavity_depth_mm=cavity_depth_mm,
        total_height_mm=total_height_mm,
        inner_rim_dia_mm=2.0 * p.R_in,
    )


# --------------------------------------------------------------------------
# PDGA legality
# --------------------------------------------------------------------------
def pdga_check(p: DiscParams, mp: dict | None = None) -> dict:
    if mp is None:
        mp = mass_properties(p)
    D_cm = p.D_out_mm / 10.0
    cavity_cm = mp["cavity_depth_mm"] / 10.0

    checks = {}

    def add(name, ok, value, limit):
        checks[name] = dict(pass_=bool(ok), value=round(float(value), 3), limit=limit)

    add("outside_diameter_cm", PDGA["diameter_min_cm"] <= D_cm <= PDGA["diameter_max_cm"],
        D_cm, f"{PDGA['diameter_min_cm']}-{PDGA['diameter_max_cm']}")

    rim_ok = (PDGA["rim_depth_pct_min"] * D_cm <= cavity_cm <= PDGA["rim_depth_pct_max"] * D_cm) \
        and (cavity_cm >= PDGA["rim_depth_abs_min_cm"])
    add("rim_depth_cm", rim_ok, cavity_cm,
        f"{PDGA['rim_depth_pct_min']*D_cm:.2f}-{PDGA['rim_depth_pct_max']*D_cm:.2f} & >= {PDGA['rim_depth_abs_min_cm']}")

    add("rim_width_cm", (p.rim_width_mm / 10.0) <= PDGA["rim_width_max_cm"],
        p.rim_width_mm / 10.0, f"<= {PDGA['rim_width_max_cm']}")

    add("inner_rim_dia_cm", (mp["inner_rim_dia_mm"] / 10.0) >= PDGA["inner_rim_dia_min_cm"],
        mp["inner_rim_dia_mm"] / 10.0, f">= {PDGA['inner_rim_dia_min_cm']}")

    add("plate_thickness_cm", (p.plate_thick_mm / 10.0) <= PDGA["plate_thickness_max_cm"],
        p.plate_thick_mm / 10.0, f"<= {PDGA['plate_thickness_max_cm']}")

    wpc_ok = mp["mass_g"] <= PDGA["weight_per_cm_max"] * D_cm
    add("weight_per_cm_g", wpc_ok, mp["mass_g"] / D_cm, f"<= {PDGA['weight_per_cm_max']}")

    add("weight_abs_g", mp["mass_g"] <= PDGA["weight_abs_max_g"],
        mp["mass_g"], f"<= {PDGA['weight_abs_max_g']}")

    add("edge_radius_mm", p.nose_radius_mm >= PDGA["edge_radius_gauge_mm"],
        p.nose_radius_mm, f">= {PDGA['edge_radius_gauge_mm']}")

    checks["ALL_PASS"] = all(c["pass_"] for c in checks.values())
    return checks


# --------------------------------------------------------------------------
# Flight estimate — calibrated linear models when data/flight_fit.json exists
# (generated by scripts/fit_flight_numbers.py; fit to manufacturer flight
# numbers of 16 pro discs, LOO-validated), else the original rough heuristic.
# Both are interim: the CFD-trained surrogate replaces them.
# --------------------------------------------------------------------------
def _load_flight_fit() -> dict | None:
    import json
    from pathlib import Path
    path = Path(__file__).resolve().parent / "data" / "flight_fit.json"
    if not path.exists():
        return None
    return json.loads(path.read_text())["fit"]


_FLIGHT_FIT = _load_flight_fit()


def flight_estimate(p: DiscParams, mp: dict | None = None) -> dict:
    if mp is None:
        mp = mass_properties(p)
    if _FLIGHT_FIT is not None:
        feat = dict(rw=p.rim_width_mm, dia=p.D_out_mm, dome=p.dome_mm,
                    h=mp["total_height_mm"], cav=mp["cavity_depth_mm"],
                    pr=p.parting_line_mm / max(p.rim_shoulder_mm, 1e-9))
        out = {}
        for k, m in _FLIGHT_FIT.items():
            v = m["intercept"] + sum(c * feat[f] for f, c in m["features"].items())
            out[k] = round(float(np.clip(v, *m["clamp"])), 1)
        return out
    # fallback: rough uncalibrated heuristic
    speed = np.clip(1.0 + (p.rim_width_mm - 6.0) / 1.65, 1, 14)
    glide = np.clip(2.0 + p.dome_mm / 1.4 + (mp["total_height_mm"] - 18) / 6.0, 1, 7)
    understab = (p.parting_line_mm / max(p.rim_shoulder_mm, 1e-9)) + p.dome_mm / 12.0
    turn = np.clip(1.0 - (understab - 0.55) * 8.0, -5, 1)
    fade = np.clip(0.5 + (p.rim_width_mm - 10) / 6.0 - p.dome_mm / 10.0, 0, 5)
    return dict(speed=round(float(speed), 1), glide=round(float(glide), 1),
                turn=round(float(turn), 1), fade=round(float(fade), 1))


def summary(p: DiscParams) -> dict:
    mp = mass_properties(p)
    return dict(params=asdict(p), mass_properties=mp,
                pdga=pdga_check(p, mp), flight_estimate=flight_estimate(p, mp))


if __name__ == "__main__":
    import json
    p = DiscParams()
    s = summary(p)
    print(json.dumps(s, indent=2, default=float))
