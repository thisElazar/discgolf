"""
Disc -> watertight STL for CFD.
========================================================================
Reproduces the app's `paramsFromCertified` (disc_designer.html) and the
`disc_model.py` solid-of-revolution profile, then revolves the CLOSED
radial cross-section about the spin (z) axis into a watertight triangle
mesh written as a binary STL for OpenFOAM's snappyHexMesh.

The CFD disc is therefore geometrically identical to the disc the
designer app models, so validated CFD coefficients map straight back
onto the same parametric geometry.

Axis convention for CFD:  spin axis = +z (vertical), disc lies in the
x-y plane, freestream travels in +x.  Angle of attack is imposed later
by rotating the freestream vector, not the mesh.

Units: the profile math is in millimetres; the STL is written in METRES
(OpenFOAM SI), i.e. all coordinates scaled by 1e-3.
"""
from __future__ import annotations
import argparse, json, math, struct
from dataclasses import dataclass
import numpy as np

# --- per-class estimates for what the PDGA doesn't publish (from the app) ----
CLASS_EST = {
    "Putter":   {"plate": 3.2, "nose": 3.0},
    "Midrange": {"plate": 2.6, "nose": 2.5},
    "Fairway":  {"plate": 2.0, "nose": 2.0},
    "Distance": {"plate": 1.7, "nose": 1.7},
}
# map the DB's free-text class to an estimate bucket
CLS_ALIAS = {"Putter": "Putter", "Putt & Approach": "Putter", "Midrange": "Midrange",
             "Mid-Range": "Midrange", "Fairway Driver": "Fairway", "Control Driver": "Fairway",
             "Distance Driver": "Distance", "Driver": "Distance"}


@dataclass
class DiscParams:
    D_out_mm: float; rim_width_mm: float; rim_shoulder_mm: float
    parting_line_mm: float; dome_mm: float; plate_thick_mm: float
    nose_radius_mm: float; density_gcc: float = 0.95
    top_rim_power: float = 2.0; dome_power: float = 2.0

    @property
    def R(self): return self.D_out_mm / 2.0
    @property
    def R_in(self): return self.R - self.rim_width_mm


def params_from_certified(d: dict) -> DiscParams:
    """Port of paramsFromCertified() in disc_designer.html."""
    est = CLASS_EST[CLS_ALIAS.get(d["cls"], "Putter")]
    shoulder = d["rd"] + est["plate"]
    dome = min(max(d["h"] - shoulder, 0.2), 8.0)
    fl = d.get("flight")
    stab = (fl[3] + fl[2]) if fl else 1          # fade + turn
    parting = min(max(shoulder * (0.60 - 0.04 * stab), 5.0), shoulder - 2.0)
    return DiscParams(D_out_mm=d["dia"], rim_width_mm=d["rw"],
                      rim_shoulder_mm=round(shoulder, 1),
                      parting_line_mm=round(parting, 1), dome_mm=round(dome, 1),
                      plate_thick_mm=est["plate"], nose_radius_mm=est["nose"])


def z_top(p: DiscParams, r: np.ndarray) -> np.ndarray:
    R, R_in = p.R, p.R_in
    z = np.empty_like(r, dtype=float)
    plate = r <= R_in
    x = np.clip(r[plate] / max(R_in, 1e-9), 0.0, 1.0)
    z[plate] = p.rim_shoulder_mm + p.dome_mm * (1.0 - x ** p.dome_power)
    rim = ~plate
    u = np.clip((r[rim] - R_in) / max(p.rim_width_mm, 1e-9), 0.0, 1.0)
    drop = p.rim_shoulder_mm - p.parting_line_mm
    z[rim] = p.parting_line_mm + drop * (1.0 - u) ** (1.0 / p.top_rim_power)
    return z


def z_bot(p: DiscParams, r: np.ndarray) -> np.ndarray:
    R, R_in = p.R, p.R_in
    z = np.empty_like(r, dtype=float)
    plate = r <= R_in
    z[plate] = z_top(p, r[plate]) - p.plate_thick_mm
    rim = ~plate
    rr = r[rim]
    zb = np.zeros_like(rr)
    le = rr >= (R - p.nose_radius_mm)
    if np.any(le):
        t = np.clip((rr[le] - (R - p.nose_radius_mm)) / max(p.nose_radius_mm, 1e-9), 0, 1)
        zb[le] = p.parting_line_mm * (1.0 - np.sqrt(np.clip(1.0 - t ** 2, 0, 1)))
    z[rim] = zb
    return z


def cross_section_polyline(p: DiscParams, n_top=240, n_bot=200):
    """Closed 2-D profile (r,z) in mm, counter-clockwise, starting at the
    top apex (r=0). Order: top apex -> outer edge along top surface ->
    back along the underside (nose rollover, flat rim bottom, inner cavity
    wall, plate underside) -> bottom-centre. r=0 endpoints sit on the axis
    so the revolve closes to cone tips (watertight, no axis cap needed)."""
    R, R_in = p.R, p.R_in
    # top surface: 0 -> R (denser toward the rim where curvature is high)
    rt = R * np.linspace(0.0, 1.0, n_top) ** 1.3
    top = np.column_stack([rt, z_top(p, rt)])
    # underside: R -> 0.  sample rim band finely (nose + inner wall), plate coarser
    r_rim = np.linspace(R, R_in, max(int(n_bot * p.rim_width_mm / R) + 8, 24))
    r_plate = np.linspace(R_in, 0.0, n_bot)[1:]     # drop dup at R_in
    rb = np.concatenate([r_rim, r_plate])
    bot = np.column_stack([rb, z_bot(p, rb)])
    poly = np.vstack([top, bot])
    # de-duplicate consecutive identical points (edge apex where top meets bottom)
    keep = np.ones(len(poly), bool)
    d2 = np.sum(np.diff(poly, axis=0) ** 2, axis=1)
    keep[1:] = d2 > 1e-12
    return poly[keep]


def revolve_to_stl(poly_mm: np.ndarray, path: str, n_theta=200, scale=1e-3):
    """Revolve the closed (r,z) profile about z into a watertight binary STL.
    Rings at r~0 collapse to a single axis point (cone tip)."""
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    ct, st = np.cos(theta), np.sin(theta)
    P = poly_mm.shape[0]
    # vertex grid V[i, j] = profile point i revolved to angle j  -> (x,y,z) metres
    r = poly_mm[:, 0] * scale
    z = poly_mm[:, 1] * scale
    X = np.outer(r, ct)          # (P, n_theta)
    Y = np.outer(r, st)
    Z = np.repeat(z[:, None], n_theta, axis=1)

    tris = []
    axis_tol = 1e-9
    for i in range(P - 1):
        r0, r1 = r[i], r[i + 1]
        for j in range(n_theta):
            k = (j + 1) % n_theta
            a = (X[i, j], Y[i, j], Z[i, j])
            b = (X[i, k], Y[i, k], Z[i, k])
            c = (X[i + 1, k], Y[i + 1, k], Z[i + 1, k])
            dd = (X[i + 1, j], Y[i + 1, j], Z[i + 1, j])
            if r0 <= axis_tol:            # top ring collapses to apex a: fan to ring i+1
                tris.append((a, dd, c))
            elif r1 <= axis_tol:          # bottom ring collapses to apex dd: fan from ring i
                tris.append((a, b, dd))
            else:
                tris.append((a, b, c))
                tris.append((a, c, dd))

    # winding: ensure outward normals by checking against profile orientation later;
    # snappyHexMesh only needs consistent orientation + watertight, which revolve gives.
    with open(path, "wb") as f:
        f.write(b"\0" * 80)
        f.write(struct.pack("<I", len(tris)))
        for (a, b, c) in tris:
            ux, uy, uz = (b[0]-a[0], b[1]-a[1], b[2]-a[2])
            vx, vy, vz = (c[0]-a[0], c[1]-a[1], c[2]-a[2])
            nx, ny, nz = (uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx)
            nrm = math.sqrt(nx*nx+ny*ny+nz*nz) or 1.0
            f.write(struct.pack("<3f", nx/nrm, ny/nrm, nz/nrm))
            for v in (a, b, c):
                f.write(struct.pack("<3f", *v))
            f.write(b"\0\0")
    return len(tris)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--disc", default="Aviar", help="disc name to find in pro_discs.json")
    ap.add_argument("--pro-discs", default="../data/pro_discs.json")
    ap.add_argument("--out", default="case/constant/triSurface/disc.stl")
    ap.add_argument("--n-theta", type=int, default=200)
    args = ap.parse_args()

    discs = json.load(open(args.pro_discs))["discs"]
    d = next(x for x in discs if x["name"].lower() == args.disc.lower())
    p = params_from_certified(d)
    poly = cross_section_polyline(p)
    import os
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    ntri = revolve_to_stl(poly, args.out, n_theta=args.n_theta)

    R = p.R / 1000.0
    # volume centroid height (z c.g., uniform density) by ring integration in mm
    rr = np.linspace(0.0, p.R, 4000)
    h = np.clip(z_top(p, rr) - z_bot(p, rr), 0.0, None)
    zmid = 0.5 * (z_top(p, rr) + z_bot(p, rr))
    ringV = 2.0 * np.pi * rr * h
    cg_z_m = float(np.sum(ringV * zmid) / np.sum(ringV)) / 1000.0

    case_params = {
        "disc": d["name"], "diameter_m": round(2 * R, 6),
        "planform_area_m2": round(math.pi * R * R, 8),
        "cofr_m": [0.0, 0.0, round(cg_z_m, 6)],
        "Uinf": 26.75, "nu": 1.5e-5, "rhoInf": 1.225,
        "Re": 378000, "turb_intensity": 0.005, "turb_length_m": round(0.07 * 2 * R, 5),
    }
    case_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(args.out))))
    cp_path = os.path.join(os.path.dirname(case_dir), "case_params.json")
    with open(cp_path, "w") as f:
        json.dump(case_params, f, indent=2)

    print(json.dumps({
        "disc": d["name"], "class": d["cls"],
        "params": {k: getattr(p, k) for k in
                   ("D_out_mm","rim_width_mm","rim_shoulder_mm","parting_line_mm",
                    "dome_mm","plate_thick_mm","nose_radius_mm")},
        "diameter_m": round(2*R, 4),
        "planform_area_m2": round(math.pi * R * R, 6),
        "cg_z_m": round(cg_z_m, 5),
        "profile_pts": int(poly.shape[0]), "triangles": ntri,
        "stl": args.out,
    }, indent=2))


if __name__ == "__main__":
    main()
