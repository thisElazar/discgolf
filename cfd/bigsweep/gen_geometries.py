"""Generate the big-sweep geometry sample.

Latin-hypercube over the disc shape parameters, filtered to PDGA-legal
shapes (shape-based rules only; weight is density-tunable so it never
constrains the aero shape). For each kept geometry writes:
    geoms/geom_XXX/constant/triSurface/disc.stl   (watertight, for meshing)
    geoms/geom_XXX/case_params.json               (per-geom area/diam/CofR)
and a manifest geometries.json listing every geometry's parameters.

Diameter is fixed at 211 mm (the dominant PDGA value; nearly constant
across real discs) to keep the surrogate's input space to 6 dims for a
first model — can be added later.
"""
import argparse, json, math, os, sys
import numpy as np
from scipy.stats import qmc

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import geometry as G   # DiscParams, z_top, z_bot, cross_section_polyline, revolve_to_stl

D_OUT = 211.0
R = D_OUT / 2.0
# param name -> (lo, hi).  parting is sampled as a fraction of shoulder.
RANGES = {
    "rim_width_mm":    (9.0, 23.0),
    "rim_shoulder_mm": (13.0, 20.0),
    "dome_mm":         (0.4, 6.0),
    "plate_thick_mm":  (1.6, 3.4),
    "nose_radius_mm":  (1.6, 3.6),
    "parting_frac":    (0.45, 0.85),
}
KEYS = list(RANGES)


def pdga_legal(p: "G.DiscParams") -> bool:
    R_in = p.R - p.rim_width_mm
    cavity = p.rim_shoulder_mm - p.plate_thick_mm         # = z_top(R_in) - plate
    inner_dia = 2.0 * R_in
    return (
        11.0 <= cavity <= 0.12 * D_OUT and                # rim depth 11mm..12%D
        p.rim_width_mm <= 26.0 and                        # rim width <=2.6cm
        inner_dia >= 158.0 and                            # inner rim dia >=15.8cm
        p.plate_thick_mm <= 5.0 and
        p.nose_radius_mm >= 1.6 and
        5.0 <= p.parting_line_mm <= p.rim_shoulder_mm - 1.5
    )


def sample_to_params(u: np.ndarray) -> "G.DiscParams":
    v = {k: RANGES[k][0] + u[i] * (RANGES[k][1] - RANGES[k][0]) for i, k in enumerate(KEYS)}
    shoulder = v["rim_shoulder_mm"]
    parting = min(max(v["parting_frac"] * shoulder, 5.0), shoulder - 1.5)
    return G.DiscParams(
        D_out_mm=D_OUT, rim_width_mm=round(v["rim_width_mm"], 2),
        rim_shoulder_mm=round(shoulder, 2), parting_line_mm=round(parting, 2),
        dome_mm=round(v["dome_mm"], 2), plate_thick_mm=round(v["plate_thick_mm"], 2),
        nose_radius_mm=round(v["nose_radius_mm"], 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=96, help="target legal geometries")
    ap.add_argument("--oversample", type=int, default=8, help="LHS pool = n*oversample")
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--start-index", type=int, default=0,
                    help="first geometry index; >0 appends to the existing manifest "
                         "(expansion batch — use a fresh --seed so samples differ)")
    ap.add_argument("--n-theta", type=int, default=480)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "geoms"))
    args = ap.parse_args()

    pool = qmc.LatinHypercube(d=len(KEYS), seed=args.seed).random(args.n * args.oversample)
    kept = []
    for u in pool:
        p = sample_to_params(u)
        if pdga_legal(p):
            kept.append(p)
        if len(kept) >= args.n:
            break

    os.makedirs(args.out, exist_ok=True)
    manifest = []
    for i, p in enumerate(kept, start=args.start_index):
        gid = f"geom_{i:03d}"
        stl = os.path.join(args.out, gid, "constant", "triSurface", "disc.stl")
        os.makedirs(os.path.dirname(stl), exist_ok=True)
        poly = G.cross_section_polyline(p)
        G.revolve_to_stl(poly, stl, n_theta=args.n_theta)

        # per-geometry normalization + moment reference (uniform-density centroid)
        rr = np.linspace(0.0, p.R, 4000)
        h = np.clip(G.z_top(p, rr) - G.z_bot(p, rr), 0.0, None)
        zmid = 0.5 * (G.z_top(p, rr) + G.z_bot(p, rr))
        ringV = 2.0 * np.pi * rr * h
        cg_z = float(np.sum(ringV * zmid) / np.sum(ringV)) / 1000.0
        Rm = p.R / 1000.0
        cp = {"disc": gid, "diameter_m": round(2 * Rm, 6),
              "planform_area_m2": round(math.pi * Rm * Rm, 8),
              "cofr_m": [0.0, 0.0, round(cg_z, 6)],
              "Uinf": 26.75, "nu": 1.5e-5, "rhoInf": 1.225, "Re": 378000}
        json.dump(cp, open(os.path.join(args.out, gid, "case_params.json"), "w"), indent=2)

        rec = {"id": gid, "index": i,
               **{k: getattr(p, k) for k in
                  ("D_out_mm", "rim_width_mm", "rim_shoulder_mm", "parting_line_mm",
                   "dome_mm", "plate_thick_mm", "nose_radius_mm")}}
        manifest.append(rec)

    man_path = os.path.join(os.path.dirname(__file__), "geometries.json")
    if args.start_index > 0 and os.path.exists(man_path):
        prev = json.load(open(man_path))["geometries"]
        manifest = [g for g in prev if g["index"] < args.start_index] + manifest
    json.dump({"n": len(manifest), "ranges": RANGES, "geometries": manifest},
              open(man_path, "w"), indent=2)
    print(f"kept {len(kept)} legal geometries (from pool {len(pool)})")
    print(f"wrote STLs + case_params under {args.out}, manifest geometries.json")


if __name__ == "__main__":
    main()
