"""Anchor-disc CFD cases: the three golf discs Kamaruddin measured in the
Manchester wind tunnel (Aviar putter, Roc mid-range, Wraith driver).

Geometry comes from Kamaruddin (2011) Table 4.1 — the actual tested discs —
mapped through the same params_from_certified() the designer app uses, so the
anchor CFD disc == the app's model of that disc. Note the tested Roc is an
older mold (d=212 mm, rim 9 mm), not the modern PDGA entry (217/12).

Writes bigsweep-style case dirs:
    geoms/anchor_<name>/constant/triSurface/disc.stl
    geoms/anchor_<name>/case_params.json
and a manifest anchors.json with the DiscParams used.
"""
import json, math, os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import geometry as G

# certified-style dicts from Kamaruddin Table 4.1 (mm); flight numbers are the
# discs' well-known ratings (only stab = turn+fade feeds the parting-line map)
ANCHORS = {
    "anchor_aviar":  {"cls": "Putter",          "name": "Aviar",  "dia": 212, "h": 20, "rd": 15, "rw": 9,  "flight": [2, 3, 0, 1]},
    "anchor_roc":    {"cls": "Midrange",        "name": "Roc",    "dia": 212, "h": 20, "rd": 12, "rw": 9,  "flight": [4, 4, 0, 3]},
    "anchor_wraith": {"cls": "Distance Driver", "name": "Wraith", "dia": 211, "h": 14, "rd": 12, "rw": 21, "flight": [11, 5, -1, 3]},
}


def main():
    base = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(base, "geoms")
    manifest = []
    for gid, d in ANCHORS.items():
        p = G.params_from_certified(d)
        stl = os.path.join(out, gid, "constant", "triSurface", "disc.stl")
        os.makedirs(os.path.dirname(stl), exist_ok=True)
        poly = G.cross_section_polyline(p)
        G.revolve_to_stl(poly, stl, n_theta=480)

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
        json.dump(cp, open(os.path.join(out, gid, "case_params.json"), "w"), indent=2)

        manifest.append({"id": gid, "disc": d["name"],
                         **{k: getattr(p, k) for k in
                            ("D_out_mm", "rim_width_mm", "rim_shoulder_mm", "parting_line_mm",
                             "dome_mm", "plate_thick_mm", "nose_radius_mm")}})
        print(gid, manifest[-1])

    json.dump({"anchors": manifest}, open(os.path.join(base, "anchors.json"), "w"), indent=2)
    print(f"wrote {len(manifest)} anchor cases under {out}")


if __name__ == "__main__":
    main()
