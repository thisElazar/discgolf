"""Stamp a copied case for one angle of attack.

Disc fixed in the x-y plane (normal +z); freestream rotated by alpha.
Nose-up disc in a +x flow  <=>  freestream = Uinf(cos a, 0, sin a):
    U        = Uinf (cos a, 0, sin a)
    liftDir  = (-sin a, 0, cos a)      # perpendicular to U, +z at a=0
    dragDir  = ( cos a, 0, sin a)      # along U
Rewrites case/0/U and case/system/forceCoeffs.
"""
import argparse, json, math, os


def w(path, text):
    with open(path, "w") as f:
        f.write(text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--alpha", type=float, required=True, help="angle of attack, degrees")
    ap.add_argument("--params", required=True)
    args = ap.parse_args()

    P = json.load(open(args.params))
    U = P["Uinf"]
    a = math.radians(args.alpha)
    ux, uz = U * math.cos(a), U * math.sin(a)
    lx, lz = -math.sin(a), math.cos(a)
    dx, dz = math.cos(a), math.sin(a)
    cx, cy, cz = P["cofr_m"]
    A = P["planform_area_m2"]
    d = P["diameter_m"]

    w(os.path.join(args.case, "0", "U"), f"""FoamFile {{ version 2.0; format ascii; class volVectorField; object U; }}
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform ({ux:.6f} 0 {uz:.6f});
boundaryField
{{
    farfield
    {{
        type            freestreamVelocity;
        freestreamValue uniform ({ux:.6f} 0 {uz:.6f});
    }}
    disc
    {{
        type            noSlip;
    }}
}}
""")

    w(os.path.join(args.case, "system", "forceCoeffs"), f"""forceCoeffs1
{{
    type            forceCoeffs;
    libs            (forces);
    writeControl    timeStep;
    writeInterval   1;
    log             yes;
    patches         (disc);
    rho             rhoInf;
    rhoInf          {P['rhoInf']};
    magUInf         {U};
    lRef            {d};
    Aref            {A:.8f};
    liftDir         ({lx:.6f} 0 {lz:.6f});
    dragDir         ({dx:.6f} 0 {dz:.6f});
    CofR            ({cx} {cy} {cz});
    pitchAxis       (0 1 0);
}}
""")
    print(f"stamped alpha={args.alpha:+.1f}  U=({ux:.3f},0,{uz:.3f})  "
          f"liftDir=({lx:.3f},0,{lz:.3f})  dragDir=({dx:.3f},0,{dz:.3f})")


if __name__ == "__main__":
    main()
