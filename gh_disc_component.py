# =====================================================================
#  Parametric Disc Golf Disc  --  Grasshopper Python component (Rhino 8)
# =====================================================================
#  HOW TO USE
#  1. In Grasshopper (Rhino 8) drop a "Python 3 Script" component.
#  2. Add Number inputs named EXACTLY (right-click component > zoom in > +):
#        D_out, rim_width, rim_shoulder, parting_line,
#        dome, plate_thick, nose_radius, density
#     Wire a Number Slider to each (ranges suggested below).
#  3. Add outputs named:  disc, mass_g, Iz_gcm2, legal, report
#  4. Paste this whole file into the component's code editor.
#
#  Suggested slider ranges (mm, except density g/cc):
#     D_out        210 .. 214     (putter/mid ~212)
#     rim_width      9 .. 26       (putter ~11, driver ~24)
#     rim_shoulder  10 .. 22
#     parting_line   5 .. 16
#     dome           0 .. 8
#     plate_thick  1.2 .. 4.5
#     nose_radius  1.6 .. 4.0
#     density     0.90 .. 1.00
# =====================================================================
import math
import Rhino
import Rhino.Geometry as rg

# ---- defaults so the component runs even before sliders are wired ----
def _d(name, val):
    return globals()[name] if name in globals() and globals()[name] is not None else val

D_out        = float(_d("D_out", 212.0))
rim_width    = float(_d("rim_width", 11.5))
rim_shoulder = float(_d("rim_shoulder", 16.0))
parting_line = float(_d("parting_line", 9.0))
dome         = float(_d("dome", 5.5))
plate_thick  = float(_d("plate_thick", 2.4))
nose_radius  = float(_d("nose_radius", 2.0))
density      = float(_d("density", 0.98))

R    = D_out / 2.0
R_in = R - rim_width
TOP_RIM_POWER = 2.0
DOME_POWER    = 2.0

# ---------- profile surfaces (must match disc_model.py) ----------
def z_top(r):
    if r <= R_in:
        x = min(max(r / max(R_in, 1e-9), 0.0), 1.0)
        return rim_shoulder + dome * (1.0 - x ** DOME_POWER)
    u = min(max((r - R_in) / max(rim_width, 1e-9), 0.0), 1.0)
    drop = rim_shoulder - parting_line
    return parting_line + drop * (1.0 - u) ** (1.0 / TOP_RIM_POWER)

def z_bot(r):
    if r <= R_in:
        return z_top(r) - plate_thick
    if r >= (R - nose_radius):
        t = min(max((r - (R - nose_radius)) / max(nose_radius, 1e-9), 0.0), 1.0)
        return parting_line * (1.0 - math.sqrt(max(1.0 - t * t, 0.0)))
    return 0.0

# ---------- build the closed cross-section in the XZ plane ----------
N = 240
rs = [R * i / (N - 1) for i in range(N)]
pts = []
# bottom, centre -> edge
for r in rs:
    pts.append(rg.Point3d(r, 0.0, z_bot(r)))
# top, edge -> centre
for r in reversed(rs):
    pts.append(rg.Point3d(r, 0.0, z_top(r)))
pts.append(pts[0])  # close

profile = rg.PolylineCurve(pts)

# ---------- revolve 360 deg about the world Z axis ----------
axis = rg.Line(rg.Point3d(0, 0, 0), rg.Point3d(0, 0, 1))
rev  = rg.RevSurface.Create(profile, axis, 0.0, 2.0 * math.pi)
disc = rg.Brep.CreateFromRevSurface(rev, True, True) if rev else None

# ---------- mass properties from the actual solid ----------
mass_g = None
Iz_gcm2 = None
cavity_mm = z_top(R_in) - plate_thick     # underside inner-rim depth
if disc and disc.IsValid:
    vmp = rg.VolumeMassProperties.Compute(disc)
    if vmp:
        vol_mm3 = vmp.Volume
        mass_g = vol_mm3 / 1000.0 * density              # cm^3 * g/cc
        # world moment of inertia about Z through origin (kg*mm^2-ish units)
        # convert: density folds in via mass; VMP returns geometric moments.
        moi = vmp.WorldCoordinatesMomentsOfInertia       # about world axes
        # geometric Iz (mm^5) -> multiply by density(g/mm^3) -> g*mm^2, then ->g*cm^2
        rho_g_mm3 = density / 1000.0
        Iz_gcm2 = moi.Z * rho_g_mm3 / 100.0

# ---------- PDGA legality ----------
D_cm = D_out / 10.0
cav_cm = cavity_mm / 10.0
inner_dia_cm = 2.0 * R_in / 10.0
checks = []
def chk(name, ok, val, lim):
    checks.append((name, bool(ok), val, lim))

chk("diameter", 21.0 <= D_cm <= 30.0, f"{D_cm:.1f}cm", "21-30")
chk("rim_depth", (0.05*D_cm <= cav_cm <= 0.12*D_cm) and cav_cm >= 1.10,
    f"{cav_cm:.2f}cm", f"{0.05*D_cm:.2f}-{0.12*D_cm:.2f} & >=1.1")
chk("rim_width", rim_width/10.0 <= 2.60, f"{rim_width/10.0:.2f}cm", "<=2.6")
chk("inner_rim_dia", inner_dia_cm >= 15.80, f"{inner_dia_cm:.1f}cm", ">=15.8")
chk("plate_thick", plate_thick/10.0 <= 0.50, f"{plate_thick/10.0:.2f}cm", "<=0.5")
if mass_g is not None:
    chk("weight_per_cm", mass_g <= 8.30*D_cm, f"{mass_g/D_cm:.2f}", "<=8.3")
    chk("weight_abs", mass_g <= 200.0, f"{mass_g:.0f}g", "<=200")
chk("edge_radius", nose_radius >= 1.60, f"{nose_radius:.1f}mm", ">=1.6")

legal = all(c[1] for c in checks)
lines = [f"{'PASS' if ok else 'FAIL'}  {n:<14} {val:<10} (limit {lim})"
         for n, ok, val, lim in checks]
if mass_g is not None:
    lines.append(f"----  mass {mass_g:.1f} g   Iz {Iz_gcm2:.0f} g*cm^2")
lines.append(f"====  {'PDGA LEGAL' if legal else 'OUT OF SPEC'}")
report = "\n".join(lines)
