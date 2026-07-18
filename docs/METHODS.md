# Methods: what is grounded vs. calibrated vs. heuristic

A "great disc" claim is only as trustworthy as the model behind it. This file
is the honest ledger of which parts of the pipeline are validated science,
which are calibrated to a plausible operating point, and which are placeholders
awaiting real data. Keep it current as the project evolves.

## Tier 1 — GROUNDED (published / physically exact)

| Component | Basis |
|---|---|
| Disc as a solid of revolution; volume, mass, moment of inertia, CG | Exact geometry integration (`disc_model.py`) |
| PDGA legality limits | PDGA Technical Standards (see `PDGA_technical_standards.md`) |
| 6-DOF rigid-body equations of motion, incl. gyroscopic precession | Standard rigid-body dynamics; matches Hummel (2003) formulation |
| Aerodynamic force/moment **functional forms** | Hummel (2003), UC Davis thesis |
| Pitch/roll/spin **moment coefficients** (CMα, CMq, CRr, CRp, CNr) | Hummel (2003) fitted values |
| Disc mass & moments of inertia fed to the sim | Our own geometry model (exact) |
| Pro-disc preset targets: diameter, height, rim depth, rim width, inside rim diameter, max weight | PDGA certified-disc database export (`data/pdga_approved_discs.csv`, 2,399 discs, fetched July 2026) |
| "Closest certified discs" matcher in the app | Nearest-neighbor over the same PDGA measurements (`disc_db.js`, regenerate with `scripts/build_disc_db.py`) |

The flight simulator's *behavior* is validated qualitatively: it produces
lift-supported glide, realistic ballistic altitude arcs, correct hyzer/anhyzer
response, gyroscopic precession, and stability-driven drift — none of it
hand-scripted; it falls out of the equations.

## Tier 2 — CALIBRATED (plausible, first-order, needs CFD to confirm)

| Component | Status |
|---|---|
| Lift/drag coefficient **magnitudes** (CL0, CD0) for the disc-golf regime | Calibrated so a driver flies a realistic distance. Hummel's raw values are for a larger, blunter Ultimate disc and would under-fly a disc-golf driver. |
| Geometry → coefficient scaling (rim width → drag, dome → lift, parting-line/dome/rim → stability) | **Superseded July 2026 by the CFD surrogate** (`cfd_surrogate.js`; see the roadmap below) — geometry→coefficients is now CFD-derived rather than hand-tuned slopes. The old `aero_from_geometry()` estimates remain only as the Python fallback. |
| Low-speed fade term (`LOWSPEED_FADE_*`) | Stands in for advance-ratio-dependent moments that a constant-coefficient model can't capture. Real, but its magnitude is tuned, not measured. |
| Throw speeds used in the preset demo (14–27 m/s) | Intermediate-amateur arm speeds → intermediate distances. Pro speeds (30+ m/s) + more spin go much farther. |
| Pro-disc preset **derivation** of unpublished dimensions | PDGA publishes diameter/height/rim depth/rim width but not plate thickness, parting-line height, or nose radius. The app estimates those per class (`CLASS_EST`), estimates parting line from the disc's printed stability, and fits plastic density so model mass lands ~0.5 g under the certified max weight. Certified numbers are reproduced exactly; the estimated ones are labeled as such in the UI. |
| Manufacturer flight numbers shown on presets | Printed ratings (marketing scale), not measurements — shown because athletes navigate by them. |
| The app's live flight estimate (`flightEstimate()` / `flight_estimate()`) | Linear models fit to the manufacturer flight numbers of the 43 pro-disc presets (`scripts/fit_flight_numbers.py` → `data/flight_fit.json`, `flight_fit.js`). Leave-one-out RMSE: speed 0.81, glide 0.74, turn 0.46, fade 0.45. Caveats: n=43; targets are a marketing scale; the parting-line feature is **circular** (preset parting lines were derived from printed stability, so the fit encodes our assumption rather than confirming it); feature sets chosen partly by LOO on the same data. Trust it near real-disc geometry; it extrapolates poorly. |

## Tier 3 — HEURISTIC / PLACEHOLDER (replace ASAP)

| Component | Status |
|---|---|
| "Derived flight numbers" from simulated trajectories | Reasonable descriptors, but they inherit the Tier-2 calibration uncertainty. |
| Flight estimate outside the fitted envelope (rim width ≳ certified range, exotic domes/parting lines) | The Tier-2 linear fit is only anchored by 43 real discs; beyond them it's extrapolation — treat as placeholder. |

## The path to remove the caveats

1. **CFD sweep.** ✅ DONE (July 2026). Steady RANS (k-ω SST) run on the
   parametric profile across angle of attack (−5°…+15°) for 96 PDGA-legal
   geometries on a Slurm cluster (OpenFOAM, `cfd/`). Validated against the
   Potts / Kamaruddin wind-tunnel data on the Aviar putter: **trim angle (7.8°
   vs 7.5°), aerodynamic-centre position (xac/d 0.055 vs 0.05), and pitching-
   moment gradient all match.** Caveat: absolute CL/CD magnitudes run ~15–20%
   high — the coarse wall-function mesh (y+≈28) trades magnitude accuracy for
   throughput. (A finer mesh dropped y+ into the invalid buffer layer and made
   things worse; the coarse mesh is the correct wall-function regime.)
2. **Surrogate.** ✅ DONE. A standardized degree-2 polynomial maps the 6 shape
   parameters → the 7 static aero coefficients (`cfd/bigsweep/`, exported to
   `cfd_surrogate.js`, ~2.6 KB). 8-fold CV R²: CMa 0.90, CM0 0.89, CLa 0.84,
   CL0 0.80, CD0 0.72, α0 0.47, CDa 0.36 — moment/lift terms strong, drag-bucket
   terms noisy (coarse-mesh force scatter). The designer app now computes a
   **real CFD-surrogate 6-DOF flight path** from it: `cfd_flight.js` is a JS
   port of `flight_sim.py` verified to match the Python integrator to 5 decimals
   on identical coefficients. The path is drawn in blue ("CFD physics") beside
   the manufacturer rating.
3. **Re-validate & anchor** — STILL TODO: (a) bias-correct the ~15–20% CL/CD
   level against the three discs Potts measured (Aviar / Roc / Wraith);
   (b) validate the 6-DOF *trajectories* against published disc-golf paths (the
   2022 Springer CFD + rigid-body paper) and field data; (c) optionally expand
   the sweep to ~300 geometries to lift the weak drag-bucket coefficients.

**Status:** geometry→coefficients is now CFD-derived, not hand-calibrated — a
real step from Tier 2 toward Tier 1. It is **not yet fully Tier-1**: the coarse-
mesh magnitude bias is uncorrected and the trajectories aren't yet validated
against real flights. The app's flight NUMBERS panel (speed/glide/turn/fade)
still uses the Tier-2 ridge fit — mapping CFD coefficients → those four numbers
is future work; for now the CFD result is expressed as the drawn trajectory,
predicted distance, and apex.
