# Methods: what is grounded vs. calibrated vs. heuristic

A "great disc" claim is only as trustworthy as the model behind it. This file
is the honest ledger of which parts of the pipeline are validated science,
which are calibrated to a plausible operating point, and which are placeholders
awaiting real data. Keep it current as the project evolves.

## Tier 1: GROUNDED (published / physically exact)

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
response, gyroscopic precession, and stability-driven drift. None of it is
hand-scripted; it falls out of the equations.

## Tier 2: CALIBRATED (plausible, first-order, needs CFD to confirm)

| Component | Status |
|---|---|
| Lift/drag coefficient **magnitudes** (CL0, CD0) for the disc-golf regime | **Superseded July 2026:** the CFD-surrogate path is now anchored to Kamaruddin's wind-tunnel measurements (kCL = 0.85, kCD = 0.98; see roadmap 3a). The hand-calibrated values remain only in the Python fallback `aero_from_geometry()`. |
| Geometry → coefficient scaling (rim width → drag, dome → lift, parting-line/dome/rim → stability) | **Superseded July 2026 by the CFD surrogate** (`cfd_surrogate.js`; see the roadmap below): geometry→coefficients is now CFD-derived rather than hand-tuned slopes. The old `aero_from_geometry()` estimates remain only as the Python fallback. |
| Low-speed fade term (`LOWSPEED_FADE_*`) | Stands in for advance-ratio-dependent moments that a constant-coefficient model can't capture. Real, but its magnitude is tuned, not measured. |
| Throw speeds used in the preset demo (14–27 m/s) | Intermediate-amateur arm speeds → intermediate distances. Pro speeds (30+ m/s) + more spin go much farther. |
| Pro-disc preset **derivation** of unpublished dimensions | PDGA publishes diameter/height/rim depth/rim width but not plate thickness, parting-line height, or nose radius. The app estimates those per class (`CLASS_EST`), estimates parting line from the disc's printed stability, and fits plastic density so model mass lands ~0.5 g under the certified max weight. Certified numbers are reproduced exactly; the estimated ones are labeled as such in the UI. |
| Manufacturer flight numbers shown on presets | Printed ratings (marketing scale), not measurements, shown because athletes navigate by them. |
| The app's live flight estimate (`flightEstimate()` / `flight_estimate()`) | **Now CFD-driven (July 2026, roadmap 4):** the numbers come from the surrogate's aero coefficients (`scripts/fit_cfd_flight_numbers.py` → `cfd_flight_fit.js`). The older geometry ridge fit (`scripts/fit_flight_numbers.py`, LOO RMSE: speed 0.81, glide 0.74, turn 0.46, fade 0.45) remains as fallback. Shared caveats: n=43; targets are a marketing scale; feature sets chosen partly by LOO on the same data; the preset parting lines were themselves derived from printed stability (circularity, which now enters through the geometry the surrogate is evaluated at rather than as a direct feature). Trust it near real-disc geometry; it extrapolates poorly. |

## Tier 3: HEURISTIC / PLACEHOLDER (replace ASAP)

| Component | Status |
|---|---|
| "Derived flight numbers" from simulated trajectories | Reasonable descriptors, but they inherit the Tier-2 calibration uncertainty. |
| Flight estimate outside the fitted envelope (rim width ≳ certified range, exotic domes/parting lines) | The Tier-2 linear fit is only anchored by 43 real discs; beyond them it's extrapolation; treat it as placeholder. |

## The path to remove the caveats

1. **CFD sweep.** ✅ DONE (July 2026). Steady RANS (k-ω SST) run on the
   parametric profile across angle of attack (−5°…+15°) for 96 PDGA-legal
   geometries on a Slurm cluster (OpenFOAM, `cfd/`). Validated against the
   Potts / Kamaruddin wind-tunnel data on the Aviar putter: **trim angle (7.8°
   vs 7.5°), aerodynamic-centre position (xac/d 0.055 vs 0.05), and pitching-
   moment gradient all match.** Caveat: absolute CL/CD magnitudes run ~15–20%
   high: the coarse wall-function mesh (y+≈28) trades magnitude accuracy for
   throughput. (A finer mesh dropped y+ into the invalid buffer layer and made
   things worse; the coarse mesh is the correct wall-function regime.)
2. **Surrogate.** ✅ DONE. A standardized degree-2 polynomial maps the 6 shape
   parameters → the 7 static aero coefficients (`cfd/bigsweep/`, exported to
   `cfd_surrogate.js`, ~2.6 KB). 8-fold CV R²: CMa 0.90, CM0 0.89, CLa 0.84,
   CL0 0.80, CD0 0.72, α0 0.47, CDa 0.36; moment/lift terms strong, drag-bucket
   terms noisy (coarse-mesh force scatter). The designer app now computes a
   **real CFD-surrogate 6-DOF flight path** from it: `cfd_flight.js` is a JS
   port of `flight_sim.py` verified to match the Python integrator to 5 decimals
   on identical coefficients. The path is drawn in blue ("CFD physics") beside
   the manufacturer rating.
3. **Re-validate & anchor** (July 2026):
   (b) ✅ **Trajectory validation.** Kamaruddin's PhD thesis (Manchester 2011,
   the full data behind the 2018 paper) publishes 6-DOF simulated flights of
   the Aviar/Roc/Wraith driven by her measured wind-tunnel coefficients:
   ranges 49 / 54 / 63 m at a fixed launch (20 m/s, pitch 15°, AdvR 0.5,
   ground level). Driving OUR integrator with HER measured coefficients
   (`cfd/verify_trajectory.py`, coefficients in
   `data/kamaruddin_wind_tunnel.json`) isolates the trajectory model from our
   CFD. Result: with model assumptions matched to hers (no aero damping
   moments, constant spin; the terms her sim omits), our ranges are **49.6 / 52.6
   / 54.2 m vs her 49 / 54 / 63 m** (putter +1%, mid −3%, driver −14%), with
   matching apex heights, lateral scale, and the driver's S-shaped lateral
   reversal. The driver gap traces to her yaw-rate ≃ spin-rate simplification
   (its effect grows with flight length) plus the figure-read CD curvature
   (±0.02 → ±2 m). With our full physics (Hummel damping + spin decay) all
   three fly 12–15% shorter than her idealized sim, the same direction the
   2022 Sports Engineering study reports for its own model vs real throws.
   Verdict: the 6-DOF integrator is sound; results live in
   `cfd/trajectory_validation.json`.
   (a) ✅ **CL/CD bias anchor.** Kamaruddin's thesis also publishes the
   absolute wind-tunnel coefficients for the discs she tested (Table 5.6) and
   their measured geometries (Table 4.1). We ran her three representative
   discs (Aviar / Roc / Wraith, built from her measured dimensions through
   the app's own `params_from_certified`) through the identical big-sweep CFD
   pipeline (`cfd/anchors/`, 21 solves on the validated coarse mesh) and
   compared levels: **kCL = 0.850, kCD = 0.980** (pooled least-squares,
   measured = k × CFD). So the coarse mesh runs lift ~15% high and drag only
   ~2% high. Per-disc agreement is strongest exactly where it matters for
   behaviour: Wraith CD@0 0.054 vs 0.055 measured, CMa 0.0061 vs 0.006/deg,
   trim 4.0° vs 4°; trim matches on all three discs. The exported surrogate
   (`cfd_surrogate.js`) now ships with these factors applied to the CL and CD
   targets (`bias_correction` field; α₀ is scale-invariant), so every
   downstream consumer (designer flight path, wind-tunnel charts, flight
   numbers) sits on the wind-tunnel level. Full comparison in
   `cfd/anchors/anchor_bias.json`.
   (c) ⏳ **Sweep expansion 96 → 300 geometries.** Generation + meshing +
   1,428 additional solves running on the cluster (`cfd/bigsweep/`, ids
   geom_096–299, fresh LHS seed). When it finishes:
   `bash cfd/bigsweep/refresh_after_expansion.sh` on the cluster, then the
   local re-export steps in that script's header. Expected to mainly lift the
   weak drag-bucket terms (CDa CV R² 0.36, α₀ 0.47).
4. **CFD flight numbers.** ✅ DONE (July 2026). The app's speed/glide/turn/
   fade panel no longer uses the Tier-2 geometry ridge fit: `flightEstimate()`
   now computes the four numbers from the **surrogate's aerodynamic
   coefficients** (drag floor, lift, pitching moment at cruise and late-flight
   α, trim angle) via small ridge models calibrated on the same 43 pro discs
   (`scripts/fit_cfd_flight_numbers.py` → `cfd_flight_fit.js`,
   `data/cfd_flight_fit.json`). LOO RMSE: speed 0.95, glide 0.86, turn 0.66,
   fade 0.67 (constant-predictor baselines: 3.7 / 1.3 / 0.9 / 1.0). Verified
   JS ≡ Python on the anchor discs (Aviar rated 2/3/0/1 → 1/2.9/−0.3/0.8;
   Wraith rated 11/5/−1/3 → 12.5/5.3/−0.9/1.9). The Tier-2 fit remains as
   fallback when the CFD files are absent. Same caveats as before: ratings
   are a marketing scale; feature sets LOO-selected on the same data.

**Status (baseline complete, July 2026):** the full chain geometry → CFD
coefficients → wind-tunnel-anchored magnitudes → 6-DOF trajectory → flight
numbers is now CFD-derived end to end. The magnitude bias is corrected against
Kamaruddin's measurements, the integrator is validated against her published
flights, and the flight-numbers panel runs on surrogate aerodynamics. Remaining
known gaps: the trajectory model is validated against a published *simulation*
driven by measured coefficients, not against instrumented field throws; the
drag-bucket surrogate terms stay noisy until the 300-geometry expansion lands;
spin-dependent (advance-ratio) effects are still parameterized, not computed;
spinning CFD (MRF) is the next tier.
