# PDGA Disc Technical Standards (design reference)

These are the measurable requirements a disc must satisfy to be **PDGA-approved**
for sanctioned play. They define the "box" our parametric model must stay inside.
Every numeric limit below is enforced live in `disc_model.py`, the Grasshopper
component, and the slider app.

Source: PDGA Technical Standards — Guidelines
(https://www.pdga.com/technical-standards/guidelines). Verify against the current
published standard before submitting a disc for certification; the PDGA updates
these periodically.

## Dimensional limits

| Property | Requirement | Notes |
|---|---|---|
| Outside diameter | **21.0 – 30.0 cm** | Real disc golf discs cluster at 21.0–21.4 cm |
| Underside inner rim depth | **5–12% of outside diameter, AND ≥ 1.1 cm** | The concave underside cavity depth |
| Rim width | **≤ 2.6 cm** | Radial width of the rim; drivers push this limit |
| Inner rim diameter | **≥ 15.8 cm** | Circular inner rim minimum |
| Flight-plate thickness | **≤ 0.5 cm** (standard areas) | See thickened-center exception below |
| Thickened center section | **≤ 1.2 cm**, over a **5–10 cm** diameter, centered, with **≤ 60° transition slope** | Allows a raised/bossed center |
| Rim-area surface elevation | **≤ 3 mm above the outermost edge** | No bumps proud of the outer edge |

## Weight limits

| Property | Requirement |
|---|---|
| Weight per cm of diameter | **≤ 8.3 g per cm** of outside diameter |
| Absolute maximum weight | **≤ 200 g** |

For a 21.2 cm disc the per-cm rule caps weight at **~176 g** — which is exactly
why "max weight" discs sit around 175 g.

## Physical / material tests

| Test | Requirement | What it checks |
|---|---|---|
| Leading-edge radius | Must pass a **1/16" (1.6 mm) radius gauge** | Edge not too sharp/pointed |
| Flexibility | Rating **≤ 27 lb (12.25 kg)** | Disc not too rigid |
| Rim configuration | Rating **≥ 26.0** | A PDGA rim-geometry gauge measurement |

## Modeling notes / caveats

- **Rim configuration rating (≥ 26.0):** the PDGA measures this with a specific
  gauge and does not publish the exact formula, so our tools report the raw rim
  geometry (depth, width) as *informational* rather than asserting a pass/fail we
  can't compute reliably. Treat it as a physical-sample check at certification.
- **Flexibility (≤ 27 lb):** a material + wall-stiffness property, not a pure
  geometry output. It depends on resin choice and plate thickness, so it becomes
  a real constraint once we pick materials.
- The 8 checks our tools *do* enforce live: outside diameter, rim depth (% and
  absolute), rim width, inner rim diameter, plate thickness, weight-per-cm,
  absolute weight, and leading-edge radius.
