# Sources & Annotated Bibliography

Research basis for the parametric disc project. Grouped by role. Links verified
July 2026.

## Directly on-point: parametric disc golf design + CFD

- **Optimal design for disc golf by computational fluid dynamics and machine
  learning** — *Structural and Multidisciplinary Optimization* (Springer, 2021).
  https://link.springer.com/article/10.1007/s00158-021-03107-7
  Essentially our idea, proven: a disc profile parameterized by ~19 control
  points fitted with cubic B-splines and revolved into 3D (12 design variables);
  hundreds of shapes run through RANS CFD (k-ω SST); regression-tree-ensemble
  surrogate models trained on ~700 designs (~2,800 CFD cases) for fast search.
  Key finding: min-drag → wide-rim "wedge/football" drivers; max-lift → rimless
  domes; the two goals are contradictory (Pareto trade-off). **Read this first.**

- **Disc golf trajectory modelling combining computational fluid dynamics and
  rigid body dynamics** — *Sports Engineering* (Springer, 2022).
  https://link.springer.com/article/10.1007/s12283-022-00390-5
  Couples CFD-derived aero coefficients to 6-DOF rigid-body dynamics to predict
  actual trajectories. The blueprint for our "flight simulator" tier.

## Foundational flying-disc aerodynamics (coefficient data)

- **Potts, J. R. — Disc-wing Aerodynamics / Aerodynamic Performance of Flying
  Discs** (Sheffield Hallam University).
  https://shura.shu.ac.uk/14521/1/Potts%20-%20Aerodynamic%20Performance%20of%20Flying%20Discs%20-%20(AM).pdf
  Landing page: https://shura.shu.ac.uk/14521/
  The most thorough wind-tunnel study of flying-disc aerodynamics: lift, drag,
  and pitching-moment coefficients vs angle of attack, spin effects, cavity
  effects. Primary source for measured aero coefficients.

- **Potts & Crowther — Frisbee™ Aerodynamics / Aerodynamic performance of flying
  discs**, *Aircraft Engineering and Aerospace Technology* (Emerald).
  https://www.emerald.com/aeat/article/90/2/390/48679/Aerodynamic-performance-of-flying-discs

- **Hummel, S. A. — Frisbee Flight Simulation and Throw Biomechanics** (M.S.
  thesis, UC Davis, 2003).
  https://morleyfielddgc.wordpress.com/wp-content/uploads/2009/04/hummelthesis.pdf
  Lab page: https://research.engineering.ucdavis.edu/biosport/sample-page/test-page-1/frisbee-flight-simulation-and-throw-biomechanics/
  The canonical 6-DOF flying-disc flight simulation with measured aero
  coefficients. The reference implementation for our trajectory model.

## Data sets

- **PDGA certified-disc database** — export of every PDGA-approved disc with
  certified measurements (max weight, diameter, height, rim depth, inside rim
  diameter, rim thickness, rim config, flexibility).
  https://www.pdga.com/technical-standards/equipment-certification/discs/export
  Snapshot: `data/pdga_approved_discs.csv` (2,401 rows, fetched July 2026).
  Note: the endpoint 403s on curl's default User-Agent — pass a browser UA.
  `scripts/build_disc_db.py` slims it into `disc_db.js` for the designer app,
  which uses it for pro-disc presets and the closest-certified-disc matcher.

## Supporting / explanatory

- **The Aerodynamics and Stability of Flying Discs** — Stanford PH210 course note.
  http://large.stanford.edu/courses/2007/ph210/scodary1/
- **Gone With the Wind: An Investigation into the Flight Dynamics of Discs** —
  Pozderac, junior thesis (Wooster).
  https://physics.wooster.edu/wp-content/uploads/2021/08/Junior-IS-Thesis-Web_2016_Pozderac.pdf
- **Disc Golf and its Relation to Aerodynamics and Mechanics** — Harp (Arkansas).
  https://bpb-us-e1.wpmucdn.com/wordpressua.uark.edu/dist/a/146/files/2025/04/DiscGolf.pdf
- **An Aerodynamic Simulation of Disc Flight** — CSB/SJU honors thesis.
  https://digitalcommons.csbsju.edu/cgi/viewcontent.cgi?article=1067&context=honors_theses
- **Disc Golf Guide: The Science of Flight** — Disc Golf United blog (advance
  ratio, rim-to-plate mass ratio, center-of-pressure vs center-of-mass).
  https://blog.discgolfunited.com/the-science-of-flight/

## Regulatory

- **PDGA Technical Standards — Guidelines** (the disc legality rules).
  https://www.pdga.com/technical-standards/guidelines
- **PDGA Equipment Certification — Discs** (overview of approval process).
  https://www.pdga.com/technical-standards/equipment-certification/discs/guide
- **Why & How Discs Get PDGA Approved** — UDisc blog (readable overview).
  https://udisc.com/blog/post/why-how-discs-get-pdga-approved
