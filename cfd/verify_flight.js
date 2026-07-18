// Node side of the flight-port verification: run cfd_flight.js on the same
// test case (written by verify_flight.py) and print endpoints to compare.
const fs = require("fs");
const path = require("path");
const ROOT = path.dirname(__dirname);

// load the surrogate as a global (cfd_flight.js reads global CFD_SURROGATE)
eval(fs.readFileSync(path.join(ROOT, "cfd_surrogate.js"), "utf8")
       .replace("const CFD_SURROGATE", "global.CFD_SURROGATE"));
const { cfdAero, simulateFlight } = require(path.join(ROOT, "cfd_flight.js"));

const c = JSON.parse(fs.readFileSync(path.join(ROOT, "cfd", "test_case.json"), "utf8"));
const res = simulateFlight(c.params, c.mp, c.opts);
const ac = cfdAero(c.params);

const r = (x) => Math.round(x * 100000) / 100000;
console.log("JS (cfd_flight.js):");
console.log("  coeffs:", Object.fromEntries(
  ["CL0", "CLa", "CD0", "CDa", "alpha0", "CM0", "CMa"].map((k) => [k, r(ac[k])])));
console.log(`  distance_ft = ${res.distance_ft.toFixed(2)}`);
console.log(`  lateral_m   = ${res.lateral_max_m >= 0 ? "+" : ""}${res.lateral_max_m.toFixed(3)}`);
console.log(`  apex_m      = ${res.apex_m.toFixed(3)}`);
console.log(`  time_s      = ${res.flight_time_s.toFixed(3)}`);
