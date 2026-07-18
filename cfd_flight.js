// CFD-backed flight: geometry -> aero coefficients (CFD surrogate) -> 6-DOF
// trajectory. Faithful port of flight_sim.py (_deriv + RK4 simulate). Works in
// the browser (uses global CFD_SURROGATE) and in Node (module.exports) so the
// integrator can be verified against the Python sim on identical coefficients.
(function (root) {
  "use strict";
  var RHO = 1.225, G = 9.81, DEG = Math.PI / 180;
  var FADE_VREF = 14.0, FADE_GAIN = 0.060;
  // dynamic derivatives not produced by steady CFD — kept at Hummel values
  var CMq = -5.0e-3, CRr = 1.4e-2, CRp = -5.5e-3, CNr = -3.4e-5;

  var SUR = (typeof CFD_SURROGATE !== "undefined") ? CFD_SURROGATE
          : (typeof require !== "undefined" ? null : null);
  function surrogate() {
    return (typeof CFD_SURROGATE !== "undefined") ? CFD_SURROGATE : SUR;
  }

  // geometry params -> 7 CFD-surrogate aero coefficients (+ fixed damping terms)
  function cfdAero(p) {
    var S = surrogate();
    var x = [p.rim_width_mm, p.rim_shoulder_mm, p.parting_line_mm,
             p.dome_mm, p.plate_thick_mm, p.nose_radius_mm];
    var xs = x.map(function (v, i) { return (v - S.mean[i]) / S.std[i]; });
    // canonical 28-term degree-2 basis: [1] + x_i + x_i*x_j (i<=j)
    var t = [1];
    for (var i = 0; i < 6; i++) t.push(xs[i]);
    for (i = 0; i < 6; i++) for (var j = i; j < 6; j++) t.push(xs[i] * xs[j]);
    function dot(w) { var s = 0; for (var k = 0; k < w.length; k++) s += w[k] * t[k]; return s; }
    var T = S.targets;
    return {
      CL0: dot(T.CL0), CLa: dot(T.CLa),
      CD0: dot(T.CD0), CDa: dot(T.CDa), alpha0: dot(T.alpha0),
      CM0: dot(T.CM0), CMa: dot(T.CMa),
      CMq: CMq, CRr: CRr, CRp: CRp, CNr: CNr
    };
  }

  function rotB2I(phi, theta, psi) {
    var cph = Math.cos(phi), sph = Math.sin(phi);
    var cth = Math.cos(theta), sth = Math.sin(theta);
    var cps = Math.cos(psi), sps = Math.sin(psi);
    return [
      [cth * cps, sph * sth * cps - cph * sps, cph * sth * cps + sph * sps],
      [cth * sps, sph * sth * sps + cph * cps, cph * sth * sps - sph * cps],
      [-sth, sph * cth, cph * cth]
    ];
  }
  function matVec(M, v) {
    return [M[0][0]*v[0]+M[0][1]*v[1]+M[0][2]*v[2],
            M[1][0]*v[0]+M[1][1]*v[1]+M[1][2]*v[2],
            M[2][0]*v[0]+M[2][1]*v[1]+M[2][2]*v[2]];
  }
  function matTVec(M, v) {   // M^T * v
    return [M[0][0]*v[0]+M[1][0]*v[1]+M[2][0]*v[2],
            M[0][1]*v[0]+M[1][1]*v[1]+M[2][1]*v[2],
            M[0][2]*v[0]+M[1][2]*v[1]+M[2][2]*v[2]];
  }
  function cross(a, b) {
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
  }
  function dot3(a, b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
  function norm3(a) { return Math.sqrt(dot3(a, a)); }

  function deriv(s, phys, ac) {
    var vx = s[3], vy = s[4], vz = s[5];
    var phi = s[6], theta = s[7], psi = s[8], p = s[9], q = s[10], r = s[11];
    var V = [vx, vy, vz], vmag = norm3(V);
    if (vmag < 1e-3) return new Array(12).fill(0);
    var Vhat = [vx/vmag, vy/vmag, vz/vmag];

    var R = rotB2I(phi, theta, psi);
    var zb = matVec(R, [0, 0, 1]);
    var alpha = -Math.asin(Math.max(-1, Math.min(1, dot3(Vhat, zb))));

    var qdyn = 0.5 * RHO * vmag * vmag;
    var Lf = qdyn * phys.A * (ac.CL0 + ac.CLa * alpha);
    var Df = qdyn * phys.A * (ac.CD0 + ac.CDa * (alpha - ac.alpha0) * (alpha - ac.alpha0));

    var dragDir = [-Vhat[0], -Vhat[1], -Vhat[2]];
    var zbV = dot3(zb, Vhat);
    var perp = [zb[0]-zbV*Vhat[0], zb[1]-zbV*Vhat[1], zb[2]-zbV*Vhat[2]];
    var pn = norm3(perp) + 1e-12;
    var liftDir = [perp[0]/pn, perp[1]/pn, perp[2]/pn];
    var Fx = Lf*liftDir[0] + Df*dragDir[0];
    var Fy = Lf*liftDir[1] + Df*dragDir[1];
    var Fz = Lf*liftDir[2] + Df*dragDir[2] - phys.m*G;
    var ax = Fx/phys.m, ay = Fy/phys.m, az = Fz/phys.m;

    var k = qdyn * phys.A * phys.d;
    var fade = FADE_GAIN * Math.max(0, (FADE_VREF - vmag) / FADE_VREF);
    var CM0_eff = ac.CM0 - fade;
    var Mp = k * (CM0_eff + ac.CMa*alpha + ac.CMq * q * phys.d/(2*vmag));
    var Rl = k * (ac.CRr*r + ac.CRp*p) * phys.d/(2*vmag);
    var Nz = k * (ac.CNr*r) * phys.d/(2*vmag);

    var pAxis = cross(zb, Vhat); var pan = norm3(pAxis)+1e-12;
    pAxis = [pAxis[0]/pan, pAxis[1]/pan, pAxis[2]/pan];
    var rAxis = [Vhat[0]-dot3(Vhat,zb)*zb[0], Vhat[1]-dot3(Vhat,zb)*zb[1], Vhat[2]-dot3(Vhat,zb)*zb[2]];
    var ran = norm3(rAxis)+1e-12; rAxis = [rAxis[0]/ran, rAxis[1]/ran, rAxis[2]/ran];
    var tau = [Mp*pAxis[0]+Rl*rAxis[0]+Nz*zb[0],
               Mp*pAxis[1]+Rl*rAxis[1]+Nz*zb[1],
               Mp*pAxis[2]+Rl*rAxis[2]+Nz*zb[2]];
    var tb = matTVec(R, tau);

    var Id_ = phys.Id, Ia_ = phys.Ia;
    var pdot = (tb[0] - (Ia_-Id_)*q*r)/Id_;
    var qdot = (tb[1] - (Id_-Ia_)*r*p)/Id_;
    var rdot = tb[2]/Ia_;
    var phidot = p + (q*Math.sin(phi)+r*Math.cos(phi))*Math.tan(theta);
    var thetadot = q*Math.cos(phi) - r*Math.sin(phi);
    var psidot = (q*Math.sin(phi)+r*Math.cos(phi))/Math.cos(theta);

    return [vx, vy, vz, ax, ay, az, phidot, thetadot, psidot, pdot, qdot, rdot];
  }

  // p: geometry params; mp: mass props (Iz_gcm2, mass_g); opts: throw settings
  function simulateFlight(p, mp, opts) {
    opts = opts || {};
    var v0 = opts.v0 != null ? opts.v0 : 24.0;
    var launch = opts.launch_deg != null ? opts.launch_deg : 8.0;
    var hyzer = opts.hyzer_deg != null ? opts.hyzer_deg : 0.0;
    var nose = opts.nose_deg != null ? opts.nose_deg : 0.0;
    var h0 = opts.h0 != null ? opts.h0 : 1.4;
    var handed = opts.handed || "RHBH";
    var dt = opts.dt != null ? opts.dt : 1e-3;
    var tmax = opts.tmax != null ? opts.tmax : 15.0;

    var d = p.D_out_mm/1000, A = Math.PI*(d/2)*(d/2);
    var Ia = mp.Iz_gcm2 * 1e-7, Id = 0.5*Ia, m = mp.mass_g/1000;
    var phys = { m: m, d: d, A: A, Ia: Ia, Id: Id };
    var ac = cfdAero(p);

    var spinRev = v0/(Math.PI*d), spin = 2*Math.PI*spinRev;
    var r0 = (handed === "RHBH") ? -spin : spin;
    var la = launch*DEG;
    var theta0 = -(launch+nose)*DEG;
    var phi0 = (handed === "RHBH") ? -hyzer*DEG : hyzer*DEG;

    var s = [0,0,h0, v0*Math.cos(la),0,v0*Math.sin(la), phi0,theta0,0, 0,0,r0];
    var X=[0], Y=[0], Z=[h0], t=0;
    while (t < tmax && s[2] > 0) {
      var k1 = deriv(s, phys, ac);
      var s2 = s.map(function(v,i){return v+0.5*dt*k1[i];});
      var k2 = deriv(s2, phys, ac);
      var s3 = s.map(function(v,i){return v+0.5*dt*k2[i];});
      var k3 = deriv(s3, phys, ac);
      var s4 = s.map(function(v,i){return v+dt*k3[i];});
      var k4 = deriv(s4, phys, ac);
      s = s.map(function(v,i){return v+(dt/6)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);});
      t += dt;
      X.push(s[0]); Y.push(s[1]); Z.push(s[2]);
    }
    var yMax = 0; for (var i=0;i<Y.length;i++) if (Math.abs(Y[i])>Math.abs(yMax)) yMax=Y[i];
    var zMax = Math.max.apply(null, Z);
    return { x:X, y:Y, z:Z, coeffs:ac,
             distance_m:X[X.length-1], distance_ft:X[X.length-1]/0.3048,
             lateral_max_m:yMax, apex_m:zMax, flight_time_s:(X.length-1)*dt };
  }

  var api = { cfdAero: cfdAero, simulateFlight: simulateFlight };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  root.CFDFlight = api;
})(typeof window !== "undefined" ? window : this);
