// Shared disc geometry model: solid-of-revolution profile, mass properties,
// archetype presets, and the certified-measurement conversion. Mirrors
// disc_model.py. Used by the designer, wind tunnel, and flight sim pages.
const clamp=(v,lo,hi)=>Math.min(Math.max(v,lo),hi);
function derived(p){const R=p.D_out_mm/2;return{R,R_in:R-p.rim_width_mm};}
function zTop(p,r){const{R,R_in}=derived(p);if(r<=R_in){const x=clamp(r/Math.max(R_in,1e-9),0,1);return p.rim_shoulder_mm+p.dome_mm*(1-Math.pow(x,2));}const u=clamp((r-R_in)/Math.max(p.rim_width_mm,1e-9),0,1);const drop=p.rim_shoulder_mm-p.parting_line_mm;return p.parting_line_mm+drop*Math.pow(1-u,1/2);}
function zBot(p,r){const{R,R_in}=derived(p);if(r<=R_in)return zTop(p,r)-p.plate_thick_mm;if(r>=R-p.nose_radius_mm){const t=clamp((r-(R-p.nose_radius_mm))/Math.max(p.nose_radius_mm,1e-9),0,1);return p.parting_line_mm*(1-Math.sqrt(Math.max(1-t*t,0)));}return 0;}
function massProperties(p,n){n=n||3000;const{R,R_in}=derived(p);const dr=R/(n-1);let V=0,Iz=0,Vz=0,rimV=0;for(let i=0;i<n;i++){const r=R*i/(n-1);const top=zTop(p,r),bot=zBot(p,r);const h=Math.max(top-bot,0);const zmid=.5*(top+bot);const ringV=2*Math.PI*r*h*dr;V+=ringV;Vz+=ringV*zmid;if(r>R_in)rimV+=ringV;const r_cm=r/10,h_cm=h/10;const dm=p.density_gcc*2*Math.PI*r_cm*h_cm*(dr/10);Iz+=dm*r_cm*r_cm;}const V_cm3=V/1000;return{volume_cm3:V_cm3,mass_g:V_cm3*p.density_gcc,Iz_gcm2:Iz,cg_height_mm:V>0?Vz/V:0,rim_mass_fraction:V>0?rimV/V:0,cavity_depth_mm:zTop(p,R_in)-p.plate_thick_mm,total_height_mm:zTop(p,0),inner_rim_dia_mm:2*R_in};}
const PRESETS={
  Putter:{D_out_mm:212,rim_width_mm:11.5,rim_shoulder_mm:16,parting_line_mm:9,dome_mm:5.5,plate_thick_mm:2.4,nose_radius_mm:2.0,density_gcc:.98},
  Midrange:{D_out_mm:212,rim_width_mm:13.5,rim_shoulder_mm:14.5,parting_line_mm:10,dome_mm:3.5,plate_thick_mm:2.2,nose_radius_mm:2.2,density_gcc:.97},
  Fairway:{D_out_mm:212,rim_width_mm:17.5,rim_shoulder_mm:13.5,parting_line_mm:9,dome_mm:2.0,plate_thick_mm:1.9,nose_radius_mm:1.8,density_gcc:.95},
  Distance:{D_out_mm:211,rim_width_mm:22,rim_shoulder_mm:13.5,parting_line_mm:8,dome_mm:1.0,plate_thick_mm:1.6,nose_radius_mm:1.7,density_gcc:.90},
};
const CLASS_EST={ // per-class estimates for what PDGA doesn't publish
  Putter:  {plate:3.2,nose:3.0},
  Midrange:{plate:2.6,nose:2.5},
  Fairway: {plate:2.0,nose:2.0},
  Distance:{plate:1.7,nose:1.7},
};
function paramsFromCertified(d){
  const est=CLASS_EST[d.cls];
  const shoulder=d.rd+est.plate;                       // cavity depth = shoulder − plate
  const dome=clamp(d.h-shoulder,0.2,8);                // total height = shoulder + dome
  const stab=d.fl?d.fl[3]+d.fl[2]:1;                   // fade+turn: net stability (neutral if unknown)
  const parting=clamp(shoulder*(0.60-0.04*stab),5,shoulder-2); // estimated: lower parting line → more overstable
  const p={D_out_mm:d.dia,rim_width_mm:d.rw,rim_shoulder_mm:+shoulder.toFixed(1),
    parting_line_mm:+parting.toFixed(1),dome_mm:+dome.toFixed(1),
    plate_thick_mm:est.plate,nose_radius_mm:est.nose,density_gcc:.95};
  // fit density so the model's mass lands on the certified max weight (clamped to real
  // plastics). PDGA lists max weight rounded up, so also respect the exact 8.3 g/cm rule.
  const target=Math.min(d.wt,8.3*d.dia/10)-0.5;
  const vol=massProperties({...p,density_gcc:1}).volume_cm3;
  p.density_gcc=+clamp(target/Math.max(vol,1e-9),.85,1.00).toFixed(3);
  return p;
}
