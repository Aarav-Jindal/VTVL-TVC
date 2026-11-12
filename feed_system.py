# Feed system sizing: lines, filters, main valves, flex hoses, throttle valves, regulator, relief valve, tank end caps.
# Follows the P&ID (diagrams/pid.png):
#   nitrous: tank -> filter -> main valve -> flex hose -> throttle valve -> injector
#   ethanol: tank -> filter -> main valve -> flex hose -> cooling jacket -> engine loop -> throttle valve -> injector
# The ethanol throttle sits after the jacket, so the jacket stays near tank pressure at every thrust.
# Pressure budget from tank to chamber at full thrust and across the throttle range.
# Run after heat_transfer.py and injector.py (uses their jacket drop and injector hole sizes).
import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.optimize import brentq
import design as d
import two_phase as tp

psiToPa = 6894.76
m3sToGpm = 15850.3 # US gallons per minute
P_tank = d.get("P_tank_psi") * psiToPa
Pc = d.get("Pc_psi") * psiToPa
mdot_ox = d.get("mdot_ox_kgps")
mdot_fu = d.get("mdot_fuel_kgps") # all of it goes through the jacket (core + film)
T_ox = d.get("T_N2O_max_K")
T_fu = d.get("T_coolant_in_K")
rho_ox = PropsSI('D', 'T', T_ox, 'P', P_tank, 'NitrousOxide') # liquid nitrous, pushed above its boiling pressure (kg/m^3)
rho_fu = PropsSI('D', 'T', T_fu, 'P', P_tank, 'Ethanol')
mu_ox = 6.0e-5 # liquid nitrous viscosity near 25 C (Pa s); CoolProp has no viscosity model for it. Friction barely depends on it
mu_fu = PropsSI('V', 'T', T_fu, 'P', P_tank, 'Ethanol')
P_vap_ox = PropsSI('P', 'T', T_ox, 'Q', 0, 'NitrousOxide') # nitrous boils below this (Pa)
roughness = d.get("tube_roughness_m")
dP_inj_fu = d.get("dP_eth_injector_psi") * psiToPa
dP_jacket = d.get("jacket_dP_psi") * psiToPa
Cd = d.get("injector_Cd")
F_L = d.get("throttle_valve_FL")

# Pressure drop of each part at a flow rate (Pa)
def pipe_dP(mdot, rho, mu, D, L, K):
    v = mdot / (rho * np.pi / 4 * D**2)
    Re = rho * v * D / mu
    f = (-1.8 * np.log10((roughness / D / 3.7)**1.11 + 6.9 / Re))**-2 # Darcy friction factor (Haaland)
    return (f * L / D + K) * rho * v**2 / 2, v
def cv_dP(mdot, rho, Cv):
    # Liquid through a valve or filter: Q (gpm) = Cv * sqrt(dP (psi) / SG)
    Q_gpm = mdot / rho * m3sToGpm
    return (rho / 1000) * (Q_gpm / Cv)**2 * psiToPa
def cv_needed(mdot, rho, dP):
    return mdot / rho * m3sToGpm * np.sqrt((rho / 1000) / (dP / psiToPa))

# One propellant line: tube + filter + main valve + flex hose (everything before the throttle valve)
lines = {
    "ox": dict(mdot=mdot_ox, rho=rho_ox, mu=mu_ox, OD=d.get("ox_line_od_m"), wall=d.get("ox_line_wall_m"), L=d.get("ox_line_length_m"),
               K=d.get("ox_fittings_K"), hose_D=d.get("ox_hose_id_m"), filt=d.get("ox_filter_Cv"), valve=d.get("ox_main_valve_Cv")),
    "fu": dict(mdot=mdot_fu, rho=rho_fu, mu=mu_fu, OD=d.get("fu_line_od_m"), wall=d.get("fu_line_wall_m"), L=d.get("fu_line_length_m"),
               K=d.get("fu_fittings_K"), hose_D=d.get("fu_hose_id_m"), filt=d.get("fu_filter_Cv"), valve=d.get("fu_main_valve_Cv")),
}
def upstream_dP(p, f):
    """Losses from the tank to the throttle valve at throttle f (fraction of full flow): parts in Pa and total."""
    m = p["mdot"] * f
    tube, v = pipe_dP(m, p["rho"], p["mu"], p["OD"] - 2 * p["wall"], p["L"], p["K"])
    hose, v_hose = pipe_dP(m, p["rho"], p["mu"], p["hose_D"], d.get("hose_length_m"), d.get("hose_K"))
    parts = {"tube": tube, "filter": cv_dP(m, p["rho"], p["filt"]), "main valve": cv_dP(m, p["rho"], p["valve"]), "hose": hose}
    return parts, sum(parts.values()), v, v_hose

# Throttle valve sizing for a liquid that boils as it passes the valve seat (ISA/IEC liquid sizing):
# past dP_choke, boiling at the seat chokes the flow, so a bigger drop adds no flow (the valve sets the flow alone).
def valve_cv(mdot, fluid, P1, h1, P2):
    rho1 = PropsSI('D', 'P', P1, 'H', h1, fluid)
    Pv = tp.vapor_pressure(fluid, P1, h1)
    F_F = 0.96 - 0.28 * np.sqrt(Pv / PropsSI('pcrit', fluid)) # where the seat pressure ends up vs Pv
    dP_choke = F_L**2 * (P1 - F_F * Pv)
    return cv_needed(mdot, rho1, min(P1 - P2, dP_choke)), P1 - P2 >= dP_choke

# Injector manifold pressure: where the holes pass the flow. After the throttle valve the liquid may be part vapor
# (same enthalpy, lower pressure), which needs more pressure to push the same flow through the holes.
def manifold_pressure(fluid, mdot, area, P_in, h, P_chamber):
    flow = lambda P: Cd * area * tp.mass_flux(fluid, P, h, P_chamber) - mdot
    if flow(P_in) < 0: return np.nan # can't pass the flow even with the throttle valve wide open
    return brentq(flow, P_chamber + 1e3, P_in, xtol=50)

# Injector hole areas from injector.py (m^2). Ethanol: the core holes and the film path act together as one hole
# sized for all-liquid ethanol at full thrust (the film restrictor and holes in series are lumped in).
A_ox = d.get("A_ox_holes_mm2") * 1e-6
T_out_full = d.get("T_coolant_out_K")
rho_fu_hot = PropsSI('D', 'T', T_out_full, 'P', Pc + dP_inj_fu, 'Ethanol')
A_fu = mdot_fu / (Cd * np.sqrt(2 * rho_fu_hot * dP_inj_fu))
h_ox = PropsSI('H', 'T', T_ox, 'P', P_tank, 'NitrousOxide') # nitrous enthalpy, kept through the valve (J/kg)
# Engine loop: chamber outlet port -> tube over the flanges -> fuel throttle valve on the injector inlet port
loop = dict(D=lines["fu"]["OD"] - 2 * lines["fu"]["wall"], L=d.get("fu_engine_loop_length_m"), K=d.get("fu_engine_loop_K"))

# Throttle points: full thrust, hover at liftoff (full tanks), hover at landing (empty tanks), lowest throttle
g = 9.80665
F = d.get("thrust_N")
points = [("Full thrust", 1.0), ("Hover, liftoff", d.get("m_wet_kg") * g / F),
          ("Hover, landing", d.get("m_dry_kg") * g / F), ("Lowest throttle", d.get("throttle_min"))]
# Thrust ~ flow at fixed mixture ratio: Pc ~ f, jacket and line drops ~ f^2.
# Coolant temperature rise ~ heat / flow ~ f^0.8 / f (Bartz: heat flux ~ Pc^0.8)
P_crit_fu = PropsSI('pcrit', 'Ethanol')
rows = []
for name, f in points:
    _, up_ox, _, _ = upstream_dP(lines["ox"], f)
    _, up_fu, _, _ = upstream_dP(lines["fu"], f)
    Pc_f = Pc * f
    # Nitrous
    P_ox_tv_in = P_tank - up_ox # nitrous pressure reaching the throttle valve (Pa)
    P_ox_man = manifold_pressure('NitrousOxide', mdot_ox * f, A_ox, P_ox_tv_in, h_ox, Pc_f)
    Cv_ox, choke_ox = valve_cv(mdot_ox * f, 'NitrousOxide', P_ox_tv_in, h_ox, P_ox_man)
    # Ethanol: jacket first, then the engine loop, then the throttle valve
    P_jacket_out = P_tank - up_fu - dP_jacket * f**2 # (Pa)
    T_out = T_fu + (T_out_full - T_fu) * f**-0.2
    T_sat = PropsSI('T', 'P', P_jacket_out, 'Q', 0, 'Ethanol') if P_jacket_out < P_crit_fu else np.inf
    rho_hot = PropsSI('D', 'T', T_out, 'P', P_jacket_out, 'Ethanol')
    P_fu_tv_in = P_jacket_out - pipe_dP(mdot_fu * f, rho_hot, PropsSI('V', 'T', T_out, 'P', P_jacket_out, 'Ethanol'), loop["D"], loop["L"], loop["K"])[0]
    h_fu = PropsSI('H', 'T', T_out, 'P', P_fu_tv_in, 'Ethanol')
    P_fu_man = manifold_pressure('Ethanol', mdot_fu * f, A_fu, P_fu_tv_in, h_fu, Pc_f)
    Cv_fu, choke_fu = valve_cv(mdot_fu * f, 'Ethanol', P_fu_tv_in, h_fu, P_fu_man)
    rows.append(dict(name=name, f=f, ox_tv=P_ox_tv_in - P_ox_man, fu_tv=P_fu_tv_in - P_fu_man, Cv_ox=Cv_ox, Cv_fu=Cv_fu,
                     choke_ox=choke_ox, choke_fu=choke_fu, ox_margin=P_ox_tv_in - P_vap_ox,
                     P_ox_man=P_ox_man, x_ox=tp.state('NitrousOxide', P_ox_man, h_ox)[1],
                     P_fu_man=P_fu_man, x_fu=tp.state('Ethanol', P_fu_man, h_fu)[1],
                     P_fu_liquid=Pc_f + dP_inj_fu * f**2, Pv_fu=tp.vapor_pressure('Ethanol', P_fu_tv_in, h_fu),
                     P_jacket_out=P_jacket_out, T_out=T_out, T_sat=T_sat))

# Tube working pressure (Barlow, thin wall)
def tube_rating(p): return 2 * d.get("tube_uts_Pa") * p["wall"] / p["OD"] / d.get("tube_design_factor")
P_relief_flow = d.get("relief_set_psi") * d.get("relief_accumulation") * psiToPa # highest pressure the tanks and lines can see (Pa)

# Helium regulator: helium flow = propellant volume flow, at tank pressure and the (cold) tank gas temperature
Q_prop = mdot_ox / rho_ox + mdot_fu / rho_fu # m^3/s leaving the tanks
rho_He_tank = PropsSI('D', 'T', d.get("He_tank_gas_T_K"), 'P', P_tank, 'Helium')
mdot_He = rho_He_tank * Q_prop
rho_He_std = PropsSI('D', 'T', 288.71, 'P', 101325, 'Helium') # 60 F, 1 atm
scfm = lambda m: m / rho_He_std * 2118.88 # kg/s -> standard cubic feet per minute
SG_He = rho_He_std / PropsSI('D', 'T', 288.71, 'P', 101325, 'Air')
def gas_flow_scfm(Cv, P1_psia, P2_psia, T_K):
    # Gas through a valve (Swagelok / ISA form); flow stops rising once P2 < P1 / 2 (choked)
    dP = min(P1_psia - P2_psia, P1_psia / 2)
    return 22.67 * Cv * P1_psia * (1 - 2 * dP / (3 * P1_psia)) * np.sqrt(dP / (P1_psia * SG_He * T_K * 1.8))
P_tank_psi = d.get("P_tank_psi"); P_bottle_end = P_tank_psi + d.get("regulator_dropout_psi"); P_bottle_full = d.get("He_bottle_full_psi")
Cv_reg_needed = scfm(mdot_He) / gas_flow_scfm(1.0, P_bottle_end, P_tank_psi, 300) # worst case: lowest bottle pressure, warm gas
Cv_reg = d.get("regulator_Cv")
# Relief valve RV-1: must dump everything a stuck-open regulator lets through from a full bottle
Q_fail = gas_flow_scfm(Cv_reg, P_bottle_full, P_relief_flow / psiToPa, 300)
Cv_relief_needed = Q_fail / gas_flow_scfm(1.0, P_relief_flow / psiToPa, 14.7, 300)
P_proof = d.get("tank_proof_factor") * P_tank

# Tank end caps: flat plate held by radial screws through the tube wall
r_out = d.get("tank_outer_radius_m"); wall = d.get("tank_wall_m"); a = r_out - wall
sigma_t = d.get("tank_yield_Pa"); SF = d.get("tank_safety_factor"); nu = 0.33
F_cap = P_tank * np.pi * a**2 # end load on each cap (N)
t_cap = a * np.sqrt(3 * (3 + nu) / 8 * P_tank * SF / sigma_t) # simply supported flat plate (m)
m_cap = np.pi * a**2 * t_cap * 2700 * 1.3 # +30% for the seal land and port boss (kg)
n_b = int(d.get("bulkhead_n_bolts")); d_b = d.get("bulkhead_bolt_d_m"); e_b = d.get("bulkhead_bolt_edge_m")
A_shear = np.pi / 4 * (0.8 * d_b)**2 # roughly the minor diameter
SF_bolt_shear = n_b * A_shear * d.get("bulkhead_bolt_shear_Pa") / F_cap # vs ultimate
SF_bearing = n_b * d_b * wall * 1.5 * sigma_t / F_cap # tube wall crushing at the holes
SF_tearout = n_b * 2 * e_b * wall * 0.577 * sigma_t / F_cap # screw tearing out the tube end
ring_left = 1 - n_b * d_b / (2 * np.pi * r_out) # tube metal left at the screw row
SF_net = 2 * np.pi * r_out * wall * ring_left * sigma_t / F_cap

# ---- print ----
print("-- Propellants (full thrust) --")
print(f"Nitrous:   {mdot_ox*1e3:.1f} g/s, {rho_ox:.0f} kg/m^3 at {T_ox-273.15:.0f} C and {P_tank_psi:.0f} psi, boils below {P_vap_ox/psiToPa:.0f} psi")
print(f"Ethanol:   {mdot_fu*1e3:.1f} g/s, {rho_fu:.0f} kg/m^3")

print("\n-- Lines to the throttle valves (full thrust) --")
for key, label in (("ox", "Nitrous"), ("fu", "Ethanol")):
    p = lines[key]; parts, total, v, v_hose = upstream_dP(p, 1.0)
    print(f"{label}:   {p['OD']/0.0254:.3f} x {p['wall']/0.0254:.3f} in tube, {v:.2f} m/s (hose {v_hose:.2f} m/s), rated {tube_rating(p)/psiToPa:.0f} psi")
    print("           " + ", ".join(f"{k} {val/psiToPa:.1f}" for k, val in parts.items()) + f" -> {total/psiToPa:.1f} psi total")
print(f"Budget used in heat_transfer.py: {d.get('feed_losses_psi'):.0f} psi (lines 30 + main valve 20 + throttle 50)")

psi = lambda P: P / psiToPa
ch = lambda c: "choked" if c else "      "
print("\n-- Nitrous: throttle valve and injector (two-phase) --")
print(f"{'':16s}{'thrust':>7s} {'before valve':>13s} {'valve drop':>11s} {'Cv':>7s} {'':6s} {'injector in':>12s} {'vapor':>6s}")
for r in rows:
    print(f"{r['name']:16s}{r['f']*100:6.0f}% {psi(r['ox_margin']):5.0f} psi sub {psi(r['ox_tv']):7.0f} psi {r['Cv_ox']:7.4f} {ch(r['choke_ox'])} "
          f"{psi(r['P_ox_man']):8.0f} psi {r['x_ox']*100:5.0f}%")
print("('before valve' = margin above the nitrous boiling pressure; 'choked' = boiling at the seat sets the flow, chamber pressure can't push back)")
print("\n-- Ethanol: jacket, throttle valve and injector --")
print(f"{'':16s}{'thrust':>7s} {'jacket out':>11s} {'coolant/boils':>14s} {'valve drop':>11s} {'Cv':>7s} {'':6s} {'injector in':>12s} {'boils below':>12s}  boils where")
for r in rows:
    boils = "no (>crit)" if np.isinf(r['T_sat']) else f"{r['T_sat']:4.0f} K"
    where = f"after the valve ({r['x_fu']*100:.1f}% vapor by mass)" if r['x_fu'] > 0.001 else ("in the injector holes" if r['Pv_fu'] > Pc * r['f'] else "nowhere")
    print(f"{r['name']:16s}{r['f']*100:6.0f}% {psi(r['P_jacket_out']):7.0f} psi {r['T_out']:5.0f} K/{boils:>8s} {psi(r['fu_tv']):7.0f} psi {r['Cv_fu']:7.4f} {ch(r['choke_fu'])} "
          f"{psi(r['P_fu_man']):8.0f} psi {psi(r['Pv_fu']):8.0f} psi  {where}")
print("(injector in = manifold pressure the holes need; boiling in the holes chokes them, so it rises above the all-liquid Pc + drop: "
      + ", ".join(f"{psi(r['P_fu_liquid']):.0f}" for r in rows) + " psi)")
Cv_ox = [r["Cv_ox"] for r in rows]; Cv_fu = [r["Cv_fu"] for r in rows]
print(f"Throttle valve range needed: nitrous Cv {min(Cv_ox):.4f} to {max(Cv_ox):.4f}, ethanol Cv {min(Cv_fu):.4f} to {max(Cv_fu):.4f}")
boiling = [r for r in rows if r["T_out"] > r["T_sat"] - 10]
flashing = [r for r in rows if r["x_fu"] > 0.01] # more than a trace of vapor reaching the injector

print("\n-- Helium regulator and relief --")
print(f"Helium flow at full thrust:    {mdot_He*1e3:.2f} g/s = {scfm(mdot_He):.1f} SCFM (tank gas at {d.get('He_tank_gas_T_K'):.0f} K)")
print(f"Regulator Cv needed:           {Cv_reg_needed:.4f} at the end of the burn ({P_bottle_end:.0f} -> {P_tank_psi:.0f} psi), have {Cv_reg:.3f} ({Cv_reg/Cv_reg_needed:.1f}x)")
print(f"Regulator stuck open:          {Q_fail:.0f} SCFM from a full bottle")
print(f"Relief RV-1 Cv needed:         {Cv_relief_needed:.3f} (set {d.get('relief_set_psi'):.0f} psi, full flow at {P_relief_flow/psiToPa:.0f} psi)")
print(f"Tank proof pressure:           {P_proof/psiToPa:.0f} psi (relief full-flow pressure must stay below it)")

print("\n-- Tank end caps (each) --")
print(f"End load:                      {F_cap/1e3:.1f} kN at {P_tank_psi:.0f} psi")
print(f"Flat cap thickness:            {t_cap*1e3:.1f} mm (6061-T6, SF {SF:.0f}), about {m_cap:.2f} kg (budget {d.get('bulkhead_mass_kg'):.2f} kg)")
print(f"Radial screws:                 {n_b} x {d_b*1e3:.2f} mm: shear SF {SF_bolt_shear:.1f}, tube bearing SF {SF_bearing:.1f}, "
      f"tear-out SF {SF_tearout:.1f}, tube at screw row SF {SF_net:.1f}")

problems = []
for key, label in (("ox", "nitrous"), ("fu", "ethanol")):
    p = lines[key]; _, total, v, v_hose = upstream_dP(p, 1.0)
    v_limit = d.get("ox_line_velocity_max_mps") if key == "ox" else d.get("line_velocity_max_mps") # nitrous has its own, higher speed limit
    if max(v, v_hose) > v_limit: problems.append(f"{label} line too fast")
    if tube_rating(p) < P_relief_flow: problems.append(f"{label} tube rated below the relief pressure")
if any(np.isnan(r["P_ox_man"]) or np.isnan(r["P_fu_man"]) for r in rows): problems.append("injector holes can't pass the flow even with a throttle valve wide open")
if min(r["ox_margin"] for r in rows) < d.get("N2O_subcool_margin_psi") * psiToPa: problems.append("nitrous could boil before the throttle valve")
if boiling: problems.append("ethanol boils (or within 10 K) in the jacket at: " + ", ".join(r["name"].lower() for r in boiling))
if flashing: problems.append("hot ethanol partly boils after its throttle valve (two-phase injection) at: " + ", ".join(r["name"].lower() for r in flashing))
if Cv_reg < 1.5 * Cv_reg_needed: problems.append("regulator too small")
if P_relief_flow > P_proof: problems.append("relief pressure above tank proof")
if min(SF_bolt_shear, SF_bearing, SF_tearout, SF_net) < SF: problems.append("end-cap screws too weak")
print("\nProblems:                      " + ("; ".join(problems) if problems else "none"))

d.save_many({
    "ox_feed_dP_psi": (upstream_dP(lines["ox"], 1.0)[1] / psiToPa, "psi"),
    "fu_feed_dP_psi": (upstream_dP(lines["fu"], 1.0)[1] / psiToPa, "psi"),
    "ox_throttle_Cv_max": (max(Cv_ox), "-"), "ox_throttle_Cv_min": (min(Cv_ox), "-"),
    "fu_throttle_Cv_max": (max(Cv_fu), "-"), "fu_throttle_Cv_min": (min(Cv_fu), "-"),
    "hover_throttle_liftoff": (points[1][1], "-"), "hover_throttle_landing": (points[2][1], "-"),
    "P_jacket_out_hover_landing_psi": (rows[2]["P_jacket_out"] / psiToPa, "psi"),
    "x_fu_injector_hover_landing": (rows[2]["x_fu"], "-"),
    "x_ox_injector_full": (rows[0]["x_ox"], "-"),
    "He_flow_scfm": (scfm(mdot_He), "SCFM"),
    "regulator_Cv_needed": (Cv_reg_needed, "-"),
    "relief_Cv_needed": (Cv_relief_needed, "-"),
    "tank_cap_thickness_mm": (t_cap * 1e3, "mm"),
    "tank_cap_mass_kg": (m_cap, "kg"),
    "tank_cap_screw_SF": (min(SF_bolt_shear, SF_bearing, SF_tearout, SF_net), "-"),
}, "feed_system.py")
