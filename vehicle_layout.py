# Vehicle layout: stacks the sections (top to bottom: avionics, helium bottle + regulator, ethanol tank, nitrous tank,
# engine bay) inside the body tube, places every part from parts_list.csv, and finds the height, center of gravity,
# gimbal lever arm and pitch inertia at liftoff (full) and landing (empty).
# Heights are measured up from the nozzle exit (x 0). Run after massIteration.py (the engine CG comes from the CAD).
import csv
import numpy as np
import design as d

g = 9.80665
mm = 1e-3
# ---- Sections, bottom to top ----
gap = d.get("section_gap_m"); cap = d.get("bulkhead_length_m") # tube length beyond the liquid + ullage, for both end caps
r_in = d.get("tank_outer_radius_m") - d.get("tank_wall_m"); A_in = np.pi*r_in**2
order = [("engine_bay", d.get("engine_bay_len_m")),
         ("n2o_tank", d.get("L_tank_N2O_m") + cap),
         ("eth_tank", d.get("L_tank_ethanol_m") + cap),
         ("he_valves", d.get("He_valve_stack_len_m")),
         ("he_bottle", d.get("He_bottle_len_m")),
         ("avionics", d.get("avionics_bay_len_m"))]
sec = {}; x = 0.0
for i, (name, length) in enumerate(order):
    if i and name != "he_bottle": x += d.get("n2o_eth_gap_m") if name == "eth_tank" else gap # the regulator stack sits right under the bottle; ethanol valves and filter sit between the tanks
    sec[name] = (x, x + length); x += length
L_stack = x # total height, nozzle exit to the top
x_bot = d.get("body_tube_bottom_m") # the tube starts just below the gimbal mount; the engine hangs out below it
L_body = L_stack - x_bot; x_body = 0.5*(L_stack + x_bot)
mid = lambda s: 0.5*(sec[s][0] + sec[s][1])

# ---- Where each kind of part sits (m above the nozzle exit) ----
L_gas = d.get("L_gas_total_mm")*mm # nozzle exit to injector face
x_face = L_gas
x_engine_top = x_face + d.get("inj_hub_thickness_m") + d.get("pintle_head_len_m") # top of the pintle head
x_pivot = x_face + d.get("gimbal_pivot_above_face_m")
spots = {"avionics": mid("avionics"), "he_valves": mid("he_valves"), "he_bottle": mid("he_bottle"),
         "eth_top": sec["eth_tank"][1] - cap/4, "eth_bottom": sec["eth_tank"][0] + cap/4,
         "n2o_top": sec["n2o_tank"][1] - cap/4, "n2o_bottom": sec["n2o_tank"][0] + cap/4,
         "engine_bay": 0.5*(x_engine_top + sec["engine_bay"][1]), "engine": x_face - d.get("x_engine_cg_mm")*mm,
         "engine_top": x_face + d.get("inj_hub_thickness_m") + 0.02, # throttle valves on the back face / pintle head
         "gimbal": x_pivot, "legs": -d.get("leg_cg_below_exit_m")}

# ---- Masses: (name, kg, x, own pitch inertia about its own center) ----
items = []
known = est = 0.0
with open("parts_list.csv", newline="", encoding="utf-8") as f:
    rows = [r for r in csv.DictReader(line for line in f if not line.startswith("#"))]
for r in rows:
    m = int(r["qty"])*float(r["each_g"])/1000
    if r["group"] == "known": known += m
    else: est += m
    w = r["where"]
    if w == "body":
        items.append((r["item"], m, x_body, m*L_body**2/12))
    elif w == "rings":
        ends = [sec["eth_tank"][1], sec["eth_tank"][0], sec["n2o_tank"][1], sec["n2o_tank"][0]]
        for e in ends: items.append((r["item"], m/len(ends), e, 0.0))
    else:
        items.append((r["item"], m, spots[w], 0.0))
R_body = d.get("body_od_m")/2
m_body = d.get("body_tube_kg_per_m")*L_body
items.append(("Body tube", m_body, x_body, m_body*(R_body**2/2 + L_body**2/12)))
m_engine = d.get("m_engine_cad_kg")
items.append(("Engine (CAD)", m_engine, spots["engine"], 0.0))
m_gimbal = d.get("m_gimbal_cad_kg")
items.append(("Gimbal (CAD)", m_gimbal, spots["gimbal"], 0.0))
for tank, key in (("n2o_tank", "m_tank_N2O_kg"), ("eth_tank", "m_tank_ethanol_kg")):
    lo, hi = sec[tank]; m = d.get(key)
    items.append((f"{tank} (tube + caps)", m, 0.5*(lo + hi), m*(r_in**2/2 + (hi - lo)**2/12)))
m_bottle = d.get("He_bottle_mass_kg")
items.append(("Helium bottle", m_bottle, mid("he_bottle"), m_bottle*d.get("He_bottle_len_m")**2/12))

# Propellant at liftoff: liquid from the bottom of each tank up; helium in the bottle. At landing: helium in the tanks.
def liquid(tank, m, rho):
    h = m/(rho*A_in); x0 = sec[tank][0] + cap/2
    return (tank + " liquid", m, x0 + h/2, m*(r_in**2/4 + h**2/12)), h
n2o, h_n2o = liquid("n2o_tank", d.get("m_N2O_kg"), d.get("rho_N2O"))
eth, h_eth = liquid("eth_tank", d.get("m_ethanol_kg"), d.get("rho_ethanol"))
m_He = d.get("m_He_kg")
V_n2o, V_eth = d.get("L_tank_N2O_m"), d.get("L_tank_ethanol_m") # same bore, so volume ~ length
he_full = [("Helium (bottle)", m_He, mid("he_bottle"), 0.0)]
he_empty = [("Helium (tanks)", m_He*V_n2o/(V_n2o + V_eth), mid("n2o_tank"), 0.0), ("Helium (tanks)", m_He*V_eth/(V_n2o + V_eth), mid("eth_tank"), 0.0)]

def props(lst):
    m = np.array([i[1] for i in lst]); xs = np.array([i[2] for i in lst]); own = np.array([i[3] for i in lst])
    M = m.sum(); xc = (m*xs).sum()/M
    return M, xc, (own + m*(xs - xc)**2).sum()
dry = items
M_full, x_full, I_full = props(dry + [n2o, eth] + he_full)
M_empty, x_empty, I_empty = props(dry + he_empty)

# ---- Checks against the mass budget ----
F = d.get("thrust_N")
budget_dry = known + m_engine + m_gimbal + est + d.get("m_body_tube_kg") + d.get("m_tank_N2O_kg") + d.get("m_tank_ethanol_kg") + m_bottle # massIteration.py, no margin
arm_full, arm_empty = x_full - x_pivot, x_empty - x_pivot

print("-- Stack (m above the nozzle exit) --")
for name, _ in order:
    lo, hi = sec[name]; print(f"{name:12s} {lo:6.3f} to {hi:6.3f}  ({hi - lo:.3f} m)")
print(f"Total height:                {L_stack:.2f} m")
print(f"Body tube:                   {x_bot:.3f} to {L_stack:.3f} m: {L_body:.2f} m long, {d.get('body_od_m')/0.0254:.1f} in OD, {m_body:.2f} kg")
print(f"Liquid heights at liftoff:   nitrous {h_n2o:.2f} m, ethanol {h_eth:.2f} m")

print("\n-- Mass check --")
print(f"Laid out (no margin):        {props(dry)[0]:.2f} kg dry; mass budget (massIteration.py, no margin) {budget_dry:.2f} kg")
print(f"Body tube:                   {L_body:.3f} m here, {d.get('L_body_tube_m'):.3f} m in the mass loop")
if abs(props(dry)[0] - budget_dry) > 0.01 or abs(L_body - d.get("L_body_tube_m")) > 1e-3:
    print("WARNING: the layout and the mass budget disagree (stacks or parts list out of step)")

print("\n-- Center of gravity and gimbal --")
print(f"Liftoff (full):              {M_full:.2f} kg (no margin), CG {x_full:.3f} m above the nozzle exit")
print(f"Landing (empty):             {M_empty:.2f} kg (no margin), CG {x_empty:.3f} m")
print(f"CG travel during the flight: {(x_empty - x_full)*1e3:+.0f} mm")
print(f"Gimbal pivot:                {x_pivot:.3f} m; lever arm to the CG {arm_full:.3f} m full, {arm_empty:.3f} m empty")
print(f"Pitch inertia about the CG:  {I_full:.2f} kg m^2 full, {I_empty:.2f} kg m^2 empty")
for label, arm, I in (("liftoff", arm_full, I_full), ("landing", arm_empty, I_empty)):
    tq = F*np.sin(np.radians(1))*arm
    print(f"Per 1 deg of gimbal ({label}): {tq:.2f} N m, {np.degrees(tq/I):.0f} deg/s^2")

print("\n-- Line lengths the layout implies (feed_system.py still uses placeholders) --")
print(f"Nitrous: tank bottom to the throttle valve ~{sec['n2o_tank'][0] - x_engine_top:.2f} m (placeholder ox_line_length_m {d.get('ox_line_length_m'):.2f})")
print(f"Ethanol: tank bottom, past the nitrous tank, to the jacket inlet ~{sec['eth_tank'][0]:.2f} m (placeholder fu_line_length_m {d.get('fu_line_length_m'):.2f})")

d.save_many({
    "L_vehicle_m": (L_stack, "m"),
    "x_body_mid_m": (x_body, "m"),
    "body_tube_mass_kg": (m_body, "kg"),
    "x_cg_liftoff_m": (x_full, "m"),
    "x_cg_landing_m": (x_empty, "m"),
    "x_gimbal_pivot_m": (x_pivot, "m"),
    "gimbal_arm_liftoff_m": (arm_full, "m"),
    "gimbal_arm_landing_m": (arm_empty, "m"),
    "I_pitch_liftoff_kgm2": (I_full, "kg*m^2"),
    "I_pitch_landing_kgm2": (I_empty, "kg*m^2"),
}, "vehicle_layout.py")
