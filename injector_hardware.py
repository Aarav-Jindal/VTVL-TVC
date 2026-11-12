# Injector body (printed 316L) and pintle (machined 316) hardware: passages, ports, seals, walls and bolts.
# Checks the sizes in design_values.json and prints every dimension needed for the CAD.
# x is measured from the injector face (x 0; the chamber is +x, the injector body is -x), r from the engine axis.
# "depth" = distance from the face into the injector body. Clock = angle around the axis, 0 deg = top (Front Plane).
# Run chamber_hardware.py and injector.py first.
#
# Ethanol: fuel throttle valve -> inlet port (back face, clock 0) -> feed passage at clock 0 -> core manifold
#          -> 8 axial holes around the pintle. An internal fuel ring joins the 3 radial feed passages at their
#          outer ends, so all three feed the core manifold. The feed passage at clock 0 also feeds the film:
#          restrictor -> film manifold -> 8 film holes. The restrictor is drilled through the inlet port
#          before the fitting goes in (the inlet flow is only ~4 m/s, so it barely pushes on the restrictor).
# Nitrous: -6 ORB port in the pintle head -> pintle bore -> 16 radial holes.
import numpy as np
from CoolProp.CoolProp import PropsSI
import design as d

psiToPa = 6894.76
mm = 1e-3
Pc_Pa = d.get("Pc_psi") * psiToPa
P_design = d.get("P_tank_psi") * psiToPa # worst static case for walls: tank pressure in every passage, chamber empty (Pa)
mdot_fuel = d.get("mdot_fuel_kgps")
mdot_film = mdot_fuel * d.get("film_fraction")
mdot_core = mdot_fuel - mdot_film
mdot_ox = d.get("mdot_ox_kgps")
rho_ox = d.get("rho_ox_manifold") # nitrous inside the pintle is part vapor (injector.py)
P_manifold = Pc_Pa + d.get("dP_eth_injector_psi") * psiToPa # ethanol manifold pressure at full thrust (Pa)
rho_eth = PropsSI('D', 'T', d.get("T_coolant_out_K"), 'P', P_manifold, 'Ethanol')
Cd = d.get("injector_Cd")
allow_body = d.get("injector_yield_Pa") / d.get("injector_safety_factor")
allow_pintle = d.get("pintle_yield_Pa") / d.get("injector_safety_factor")

# Depth layers, face -> back: face plate, film manifold, restrictor wall, feed passages
t_face = d.get("face_plate_m")
h_fm = d.get("film_manifold_height_m")
L_restr = d.get("film_restrictor_len_m")
d_feed = d.get("feed_passage_d_m")
n_feed = int(d.get("n_feed_passages"))
depth_feed = t_face + h_fm + L_restr + d_feed / 2 # feed passage centre
depth_cavity = depth_feed + d_feed / 2 # collector groove and core manifold reach the feed passage bottom
T_hub = d.get("inj_hub_thickness_m")
r_hub = d.get("inj_hub_outer_r_m")
t_flange = d.get("inj_flange_thickness_m")

# Radii
R_p = d.get("D_pintle_mm") / 2 * mm # pintle radius = its bore in the body (sliding fit)
r_core_in = R_p + d.get("core_manifold_wall_m")
r_core_out = r_core_in + d.get("core_manifold_width_m")
d_eth = d.get("d_eth_axial_hole_mm") * mm
r_eth = r_core_in + 0.1*mm + d_eth / 2 # ethanol hole circle: hole edge 0.1 mm inside the manifold's inner wall
d_film = d.get("d_film_hole_mm") * mm
d_restr = d.get("d_film_restrictor_mm") * mm
r_film = d.get("R_film_ring_mm") * mm
r_fm_in = r_film - d.get("film_manifold_width_m") / 2
r_fm_out = r_film + d.get("film_manifold_width_m") / 2
r_ring = d.get("fuel_ring_d_m") / 2
w_ring = d.get("fuel_ring_width_m")
r_ring_in, r_ring_out = r_ring - w_ring / 2, r_ring + w_ring / 2 # closed ring from depth t_face to depth_cavity
r_inner_seal = d.get("inner_groove_id_m") / 2 + d.get("oring_groove_width_m") # chamber's inner O-ring groove OD
r_outer_seal_id = d.get("outer_groove_id_m") / 2

# Flow velocities: manifolds well below the jets they feed (same rule as the chamber manifolds)
v_eth_jet = d.get("v_eth_injector_mps")
v_film_jet = Cd * np.sqrt(2 * d.get("dP_film_holes_psi") * psiToPa / rho_eth)
v_feed = mdot_fuel / (rho_eth * n_feed * np.pi / 4 * d_feed**2)
v_ring = (mdot_fuel / n_feed) / (rho_eth * w_ring * (depth_cavity - t_face)) # inlet passage feeds the others both ways round
v_core = (mdot_core / n_feed / 2) / (rho_eth * (r_core_out - r_core_in) * (depth_cavity - t_face))
v_fm = (mdot_film / 2) / (rho_eth * (r_fm_out - r_fm_in) * h_fm) # one restrictor, flow splits both ways round

# Walls at the design pressure
def strip(P, w, t): # plate strip over a long narrow cavity, both edges fixed: max stress = P w^2 / (2 t^2)
    return P * w**2 / (2 * t**2)
s_face_core = strip(P_design, r_core_out - r_core_in, t_face)
s_face_ring = strip(P_design, w_ring, t_face)
s_face_film = strip(P_design, r_fm_out - r_fm_in, t_face)
s_restr_wall = strip(P_design, r_fm_out - r_fm_in, L_restr)
s_bore_wall = P_design * r_core_in / d.get("core_manifold_wall_m") # thin ring between core manifold and pintle bore
s_hub = 3 * (3 + 0.3) / 8 * P_design * (r_inner_seal / T_hub)**2 # hub as a simply supported plate inside the inner seal

# Flange: same ring bending as the chamber flange (same bolts and seals)
n_bolts = int(d.get("n_bolts"))
preload = 0.75 * d.get("bolt_proof_Pa") * d.get("bolt_stress_area_m2")
r_bc = d.get("bolt_circle_d_m") / 2
arm = r_bc - (r_outer_seal_id + d.get("oring_groove_width_m"))
m_ring = n_bolts * preload * arm / (2 * np.pi * r_bc)
s_flange = 6 * m_ring / t_flange**2
t_flange_needed = np.sqrt(6 * m_ring / allow_body)

# Pintle
r_pb = d.get("pintle_bore_d_m") / 2
t_pw = R_p - r_pb # pintle wall
d_ox = d.get("d_ox_radial_hole_mm") * mm
L_tip = d.get("L_pintle_tip_mm") * mm
x_bore_end = d.get("L_skip_mm") * mm + d_ox / 2 # bore stops at the top edge of the nitrous holes
t_end = L_tip - x_bore_end
t_end_needed = r_pb * np.sqrt(0.75 * P_design / allow_pintle) # clamped round plate: max stress = 3/4 P (a/t)^2
s_pintle_hoop = P_design * r_pb / t_pw
v_pintle_bore = mdot_ox / (rho_ox * np.pi * r_pb**2)
head_r = d.get("orb6_spotface_d_m") / 2
r_pbolt = d.get("pintle_bolt_r_m")
r_ear = r_pbolt + d.get("pintle_bolt_clearance_d_m") / 2 + 2*mm # ear end radius: 2 mm metal past the bolt hole
r_seal_od = d.get("pintle_seal_groove_od_m") / 2
F_pintle = P_design * np.pi * r_seal_od**2 # pressure inside the face seal pushes the pintle out
pintle_preload = 2 * 0.75 * d.get("bolt_proof_Pa") * d.get("pintle_bolt_stress_area_m2")
pintle_SF = pintle_preload / F_pintle

# Back-face ports: (radius, clock deg, passage dia). Film access sits over the restrictor at clock 0.
ports = {"Igniter": (d.get("igniter_port_r_m"), d.get("igniter_clock_deg"), d.get("igniter_passage_d_m")),
         "Pc tap": (d.get("pc_tap_r_m"), d.get("pc_tap_clock_deg"), d.get("pc_tap_passage_d_m")),
         "Ethanol inlet": (r_film, 0.0, 4*mm)}
def xy(r, clock): return np.array([r * np.sin(np.radians(clock)), r * np.cos(np.radians(clock))])
hex_r = d.get("orb4_hex_corner_r_m")
tap_r = d.get("orb4_tap_drill_d_m") / 2
feed_clocks = [i * 360 / n_feed for i in range(n_feed)]
ear_clocks = [90, 270]
def dist_to_ear(p, clock): # distance from point p to the ear rectangle at this clock
    u = np.array([np.sin(np.radians(clock)), np.cos(np.radians(clock))]) # along the ear
    w = np.array([u[1], -u[0]])
    a, b = p @ u, p @ w
    da = max(head_r - a, 0, a - r_ear); db = max(abs(b) - d.get("pintle_ear_width_m") / 2, 0)
    return np.hypot(da, db)
gaps = {} # metal (or clearance) between neighbouring features, m
gaps["Port bottom -> feed passages"] = (T_hub - d.get("orb4_tap_depth_m")) - depth_cavity
for name, (r, clock, dp) in ports.items():
    p = xy(r, clock)
    if name != "Ethanol inlet":
        gaps[f"{name} passage -> film manifold"] = r_fm_in - (r + dp / 2)
        gaps[f"{name} passage -> feed passage"] = min(np.linalg.norm(p - xy(r, c)) for c in feed_clocks) - dp / 2 - d_feed / 2
    gaps[f"{name} fitting -> pintle head"] = r - hex_r - head_r
    gaps[f"{name} fitting -> pintle ears"] = min(dist_to_ear(p, c) for c in ear_clocks) - hex_r
    gaps[f"{name} spotface inside hub"] = r_hub - (r + d.get("orb4_spotface_d_m") / 2)
    gaps[f"{name} port -> M4 holes"] = min(np.linalg.norm(p - xy(r_pbolt, c)) for c in ear_clocks) - tap_r - d.get("pintle_bolt_tap_drill_d_m") / 2
gaps["Pintle seal groove -> pintle bore"] = d.get("pintle_seal_groove_id_m") / 2 - R_p
gaps["Hub -> bolt heads"] = r_bc - d.get("bolt_head_d_m") / 2 - r_hub
names = list(ports)
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        (r1, c1, _), (r2, c2, _) = ports[names[i]], ports[names[j]]
        gaps[f"{names[i]} <-> {names[j]} fittings"] = np.linalg.norm(xy(r1, c1) - xy(r2, c2)) - 2 * hex_r

print("-- Ethanol Passages --")
print(f"Feed passages:                 {n_feed} x {d_feed*1e3:.1f} mm, {v_feed:.2f} m/s")
print(f"Fuel ring:                     {v_ring:.2f} m/s")
print(f"Core manifold:                 {v_core:.2f} m/s (jets {v_eth_jet:.0f} m/s, ratio {v_eth_jet/v_core:.0f}x)")
print(f"Film manifold:                 {v_fm:.2f} m/s (jets {v_film_jet:.1f} m/s, ratio {v_film_jet/v_fm:.1f}x)")
print(f"Hole L/D:                      ethanol {t_face/d_eth:.1f}, film {t_face/d_film:.1f}, restrictor {L_restr/d_restr:.1f}, nitrous {t_pw/d_ox:.1f}")

print(f"\n-- Walls at {P_design/psiToPa:.0f} psi (allow {allow_body/1e6:.0f} MPa body, {allow_pintle/1e6:.0f} MPa pintle) --")
print(f"Face plate over core manifold: {s_face_core/1e6:.1f} MPa")
print(f"Face plate over fuel ring:     {s_face_ring/1e6:.1f} MPa")
print(f"Face plate over film manifold: {s_face_film/1e6:.1f} MPa")
print(f"Restrictor wall:               {s_restr_wall/1e6:.1f} MPa")
print(f"Pintle bore wall (0.4 mm):     {s_bore_wall/1e6:.0f} MPa")
print(f"Hub plate:                     {s_hub/1e6:.0f} MPa")
print(f"Flange bending:                {s_flange/1e6:.0f} MPa, needs t >= {t_flange_needed*1e3:.1f} mm, have {t_flange*1e3:.1f} mm")
print(f"Pintle hoop:                   {s_pintle_hoop/1e6:.0f} MPa")
print(f"Pintle tip wall:               {t_end*1e3:.2f} mm (needs {t_end_needed*1e3:.2f})")
print(f"Pintle bolts:                  2 x M4, joint SF {pintle_SF:.1f} (need >= {d.get('joint_safety_factor'):.1f})")
print(f"Pintle nitrous bore:           {v_pintle_bore:.1f} m/s")

print("\n-- Geometry Checks (metal or clearance, mm) --")
for name, g in gaps.items():
    print(f"{name + ':':38s} {g*1e3:6.2f}")

problems = []
if min(v_eth_jet / v_core, v_film_jet / v_fm) < 3: problems.append("manifold too small for even flow")
if max(s_face_core, s_face_ring, s_face_film, s_restr_wall, s_bore_wall, s_hub) > allow_body: problems.append("body wall overstressed")
if t_flange < t_flange_needed: problems.append("injector flange too thin")
if s_pintle_hoop > allow_pintle or t_end < t_end_needed: problems.append("pintle wall too thin")
if pintle_SF < d.get("joint_safety_factor"): problems.append("pintle bolts too weak")
if max(t_face/d_eth, t_face/d_film) > 5.5: problems.append("small holes too long to drill")
if min(gaps.values()) < 1e-3: problems.append("less than 1 mm between features")
print("\nProblems:                      " + (", ".join(problems) if problems else "none"))

print("\n-- CAD Dimensions: Injector Body (mm, depth from the face) --")
print(f"Body:            hub r {r_hub*1e3:.1f} x {T_hub*1e3:.1f} deep, flange OD {d.get('flange_outer_d_m')*1e3:.0f} x {t_flange*1e3:.1f}, pintle bore dia {2*R_p*1e3:.2f} thru")
print(f"Bolt holes:      {n_bolts} x dia {d.get('bolt_clearance_d_m')*1e3:.1f} thru flange, on dia {2*r_bc*1e3:.0f}, clock 22.5 + 45n")
print(f"Fuel ring:       r {r_ring_in*1e3:.2f} to {r_ring_out*1e3:.2f}, depth {t_face*1e3:.2f} to {depth_cavity*1e3:.2f} (closed at the face)")
print(f"Core manifold:   r {r_core_in*1e3:.3f} to {r_core_out*1e3:.3f}, depth {t_face*1e3:.2f} to {depth_cavity*1e3:.2f}")
print(f"Ethanol holes:   8 x dia {d_eth*1e3:.3f} on r {r_eth*1e3:.3f}, clock 0 + 45n, face to manifold")
print(f"Film manifold:   r {r_fm_in*1e3:.3f} to {r_fm_out*1e3:.3f}, depth {t_face*1e3:.2f} to {(t_face+h_fm)*1e3:.2f}")
print(f"Film holes:      {int(d.get('n_film_holes'))} x dia {d_film*1e3:.3f} on r {r_film*1e3:.2f}, clock 22.5 + 45n (off the restrictor jet)")
print(f"Feed passages:   {n_feed} x dia {d_feed*1e3:.1f} radial at depth {depth_feed*1e3:.2f}, r {0.5*(r_core_in+r_core_out)*1e3:.2f} to {r_ring*1e3:.1f}, clock 0/120/240")
print(f"Restrictor:      dia {d_restr*1e3:.3f} at r {r_film*1e3:.2f} clock 0, depth {(t_face+h_fm)*1e3:.2f} to {(t_face+h_fm+L_restr)*1e3:.2f} (drill through the ethanol inlet port)")
for name, (r, clock, dp) in ports.items():
    print(f"{name + ':':16s} -4 ORB on back face at r {r*1e3:.2f} clock {clock:.0f}, passage dia {dp*1e3:.1f} to " + ("the feed passage" if name == "Ethanol inlet" else "the face"))
print(f"Pintle seal:     face groove r {d.get('pintle_seal_groove_id_m')/2*1e3:.1f} to {r_seal_od*1e3:.1f}, {d.get('pintle_seal_groove_depth_m')*1e3:.1f} deep, on the back face")
print(f"Pintle bolts:    2 x M4 tapped {d.get('pintle_bolt_tap_depth_m')*1e3:.0f} deep at r {r_pbolt*1e3:.1f}, clock 90/270")
print("\n-- CAD Dimensions: Pintle (mm, x from the face) --")
print(f"Tube:            dia {2*R_p*1e3:.3f} from x {-T_hub*1e3:.1f} to tip at x {L_tip*1e3:.3f}")
print(f"Nitrous bore:    dia {2*r_pb*1e3:.1f}, flat bottom at x {x_bore_end*1e3:.3f}")
n_ox = int(d.get("n_ox_radial_holes"))
print(f"Nitrous holes:   {n_ox} x dia {d_ox*1e3:.3f} radial at x {d.get('L_skip_mm'):.3f}, clock 0 + {360/n_ox:g}n (lined up with the ethanol holes)")
print(f"Head:            dia {2*head_r*1e3:.0f} from x {-T_hub*1e3:.1f} to {-(T_hub+d.get('pintle_head_len_m'))*1e3:.1f}, -6 ORB port in the top")
print(f"Ears:            {2*r_ear*1e3:.1f} long x {d.get('pintle_ear_width_m')*1e3:.0f} wide x {d.get('pintle_ear_thickness_m')*1e3:.0f} thick at clock 90/270, 2 x dia {d.get('pintle_bolt_clearance_d_m')*1e3:.1f} holes at r {r_pbolt*1e3:.1f}")

d.save_many({
    "inj_depth_cavity_mm": (depth_cavity*1e3, "mm"),
    "inj_depth_feed_mm": (depth_feed*1e3, "mm"),
    "core_manifold_r_in_mm": (r_core_in*1e3, "mm"),
    "core_manifold_r_out_mm": (r_core_out*1e3, "mm"),
    "r_eth_hole_circle_mm": (r_eth*1e3, "mm"),
    "film_manifold_r_in_mm": (r_fm_in*1e3, "mm"),
    "film_manifold_r_out_mm": (r_fm_out*1e3, "mm"),
    "pintle_bore_end_x_mm": (x_bore_end*1e3, "mm"),
    "pintle_ear_r_mm": (r_ear*1e3, "mm"),
    "inj_flange_stress_Pa": (s_flange, "Pa"),
    "pintle_joint_SF": (pintle_SF, "-"),
}, "injector_hardware.py")
